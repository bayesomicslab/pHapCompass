#!/usr/bin/env python3
"""Takes an indexed reference genome, SNP frequencies, Indel frequencies, or alternatively SNP/Indel 
positions as input and produces modified haplotypes for the specified ploidy level, by making variations 
to the given reference genome according to the chosen stochastic model. Output could be used to simulate NGS reads 
with desired characteristics. Alternatively in random mode, the software could generate random fasta files 
for the specified ploidy, with the specified contig lengths and contig content. The (polyploid) variant 
haplotypes, i.e. the ones containing just the indel/mutations, are also coded numerically and written 
columnwise in a single file at the mutation mode, along with the position of each variant along its contig 
and the name of its contig, and the reference/alternative alleles.
Written by Ehsan Motazedi, 19-10-2015, Wageningen UR
Last updated 27-03-2018"""
import sys
import ast
import copy
import getopt
import numpy as np
import os 
import os.path
import random
import re 
import tempfile
import traceback
from Bio import SeqIO
# from Bio.Alphabet import *  # Deprecated in Biopython >=1.78
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from bisect import bisect
from collections import deque
from collections import OrderedDict
from math import ceil
from math import exp
from math import sqrt
from mugen_3 import *
from scipy.stats import binom
from shutil import move
from textwrap import wrap
from typing import List, Tuple, Optional, Dict, Any, Sequence
from collections import OrderedDict
from bisect import bisect


class CustomException(Exception):
    def __init__(self, value):
        self.parameter = value
    def __str__(self):
        return repr(self.parameter)

def get_empirical_dist(emp_filename=''):
    """Extract the empirical distribution from the given file"""
    sample=deque()
    lengths=deque()
    _line_num=0
    with open(emp_filename,'r') as empirical:
        for _n, _line in enumerate(empirical):
            lengths.append(len(_line.strip().split()))
            _line_num=_n+1
    if not _line_num:
        raise CustomException("ERROR: the empirical distribution file was empty!")
    with open(emp_filename,'r') as empirical:    
        if _line_num>1:
            while lengths:
                _x=lengths.pop()
                if _x>1:
                    raise CustomException("ERROR: the empirical distribution file must contain just one line or just one column of positive integer values!")
            for _line in empirical:
                sample.append(_line.strip())
                sample=list(sample)
        else:    
            sample=empirical.readlines().strip().split()
    try:
        _test=[float(_val)==int(float(_val)) for _val in sample]
    except ValueError as err:
        raise ValueError('The empirical distribution file must only contain positive integers for the SNP distances!')
    if False in _test:
        raise ValueError('The empirical distribution file must only contain positive integers for the SNP distances!')
    _test=[float(_val)>0 for _val in sample]
    if False in _test:
        raise ValueError('The empirical distribution file must only contain positive integers for the SNP distances!')
    return([int(float(_val)) for _val in sample])

def getinput(argv):
    """Parser for the command line input arguments""" 
    emp_filename = ""
    fastafile = ""              
    outputfile = ""      				# Other parameters default to None 
    muprob=dict()       				# Deault is no change compared to the reference file
    insertprob=dict()
    deleteprob=dict()
    rates=list()
    sdlogs=list()
    mupos=list()
    inpos=list()
    delpos=list()
    got_model=False
    got_file=False
    mualphabet=dict()
    verbose=False 						# Deafult to no verbose
    _random=False 						# Mutation is the default mode of operation
    ploidy=1						# Default to produce just one haplotype  
    length=list()
    alphabet=list()
    content=dict()
    names=list()
    wrapsize=60
    genomic=False
    model=3
    dosage=None
    seed=None
    region_start=None
    region_end=None
    try:
        opts, args = getopt.getopt(argv,"hf:s:i:d:S:I:D:m:o:vrp:l:a:c:n:w:g",["help","reference=","subfreq=",
            "insfreq=","delfreq=","subpos=","inspos=","delpos=","mutations=","output=","verbose","random",
            "ploidy=","length=","alphabet=","content=","names=","wrapsize=","genomic","model=","dosage=","sdlog=","empirical=","rndseed=","region-start=","region-end="])
        hlpmsg='\nUsage:\n\n\
Mutation mode (default): haplogenerator.py [options] <-f reference file> <-o output base name> [model argument]\n\
Generates homologues by manipulating the given reference sequences. Writes the ouput to output_hap1.fa,...\n\
output_hap[k].fa where [k] is the specified ploidy. Variant haplotypes saved in output_varianthaplos.txt\n\n\
Main arguments:\n\
-f, --reference STR  the indexed reference fasta file\n\
-o, --output    STR  the base name of the output file(s) containing haplotypes\n\n\
Options:\n\
-s, --subfreq   DICT substitution probability for each letter in the sequence alphabet used in the scan-catalyze model\n\
                LIST mutation, insertion and deletion rates, respectively, to be used in Poisson and Gaussian models\n\
                LIST mean of the log-distance between successive mutations, insertios and deletions to be used by the lognormal model.\n\
                LIST probabilities for a varinat to be a mutation, or an insertion, or a deletion, used by the empirical model\n\
-i, --insfreq   DICT insertion probability for each letter in the sequence alphabet used in scan catalyze model\n\
-d, --delfreq   DICT deletion probability for each letter in the sequence alphabet used in scan catalyze model\n\
-S, --subpos    LIST position of substituions over the genome\n\
-I, --inspos    LIST position of insertions over the genome\n\
-D, --delpos    LIST position of deletions over the genome\n\
-m, --mutations DICT possible substituations for each letter of the alphabet\n\
-v, --verbose        print out warning and status messages\n\
-p, --ploidy    INT  the number of homologues to be produced\n\
-g, --genomic        genomic coordinates will be used for indel/mutation positions instead of the default fasta coordinates\n\
--model  CONST poisson   apply a Poisson process for modelling indel/mutations, gets dosage distribution as parameter\n\
               distance  apply an independent Guassian random process for modelling indel/mutaions, gets dosage distribution as parameter\n\
               scancat   apply a scan catalyze approach, i.e. change each position independently with the given probability, for modelling indel/mutations (default), acts independently on each homologue\n\
               lognormal apply a lognormal process for modelling indel/mutations, gets dosage distribution and sdlog as parameter\n\
--empirical     STR  the name of the text file contaning numerical values corresponding to the distances between successive variants,\n\
                     in a single row or a single column, from which the SNP poisition will be sampled with replacement. NOT compatible\n\
                     with --model option!\n\
--dosage LIST  percentage of simplex, duplex, ..., ploidy-plex variants among all the variants, used only with poisson, lognormal and distance models.\n\
--sdlog  LIST  standard deviation on the log scale, for the distance between neighbour mutation, insertion and deletion points. Only for lognormal model\n\
--rndseed FLOAT   set the seed for the generation of random variation sites\n\n\
Random mode: haplogenerator.py -r [options] <-o output base name> <-a sequence alphabet> <-l length of the desired contigs within each homologue>\n\
Generates random homologues with the specified lengths. Sequence letters will be randomly chosen from the specified alphabet.\n\n\
Main arguments:\n\
-o, --output    STR  the base name of the output file(s) containing haplotypes\n\
-r, --random         operate in the random mode\n\
-a, --alphabet  LIST the set of the letters forming the sequences\n\
-l, --length    LIST list of values specifying the length of each contig sequence to be produced\n\n\
Options:\n\
-v, --verbose        print out warning and status messages\n\
-c, --content   DICT the frequency of each alphabet letter within the contigs\n\
-n, --names     LIST the names of the contigs to be written, should have equal elements as length\n\
-w, --wrapsize	LIST list of integers specifying the line length of the fasta file for each contig,\
 should have equal elements as length\n\
-p, --ploidy    INT  the number of homologues to be produced\n\n\
Application notes:\n\n\
1) Ploidy level is default to 1.\n\
2) Double quotes should always be put around the lists and dictionaries specified on the command line. Inside lists and dictionaries, alphabet letters should be single quoted, real and integer numbers should NOT be quoted.\n\
3) Contig coordinates (for specifying indel/mutation positions) start from the offsets in the fasta file, unless <-g, --genomic> option is set in which\n\
   case genomic coordinates will be considered.\n\
4) In mutation mode, mutaion probabilities are deafult to zero. Mutation positions default to none. If no mutation alphabet is specified, a purely random\n\
   mutation would happen.\n\
5) Multiletter values in the mutation alphabet dictionary are allowed. The mutation of the corresponding key will be each time randomly to one of those letters.\n\
6) In random mode, each letter will have equal frequency if no content is given. If no names are specified, contig names will be set to Chr1, Chr2,... . If no\n\
   wrap size specified for a contig, deafult of 60 will be used for wrapping the lines in the fasta file.\n\
7) Fasta coordinates for each contig start from its corresponding offset in the fasta index! Genomic coordinates start from 1 at the beginning of the genome\n\
   and increase continuously!\n\
8) The default model for indel/mutation is scan catalyze model. The input sequence is scanned through a hypothetical polymerase, and an indel\n\
   or mutation occurs independently at each point with probabilities 0.25,0.25 and 0.5, respectively. This catalyst error is corrected with probability\n\
   1-P(letter) at each point independently, with the probabilities being specified as input dictionaries. Alternative models are: 1) Poisson, in which\n\
   indel/mutation positions are randomly selected according to a Poisson model. Poisson process has a dynamic random rate at each event, chosen from\n\
   (rate/4, rate, rate/4+3/4) values with weights (1/6,2/3,1/6) where the "rates" are the given indel/mutation rates. 2) Distance model, in which equally\n\
   distant indel/mutation points are selected with random start and end positions, so that the mean distance is equal to 1/rate (for different start and\n\
   end points). The points are then shifted with a Gaussian distribution of mean zero and variance equal to 0.3/rate. 3) The lognormal\n\
   model, which acts similar to Poisson model but uses a lognormal distribution instead of Poisson and has a fixed rate and variance. With this model,\n\
   the given rates list should specify the mean log distances between neighbour mutation/indel points, and sdlog could be specified to determine the\n\
   sd of the log distance. Zero rates are considered as "NO mutation/indel" with these models.\n\
9) Dosage percentages must be given in a list, for simplex, duplex, ... up to the dosage equal to the ploidy level, respectively.'
        if not (opts or args):
            opts=[('-h','')]
        elif args:
            raise getopt.GetoptError('ERROR: Invalid parameter specification!' , None )
        else:
            pass
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                print(hlpmsg)
                sys.exit(0)
            elif opt in ("-f", "--reference"):
                fastafile = arg
            elif opt in ("-o", "--output"):
                outputfile=arg
                # make sure if output directory is given, it exists
                outdir=os.path.dirname(outputfile)
                if outdir and not os.path.exists(outdir):
                    # make directory
                    os.makedirs(outdir, exist_ok=True)
            elif opt in ("-r", "--random"):
                _random=True	 
            elif opt in ("-v", "--verbose"):
                verbose=True
            elif opt in ("-g", "--genomic"):
                genomic = True
            elif opt in ("-p", "--ploidy"):
                try:
                    ploidy=int(arg)
                except TypeError:
                    print("ERROR: ploidy level must be a positive integer!")
                    raise 
                except ValueError:
                    print("ERROR: ploidy level must be a positive integer!")
                    raise
            elif opt=="--model":
                got_model=True
                if str(arg)=='poisson':
                    model=1
                elif str(arg)=='distance':
                    model=2
                elif str(arg)=='scancat':
                    model=3
                elif str(arg)=='lognormal':
                    model=4
                else:
                    raise CustomException("ERROR: the specified model is not recognized! Only, poisson, lognormal, distance and scancat are allowed!")
            elif opt=="--empirical":
                got_file=True
                emp_filename = arg
                model=-1
            elif opt=="--rndseed":
                try:
                    seed=float(arg)
                except TypeError:
                    print("ERROR: the random seed must be a real number!")
                    raise 
                except ValueError:
                    print("ERROR: the random seed must be a real number!")
                    raise
            elif opt=="--region-start":
                try:
                    region_start=int(arg)
                except (TypeError, ValueError):
                    print("ERROR: region-start must be a positive integer!")
                    raise
            elif opt=="--region-end":
                try:
                    region_end=int(arg)
                except (TypeError, ValueError):
                    print("ERROR: region-end must be a positive integer!")
                    raise
            else:
                pass
        for opt, arg in opts:
                try:
                    if opt in ("-s", "--subfreq"):
                        if model==3:
                            muprob = ast.literal_eval(arg)						# A dictionary is expected for scan catalyze model
                        else:
                            rates = ast.literal_eval(arg)						# A list of rates is expected for Poisson and Gaussian models
                    elif opt in ("-i", "--insfreq") and model==3:
                        insertprob = ast.literal_eval(arg)
                    elif opt in ("-d", "--delfreq") and model==3:
                        deleteprob = ast.literal_eval(arg)
                    elif opt in ("-S", "--subpos"): 							# SNPpos to be implemented later
                        mupos = ast.literal_eval(arg)
                    elif opt in ("-I", "--inspos"): 							# Indelpos to be implemented later
                        inpos = ast.literal_eval(arg)
                    elif opt in ("-D", "--delpos"):
                        delpos = ast.literal_eval(arg)
                    elif opt in ("-m", "--mutations"):
                        mualphabet = ast.literal_eval(arg)
                    elif opt in ("-l", "--length"):
                        length = ast.literal_eval(arg)
                    elif opt in ("-a", "--alphabet"):
                        alphabet = ast.literal_eval(arg)
                    elif opt in ("-c", "--content"):
                        content = ast.literal_eval(arg)
                    elif opt in ("-n", "--names"):
                        names = ast.literal_eval(arg)
                    elif opt in ("-w", "--wrapsize"):
                        wrapsize = ast.literal_eval(arg)
                    elif opt=="--dosage":
                        dosage = ast.literal_eval(arg)
                    elif opt=="--sdlog":
                        sdlogs = ast.literal_eval(arg)
                    else:
                        pass
                except ValueError as err:
                    if not err.args: 
                        err.args=('',)
                    err.args =err.args + ("Input parameter values in incorrect format!\n\
Please check the format, values and quotations! Type -h, --help for help!",)			
                    raise
                except SyntaxError:
                    print('\n'.join(["Input parameter values in incorrect format!",
                        "Please check the format, values and quotations!","Type -h, --help for help!"]))
                    raise			
        if not outputfile:
            raise CustomException("ERROR: Must specify the base name of the output file(s)!")
        elif ploidy<=0:
            raise ValueError("Ploidy level must be a positive integer!")  
        elif not _random:
            if not fastafile:
                raise CustomException("ERROR: You must specify a fatsa file for the reference!")
            elif not os.path.isfile(fastafile):
                raise CustomException("ERROR: Could not find the fasta file!")
            elif not os.path.isfile(fastafile+".fai"):
                raise CustomException("ERROR: Could not find the index of the fasta file! \
The index must have the name fastafile.fai!")
            elif got_model and got_file:
                raise ValueError("--model and --empirical options are incompatible with each other!")
            elif got_file:
                if not emp_filename:
                    raise ValueError("You must specify the file name containing the empirical values for the successive variants distance!")
                elif not os.path.isfile(emp_filename):
                    raise ValueError("No file exists with the file name specified for the empirical varinats distances!")
            elif not isinstance(mualphabet,dict) and mualphabet:
                raise ValueError("-m,--mutations needs a quoted dictionary!")
            elif not isinstance(muprob,dict) and muprob:
                raise ValueError("-s, --subfreq needs a quoted dictionary for the specified model!")
            elif not isinstance(rates,list) and rates:
                raise ValueError("-s, --subfreq needs a quoted list for the specified model!")
            elif not isinstance(dosage,list) and dosage:
                raise ValueError("--dosage needs a quoted list for the specified model!")
            elif not isinstance(sdlogs,list) and sdlogs:
                raise ValueError("--sdlog needs a quoted list for the specified model!")
            elif not isinstance(insertprob,dict) and insertprob:
                raise ValueError("-i, --insfreq needs a quoted dictionary!")
            elif not isinstance(deleteprob,dict) and deleteprob:
                raise ValueError("-d, --delfreq needs a quoted dictionary!")
            elif not isinstance(mupos,list) and mupos:
                raise ValueError("-S, --subpos needs a quoted list!")
            elif not isinstance(inpos,list) and inpos:
                raise ValueError("-I, --inspos needs a quoted list!")
            elif not isinstance(delpos,list) and delpos:
                raise ValueError("-D, --delpos needs a quoted list!")
            else:
                pass														   # Other ecxeptions handeled in the main module
        elif fastafile:
            raise CustomException("ERROR: fasta file should not be specified with the random mode!")
        elif not alphabet:
            raise CustomException("ERROR: you must specify an alphabet for random sequence generation!")
        elif not length:
            raise CustomException("ERROR: contig lengths must be specified for random sequence generation!")
        elif not isinstance(alphabet,list):
            raise ValueError("-a, --alphabet needs a quoted list!")
        elif not isinstance(length,list):
            raise ValueError("-l, --length needs a quoted list!")
        elif not isinstance(content,dict) and content:
            raise ValueError("-c, --content needs a quoted dictionary!")
        else: 
            if names:
                if not isinstance(names,list):
                    raise ValueError("-n, --names needs a quoted list of strings!")
                elif len(names)!=len(length):
                    raise CustomException("ERROR: contig names and number of contigs must match eachother!")
                else:
                    pass
            if not isinstance(wrapsize,int):
                if not (isinstance(wrapsize,list)):
                    raise ValueError("-w, --wrapsize needs a quoted list!")
                elif len(wrapsize)!=len(length):
                    raise CustomException("ERROR: wraping size and number of fasta records (contigs) must match eachother!")
                else:
                    pass
            if muprob or insertprob or deleteprob or mupos or inpos or delpos or mualphabet:
                if verbose:
                    print('WARNING: you have specified Mutation mode options for the Random mode! These will be ignored!', file=sys.stderr)
            else:
                pass																# Other exceptions dealt with in th main module
        return (fastafile,outputfile,muprob,insertprob,deleteprob,mupos, 			# MuGen class used to procude haplotypes
        inpos,delpos,mualphabet,verbose,_random,
        length,alphabet,content,ploidy,names,wrapsize,genomic,
        model,rates, dosage, sdlogs, emp_filename, seed, region_start, region_end)
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
    except CustomException as instance:                                             
        print(instance) 
        sys.exit(2)
    except ValueError as e:
        print("ERROR: "+str(type(e))+ '\n' + '\n'.join(e.args))
        sys.exit(2) 
    except SyntaxError:
        raise
        sys.exit(2)
    except SystemExit as e:
        if e.code==0:								# No error has caused the exit
            os._exit(0)
        else:
            raise
    except:
        print("Unexpected Error: %s" %sys.exc_info()[0])
        raise
        sys.exit(2) 

def print_info():					
    """Procedure to print information about the input values to stderr"""
    try:
        if _random:
            print("Random mode parameters:", file=sys.stderr)
            print("output base name=", outputbase, file=sys.stderr)
            print("Sequence alphabet=", randalphabet, file=sys.stderr)
            print("Contig lengths=", length, file=sys.stderr)
            print("Contig content=", content, file=sys.stderr)
            print("Contig names=", contig_names, file=sys.stderr)
            print("Fasta record line lenghts=", linelens, file=sys.stderr)
        else:
            print("Mutation mode parameters:", file=sys.stderr)
            print("reference=", fastafile, file=sys.stderr)
            print("output base name=", outputbase, file=sys.stderr)
            print("mupos=", muposlst, file=sys.stderr)
            print("inpos=", inposlst, file=sys.stderr)
            print("delpos=", delposlst, file=sys.stderr)
            print("mualphabet=", mualphabet, file=sys.stderr)
            print("dosage percentages=", dosage, file=sys.stderr)
            if emp_filename:
                print("Distances between successive variants will be sampled from the values in {}".format(emp_filename), file=sys.stderr)  
            elif model==3:
                print("Sequence evolution model= Scan-Catalyze:", file=sys.stderr)
                print("	muprob=", muprob, file=sys.stderr)
                print("	insertprob=", insertprob, file=sys.stderr)
                print("	deleteprob=", deleteprob, file=sys.stderr)
                if dosage:
                    print("WARNING: specified dosages will be ignored for scancat model!", file=sys.stderr)
            elif model==2:
                print("Sequence evolution model= Gaussian", file=sys.stderr)
                print("	mutation/insertion/deletion rates=", rates, file=sys.stderr)
                print("WARNING: if specified, scan-cat indel/mutation probabilites will be ignored!", file=sys.stderr)
            elif model==4:
                print("Sequence evolution model= lognormal:", file=sys.stderr)
                print("	mean of mutation/insertion/deletion log distance=", rates, file=sys.stderr)
                print("	standard deviation of mutation/insertion/deletion log distance=", sdlog, file=sys.stderr)
                print("WARNING: if specified, scan-cat indel/mutation probabilites will be ignored!", file=sys.stderr)
            else:
                print("Sequence evolution model= poisson:", file=sys.stderr)
                print("	mutation/insertion/deletion rates=", rates, file=sys.stderr)
                print("WARNING: if specified, scan-cat indel/mutation probabilites will be ignored!", file=sys.stderr)
            if sdlog and model!=4:
                print("WARNING: --sdlog should not be used with this chosen model! The given values will be ignored!", file=sys.stderr)
            if genomic:
                print("WARNING: indel/mutation positions will be compared to genomic coordinates!", file=sys.stderr)
            else:
                print("WARNING: indel/mutation positions will be compared to fasta coordinates!", file=sys.stderr)
            print("WARNING: out of the range indel/mutation positions will be ignored!", file=sys.stderr)
        print("ploidy=", ploidy, file=sys.stderr)
        if seed:
            print("WARNING: the random seed is set to {}!".format(seed), file=sys.stderr)
        return(None)
    except:
        raise

def RemoveDuplicate(Sortedlist):	  
	"""Removes duplicates from a sorted list"""
	_Newlist=[]
	if not Sortedlist:
		return(_Newlist)
	else:
		iterable=(x for x in Sortedlist)
		prev=None
		for element in iterable:  
			if (element == prev):
				pass
			else:
				_Newlist.append(element)
				prev=element
		return(_Newlist)

def rnd_end(seq,murate,inrate,delrate):  
	"""Random mutation/indel end points for Poisson and Guassian processes"""
	_l=len(seq)
	if murate>0:
		if _l>(2/murate):
			_muend=min(_l-1, round(_l-1/murate+(1/murate)*random.random()))
		else:
			_muend=_l-1
	else:
		_muend=float('Nan')
	if inrate>0:		
		if _l>(2/inrate):
			_inend=min(_l-1, round(_l-1/inrate+(1/inrate)*random.random()))
		else:
			_inend=_l-1
	else:
		_inend=float('Nan')
	if delrate>0:	
		if _l>(2/delrate):
			_delend=min(_l-1, round(_l-1/delrate+(1/delrate)*random.random()))
		else:
			_delend=_l-1
	else:
		_delend=float('Nan')
	return([_muend,_inend,_delend])

def rnd_position(
    seq: Sequence,
    rates: List[float],
    model: int,
    mean_logs: List[float],
    sd_logs: List[float],
    sample: Optional[List[int]] = None
) -> Tuple[List[int], List[int], List[int]]:
    """Generates lists of mutation/insertion/deletion sites by applying a random process specified by the model parameter."""
    contig_length = len(seq)
    rndmupos, rndinpos, rnddelpos = [], [], []

    if model == 1:  # Poisson process
        for idx, rate in enumerate(rates):
            if rate < 0 or rate >= 1:
                raise ValueError("Rates must be >=0 and <1 with Poisson model!")
            pos = []
            n = 0
            while rate > 0 and n < contig_length:
                # Weighted random rate
                rate_val = weighted_choice([
                    (rate / 4., 1. / 6),
                    (rate, 2. / 3),
                    (rate / 4. + 0.75, 1. / 6)
                ])
                n += max(1, int(round(random.expovariate(rate_val))))
                if n < contig_length:
                    pos.append(n)
            if idx == 0:
                rndmupos = pos
            elif idx == 1:
                rndinpos = pos
            else:
                rnddelpos = pos

    elif model == 2:  # Distance model
        starts = rnd_start(seq, *rates)
        ends = rnd_end(seq, *rates)
        for idx, rate in enumerate(rates):
            if rate < 0 or rate >= 1:
                raise ValueError("Rates must be >=0 and <1 with distance model!")
            if rate == 0:
                continue
            n_points = int(rate * contig_length)
            if n_points == 0:
                continue
            lin = np.linspace(starts[idx], ends[idx], n_points)
            adjust = np.sqrt(0.3 / float(rate)) * np.random.standard_normal(size=len(lin))
            pos = [min(max(0, int(round(x + y))), contig_length - 1) for x, y in zip(lin, adjust)]
            if idx == 0:
                rndmupos = pos
            elif idx == 1:
                rndinpos = pos
            else:
                rnddelpos = pos

    elif model == -1:  # Empirical
        if not sample:
            raise ValueError("Empirical sample is required for empirical model.")
        _rate_for = ['Mutation', 'Insertion', 'Deletion']
        if any(r < 0 or r > 1 for r in rates):
            raise ValueError("Empirical rates must be between 0 and 1.")
        starts = np.array(rnd_start(seq, *rates))
        ends = np.array(rnd_end(seq, *rates))
        _start = int(np.mean(starts[~np.isnan(starts)]))
        _end = int(np.mean(ends[~np.isnan(ends)]))
        choices = [(val, 1) for val in sample]
        current_rnd = []
        while _start < _end:
            current_rnd.append(_start)
            _start = min(_end, _start + weighted_choice(choices))
        current_rnd.append(_end)
        for x in current_rnd:
            event = weighted_choice(list(zip(_rate_for, rates)))
            if event == 'Mutation':
                rndmupos.append(x)
            elif event == 'Insertion':
                rndinpos.append(x)
            else:
                rnddelpos.append(x)

    elif model == 4:  # Lognormal
        for idx, (mean_log, sd_log) in enumerate(zip(mean_logs, sd_logs)):
            if mean_log < 0 or sd_log <= 0:
                raise ValueError("Mean log must be >=0 and sd_log >0 for lognormal model!")
            pos = []
            n = 0
            while mean_log > 0 and n < contig_length:
                n += max(1, int(round(np.random.lognormal(mean_log, sd_log))))
                if n < contig_length:
                    pos.append(n)
            if idx == 0:
                rndmupos = pos
            elif idx == 1:
                rndinpos = pos
            else:
                rnddelpos = pos

    # Remove duplicates and sort
    rndmupos = sorted(set(rndmupos))
    rndinpos = sorted(set(rndinpos))
    rnddelpos = sorted(set(rnddelpos))
    return rndmupos, rndinpos, rnddelpos

def rnd_position_perhomolo(
    dosages: List[float],
    rndpos: List[int],
    ploidy: int
) -> List[List[int]]:
    """Assigns mutation positions to each homologue so that dosages follow the given distribution."""
    rnd_pos_perhomo = [[] for _ in range(ploidy)]
    if rndpos:
        dos_perc = OrderedDict((i + 1, perc) for i, perc in enumerate(dosages))
        dos_per_pos = OrderedDict((pos, weighted_choice(list(dos_perc.items()))) for pos in rndpos)
        homolos_per_pos = [
            random.sample(range(ploidy), dos_per_pos[pos]) for pos in rndpos
        ]
        for i in range(1, len(rndpos)):
            coin = weighted_choice([('heads', 0.4), ('tails', 0.6)])
            if (dos_per_pos[rndpos[i]] == dos_per_pos[rndpos[i - 1]] and coin == 'heads'):
                homolos_per_pos[i] = homolos_per_pos[i - 1][:]
            elif (dos_per_pos[rndpos[i]] > dos_per_pos[rndpos[i - 1]] and coin == 'heads'):
                diff = set(homolos_per_pos[i]) - set(homolos_per_pos[i - 1])
                homolos_per_pos[i] = homolos_per_pos[i - 1] + random.sample(list(diff), dos_per_pos[rndpos[i]] - dos_per_pos[rndpos[i - 1]])
            elif coin == 'heads':
                homolos_per_pos[i] = random.sample(homolos_per_pos[i - 1], dos_per_pos[rndpos[i]])
        for i, homolos in enumerate(homolos_per_pos):
            for homolo in homolos:
                rnd_pos_perhomo[homolo].append(rndpos[i])
    return rnd_pos_perhomo

def rnd_start(seq: Sequence, murate: float, inrate: float, delrate: float) -> List[float]:
    """Random mutation/indel start points for Poisson and Gaussian processes."""
    _l = len(seq)
    def start(rate):
        if rate > 0:
            return max(0, round((1 / rate) * random.random())) if _l > (1 / rate) else 0
        return float('nan')
    return [start(murate), start(inrate), start(delrate)]

def sorted_union(a: List[Any], b: List[Any]) -> List[Any]:
    """Return the sorted union of two lists."""
    return sorted(set(a) | set(b))

def weighted_choice(choices: List[Tuple[Any, float]]) -> Any:
    """Random choice from a weighted set."""
    if not choices:
        raise ValueError("weighted_choice() needs a non-empty list of (value, weight) pairs!")
    values, weights = zip(*choices)
    total = sum(weights)
    cum_weights = []
    cumsum = 0
    for w in weights:
        cumsum += w
        cum_weights.append(cumsum)
    x = random.random() * total
    i = bisect(cum_weights, x)
    return values[i]

if __name__ == "__main__":
    args = getinput(sys.argv[1:])
    try:
        sample = None
        seed = None
        handle = None
        out_handle = None
        index = None
        out_index = None
        varhaplo = None
        varhaplo2 = None
        suffix = 'hap'
        tmp = ('',)
        tmpindex = ('',)
        final = ('',)
        finalindex = ('',)
        fastafile = args[0]
        indexfile = fastafile + ".fai"
        outputbase = args[1]
        muprob = args[2]
        insertprob = args[3]
        deleteprob = args[4]
        muposlst = args[5]
        if muposlst and isinstance(muposlst[0], int):
            muposlst = [muposlst]
        inposlst = args[6]
        if inposlst and isinstance(inposlst[0], int):
            inposlst = [inposlst]
        delposlst = args[7]
        if delposlst and isinstance(delposlst[0], int):
            delposlst = [delposlst]
        mualphabet = args[8]
        verbose = args[9]
        _random = args[10]
        length = args[11]
        randalphabet = args[12]
        content = args[13]
        ploidy = args[14]
        contig_names = args[15]
        linelens = args[16]
        genomic = args[17]
        model = args[18]
        rates = args[19]
        dosage = args[20]
        sdlog = args[21]
        emp_filename = args[22]
        seed = args[23]
        region_start = args[24]
        region_end = args[25]
                
        if seed:
            np.random.seed(int(seed))
            random.seed(seed)
        Contiglist = []
        ContigLengths = []
        if _random:
            n_contig = 0
            linelen = linelens
            if content:
                content_sum = sum(content.values())
                if not all(letter in set(randalphabet) for letter in content.keys()):
                    raise CustomException('ERROR: content letters not in the given alphabet!')
                elif not all(0 <= value <= 1 for value in content.values()):
                    raise CustomException('ERROR: content weights could not be negative or >1')
                elif round(content_sum, 4) > 1:
                    raise CustomException('ERROR: content weights must sum up to 1!')
                elif not all(letter in set(content.keys()) for letter in randalphabet):
                    content_lack = list(set(randalphabet) - set(content.keys()))
                    for dummy in content_lack:
                        content[dummy] = (1 - content_sum) / len(content_lack)
                    if verbose:
                        print("WARNING: weights were not specified for some alphabet letters! The unspecified content weights will be equally set so that the total weight becomes one.", file=sys.stderr)
                elif content_sum != 1:
                    raise CustomException('ERROR: content weights must sum up 1!')
        else:
            if inposlst:
                if len(inposlst) > ploidy:
                    if verbose:
                        print("WARNING: You have specified more insertion position lists than the ploidy!\
	The extra's will be ignored!", file=sys.stderr)
                    inposlst = list(list(poslst) for poslst in inposlst[0:ploidy])
                elif len(inposlst) < ploidy:
                    if verbose:
                        print("WARNING: You have specified fewer insertion position lists than the ploidy!\
	The last position list will be used for the next coming homologues!!", file=sys.stderr)
                    for i in range(0, ploidy - len(inposlst)):
                        inposlst = inposlst + [list(inposlst[-1])]
                else:
                    inposlst = list(list(poslst) for poslst in inposlst)
            
            if delposlst:	
                if len(delposlst) > ploidy:
                    if verbose:
                        print("WARNING: You have specified more deletion position lists than the ploidy!\
	The extra's will be ignored!", file=sys.stderr)
                    delposlst = list(list(poslst) for poslst in delposlst[0:ploidy])
                elif len(delposlst) < ploidy:
                    if verbose:
                        print("WARNING: You have specified fewer deletion position lists than the ploidy!\
	The last position list will be used for the next coming homologues!", file=sys.stderr)
                    for i in range(0, ploidy - len(delposlst)):
                        delposlst = delposlst + [list(delposlst[-1])]
                else:
                    delposlst = list(list(poslst) for poslst in delposlst)
            
            if muposlst:
                if len(muposlst) > ploidy:
                    if verbose:
                        print("WARNING: You have specified more mutation position lists than the ploidy!\
	The extra's will be ignored!", file=sys.stderr)
                    muposlst = list(list(poslst) for poslst in muposlst[0:ploidy])
                elif len(muposlst) < ploidy:
                    if verbose:
                        print("WARNING: You have specified fewer mutation position lists than the ploidy!\
	The last position list will be used for the next coming homologues!", file=sys.stderr)
                    for i in range(0, ploidy - len(muposlst)):
                        muposlst = muposlst + [list(muposlst[-1])]
                else:
                    muposlst = list(list(poslst) for poslst in muposlst)		
            if not emp_filename: # if no emp_filename has been specified, the chosen stochastic model if used to generate variants
                if model!=3 and len(rates)<3:
                    rates=rates+list(0 for x in range(3-len(rates)))
                    if verbose:
                        print >> sys.stderr, ('WARNING: non-specified (mean log) mutation/insertion/deletion rates are assumed zero!')	
                elif model!=3 and len(rates)>3:
                    rates=rates[0:3]
                    if verbose:
                        print >> sys.stderr, ('WARNING: (mean log) rates list should get only three elements mutation/insertion/deletion, respectively! Only the first three elements are used and the rest ignored!')
                if model==4 and len(sdlog)<3:
                    sdlog=sdlog+list(1 for x in range(3-len(sdlog)))
                    if verbose:
                        print >> sys.stderr, ("WARNING: non-specified mutation/insertion/deletion sd_log's are assumed 1!")	
                elif model==4 and len(sdlog)>3:
                    sdlog=sdlog[0:3]
                    if verbose:
                        print >> sys.stderr, ("WARNING: sd_log's list should get only three elements mutation/insertion/deletion, respectively! Only the first three elements are used and the rest ignored!")			
                if model!=3:
                    for _rate in rates:
                        if not (isinstance(_rate,int) or isinstance(_rate,float)):
                            raise ValueError('Unallowed values in mutation/indel rates! Only non-negative real values are allowed!')
                        elif _rate<0:
                            raise ValueError('Unallowed values in mutation/indel rates! Only non-negative real values are allowed!')
                        elif _rate>=1 and model!=4:
                            raise ValueError('Unallowed values in mutation/indel rates! All rates must be <1 within the chosen model!')
                if model==4:
                    for _n, _sdlog in enumerate(sdlog):
                        if not (isinstance(_sdlog,int) or isinstance(_sdlog,float)):
                            raise ValueError('Unallowed values in mutation/indel sdlog! Only positive real numbers are allowed!')
                        elif _sdlog<=0 and rates[_n]>0:
                            raise ValueError('Unallowed values in mutation/indel sdlog! Only positive real numbers are allowed!')
            else:
                model=-1 # -1 just means that the empirical distribution is used instead of a stochastic model
                if len(rates)>3:
                    if verbose:
                        print >>sys.stderr, "WARNING: you have specified more than 3 probabilities for mutation, insertion and deletion! The extras will be ignored!"
                elif len(rates)<3:
                    if verbose:
                        print >>sys.stderr, "WARNING: you have specified less than 3 probabilities for mutation, insertion and deletion! The missing probabilities\
will be set equally so that the sum of all probabilities becomes one!"
                if True in [_x<0 for _x in rates]:
                    raise CustomException("ERROR: negative probabilities are not allowed for mutation, insertion and deletion!")
                if (round(sum(rates[0:min(3, len(rates))]),4)>1 or (len(rates)>=3 and round(sum(rates[0:min(3, len(rates))]),4)<1)):
                    raise CustomException("ERROR: the probabilities for a variant to be a mutation, an insertion or a deletion must sum up to 1 for the empirical model!")
                elif len(rates)<3:
                    rates.extend([round(1-sum(rates),5)/(3-len(rates)) for _x in range(0,3-len(rates))])
            if model!=3 and dosage:
                if len(dosage)>ploidy:
                    if verbose:
                        print >> sys.stderr, ("WARNING: You have specified probabilities for dosages which are larger than\
the ploidy level!\nThese will be ignored!")
                    dosage=dosage[0:ploidy]
                dosages_sum=sum(dosage)
                if not all(_perc<=1 and _perc>=0 for _perc in dosage):
                    raise CustomException("ERROR: Dosage probabilities for simplex to ploidy-plex must be between 0 and 1!")
                elif round(dosages_sum,4)>1 or (round(dosages_sum,4)<1 and len(dosage)==ploidy):
                    raise CustomException("ERROR: Dosage probabilities for simplex to ploidy-plex must sum up to 1!")
                elif round(dosages_sum,4)<1:
                    _perc=(1-dosages_sum)/float(ploidy-len(dosage))
                    if verbose:
                        print >> sys.stderr, "WARNING: You have not specified probabilities for dosages more than %i and the sum of probabilities is less than 1!\n\
Missing probabilities will be set equally, so that the total sum of dosage probabilities becomes one!" %len(dosage)
                    for _i in range(len(dosage),ploidy):
                        dosage.append(_perc)
                elif len(dosage)<ploidy:
                    if verbose:
                        print >> sys.stderr, "WARNING: You have not specified probabilities for dosages more than %i!\
Missing probabilities will be set to zero!" %len(dosage)
                    for _i in range(len(dosage),ploidy):
                        dosage.append(0)
                else:
                    pass
            elif model!=3 and not dosage:       # Binomial dosage propabilities will be used if no dosage is specified!
                dosage=[]
                if model==1 or model==2:
                    _p=max(0.05,sum(rates)/(1e-60+sum(1 for x in rates if x>0)))  # Mean change rate for indel/mutation rates which are not zero, is used in the binomial distribution
                else:															 # to determine dosage percentages. A minimum of 0.05 is, however, considered for this probability.
                    _p=sum(exp(y) for y in rates if y>0)/(1e-60+sum(1 for x in rates if x>0))
                    _p=max(0.05, 1./(1e-60+_p))
                for _k in range(1,ploidy+1):	
                    dosage.append(binom.pmf(_k, ploidy, _p))
                dosage=list(_perc/float(sum(dosage)) for _perc in dosage) # Convert the probabilities to percentages
            else:
                pass 
        if verbose:
            print_info()					# Print the input values if verbose is true	
            print(f"Random = {_random}", file=sys.stderr)
        else:
            pass	
        
        new_insertprob=copy.deepcopy(insertprob)				# Start with the given indel/mutation probabilities and \
        new_deleteprob=copy.deepcopy(deleteprob)				# sequence alphabets, which night be updated by each contig \		
        new_muprob=copy.deepcopy(muprob)						# within the fasta file in case it is not completely specified.
        new_mualphabet=copy.deepcopy(mualphabet)				# In such cases, random mutations will be considered and set \
        variant_homolos=list([] for x in range(0,ploidy))		# Homologues of only heterozygous sites, a list of lists\
        rndmupos, rndinpos, rnddelpos = ([],[],[])				# once in the begin for all of the homologues and contigs.
        if model==-1:
            sample = get_empirical_dist(emp_filename)
        for homolo in range(1,ploidy+1):					      
            if verbose:
                print(f"Generating haplotype {homolo} .........................", file=sys.stderr)
            
            outputfile=outputbase+"_"+suffix+str(homolo)+".fa"
            outindexfile=outputfile+".fai"				 		
            final+=(outputfile,)		           	# First writing on a temporary file will protect the already existing files in case of error!
            finalindex += (outindexfile,)  # A temporary output is first made and the results are written to the output only if everything goes smooth!
            with tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8') as outputfile, tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8') as outindexfile:
                pass
            tmp += (outputfile.name,)
            tmpindex += (outindexfile.name,)
            ref_offset = 0  # reference to write output index
            ref_linelen = 0  # reference to write output index
            ref_seqlen = 0  # reference to write output index
            if _random:                                         # The code for random mode starts here
                with open(outputfile.name, 'w') as out_handle, open(outindexfile.name, 'w') as out_index:
                    for i, l in enumerate(length):
                        l = int(l)                              # If l is not integer, ValueError is thrown
                        if content:
                            contig = ''.join(weighted_choice(list(content.items())) for _ in range(l))
                        else:
                            contig = ''.join(random.choice(randalphabet) for _ in range(l))
                        if not contig_names:
                            contig_id = 'Chr' + str(i)          # Automatically generate contig name
                        else:
                            contig_id = str(contig_names[i])    # Read contig names from contig_names input
                        if isinstance(linelens, list):          # Read line lengths from the input
                            linelen = linelens[i]
                        else:
                            linelen = linelens
                        if verbose:
                            print(f"Writing contig {contig_id} to tmp file.") # Write the haplotype onto the output fasta file and make a new index file therefor
                        if ref_offset == 0:
                            new_offset = len(contig_id) + len(os.linesep) + 1
                        else:
                            new_offset = (ref_offset + ref_seqlen +
                                        len(os.linesep) * int(ceil(float(ref_seqlen) / float(ref_linelen))) +
                                        len(contig_id) + len(os.linesep) + 1)
                        ref_seqlen = l                          # update ref_seqlen
                        ref_offset = new_offset                 # update ref_offset
                        ref_linelen = linelen                   # update ref_seqlen
                        out_handle.write('>' + contig_id + os.linesep)
                        out_handle.write((os.linesep).join(wrap(str(contig), width=linelen)) + os.linesep)
                        if verbose:
                            print(f"Indexing contig {contig_id}.")
                        out_index_line = '\t'.join([contig_id, str(l), str(new_offset),
                                                    str(linelen), str(linelen + len(os.linesep))])
                        out_index.write(out_index_line + os.linesep)  # write the index line for the haplotype onto the output index
                        if homolo == 1:
                            n_contig += 1                        # save the number of the contigs
                    if verbose:
                        print(f"{n_contig} generated contigs ready to be written to {final[homolo]}.")
            else:							  				# The code for mutation mode starts here 
                if not inposlst:
                    end_of_inpos=True							# no insertion sites given
                else:
                    inpos=inposlst[homolo-1]
                    if not inpos:
                        end_of_inpos=True						# no insertion sites given
                    else:				
                        inpos.sort()							# sort insertion sites
                        for _i, _inpos in enumerate(reversed(inpos)):
                            if not isinstance(_inpos,int):
                                raise ValueError('Unallowed values in insertion positions! Only integers are allowed!')
                        _i=None
                        del(_inpos)
                        end_of_inpos=False	
                        inpos_ind=0								# define index for inpos 
                if not delposlst:
                    end_of_delpos=True                            # sort deletion sites
                else:
                    delpos=delposlst[homolo-1]
                    if not delpos:
                        end_of_delpos=True    
                    else:
                        delpos.sort()                            # sort deletion sites
                        for _i, _delpos in enumerate(reversed(delpos)):
                            if not isinstance(_delpos,int):
                                raise ValueError('Unallowed values in deletion positions! Only integers are allowed!')
                        _i=None
                        del(_delpos)
                        end_of_delpos=False
                        delpos_ind=0                            # define index for delpos
				
                if not muposlst:
                    end_of_mupos = True                        # no mutation sites given
                else:
                    mupos = muposlst[homolo - 1]
                    if not mupos:
                        end_of_mupos = True
                    else:
                        mupos.sort()                           # sort mutation sites given
                        for _i, _mupos in enumerate(reversed(mupos)):
                            if not isinstance(_mupos, int):
                                raise ValueError('Unallowed values in mutation positions! Only integers are allowed!')
                        del (_i)
                        del (_mupos)
                        end_of_mupos = False
                        mupos_ind = 0                          # define index for mupos
                prev_offset = 0                                 # offset of the previous record, stored for checking purposes
                prev_seqlen = 0                                 # len of the previous record's seq, stored for checking
                prev_linelen = 0                                # line length of the previous record in the fast file, stored for checking
                _new_hap_pos = []                               # variant positions for this haloptype (contig coordinate as in VCF)
                _new_hap_contig = []                            # contig name for each variant
                _new_hap_ref = []                               # reference alleles at each variant positions, - for insertions
                _new_hap_alt = []                               # alternative allele at each variant positions, - for deletion  
                _mu_alleles = []                                # mutation alleles (to be part of _new_hap_alt)
                _ins_alleles = []                               # inserted alleles (to be part of _new_hap_alt)
                _change_list = []                               # The positions of a contig whereon a change has occured
                _variant_homolos = None                         # Dummy variable to be sorted first and then written to variant_homolos[homolo-1]
                if genomic:                                     # Genomic coordinate starts from 1
                    _genomic_offset = 1                        # 1 as python sequence starts from 0, but genomic positions start from 1
                else:
                    pass
                n_contig = 0                                   # Keep track of the number of the contig being processed 
                with open(fastafile, 'r') as handle, open(outputfile.name, 'w') as out_handle, open(indexfile, 'r') as index, open(outindexfile.name, 'w') as out_index:
                    for record in SeqIO.parse(handle, "fasta"):
                        _seq = record.seq
                        _id = record.id
                        _descrip = record.description
                        if verbose:
                            print(f"Processing contig {_id}... ")
                        line_fai = index.readline()
                        columns = line_fai.rstrip().split('\t')  # <-- FIXED HERE
                        try:
                            contig_id = columns[0]                # id defined from .fai
                            if homolo == 1:
                                Contiglist.append(contig_id)       # contig id added to the cotig list 
                            n_contig += 1                         # save the number of the current contigs
                            contig_length = int(columns[1])       # seq length defined from .fai
                            offset = int(columns[2])              # offset defined from .fai
                            linelen = int(columns[3])             # line length in the index
                            lineleneol = int(columns[4])          # line length plus end of the line in the index
                        except ValueError as err:
                            if not err.args: 
                                err.args = ('',)
                            err.args = err.args + ("Error in the index file! Could not \
        convert coordinate\nand length strings to an integer!",)
                            raise 
                        if len(columns) != 5:                     # check the index file
                            raise CustomException("ERROR: fasta index is not formatted correctly!")
                        elif offset <= 0 or contig_length <= 0 or linelen <= 0 or lineleneol <= 0:
                            raise CustomException("ERROR: fasta index contains zero or negative values!")
                        elif contig_id != _id:
                            raise CustomException("ERROR: fasta index and fasta reference contain different id's!")
                        elif contig_length != len(_seq):
                            raise CustomException("ERROR: fasta index and fasta reference do not match!")
                        elif prev_offset == 0 and offset != len(_descrip) + 1 + lineleneol - linelen:
                            raise CustomException("ERROR: fasta index has wrong offset(s)!")
                        elif (prev_offset > 0) and offset != (prev_offset + prev_seqlen +
                            (prev_lineleneol - prev_linelen) * int(ceil(float(prev_seqlen) / float(prev_linelen))) +
                            len(_descrip) + 1 + lineleneol - linelen):
                            raise CustomException("ERROR: fasta index has wrong offset(s)!")
                        else:
                            prev_offset = offset                  # update prev_offset
                            prev_linelen = int(columns[3])        # update prev_linelen
                            prev_lineleneol = int(columns[4])     # update prev_lineleneol
                            prev_seqlen = contig_length           # update prev_seqlen
                            if verbose:
                                print(f"...Contig {_id} matched the index file.")
                        # Apply region slicing if --region-start and --region-end were provided
                        if region_start is not None and region_end is not None:
                            if verbose:
                                print(f"...Contig {_id} is being processed with region slicing {region_start}-{region_end}.")
                            # Restrict processing to a specific region
                            # region_start and region_end are 1-based coordinates
                            REGION_START = region_start
                            REGION_END = region_end
                            # Ensure region bounds are within the contig
                            seq_len = len(_seq)
                            if REGION_START < 1:
                                REGION_START = 1
                            if REGION_END > seq_len:
                                REGION_END = seq_len
                            # Skip contigs that are too short to have any region within requested bounds
                            if REGION_START > seq_len:
                                if verbose:
                                    print(f"...Skipping contig {_id} (length {seq_len}) - too short for region {region_start}-{region_end}.")
                                continue
                            # Convert to 0-based python slice indices
                            slice_start = REGION_START - 1
                            slice_end = REGION_END
                            # Slice the sequence to the requested region and adjust offset and contig_length
                            # Adjusting offset ensures that downstream position-mapping logic (which uses offset)
                            # still works when comparing input positions to the contig region.
                            _seq = _seq[slice_start:slice_end]
                            # Update contig_length and offset to refer to the sliced region
                            contig_length = len(_seq)
                            offset = offset + slice_start
                            # record the processed contig length (only for the first homologue to avoid duplicates)
                            if homolo == 1:
                                ContigLengths.append(contig_length)
                            # Update the record description/header to show the region processed
                            try:
                                _descrip = f"{_descrip} | region={_id}:{REGION_START}-{REGION_END}"
                            except Exception:
                                # fall back if _descrip is not a string
                                _descrip = str(_descrip) + f" | region={_id}:{REGION_START}-{REGION_END}"
                        else:
                            # No region slicing - record full contig length for density calculation
                            if homolo == 1:
                                ContigLengths.append(contig_length)
                        
                        # Generate random mutation positions (for both region-sliced and full contigs)
                        if homolo == 1:                           # Add random positions to the predefined model for poisson and distance model
                            _rndmupos, _rndinpos, _rnddelpos = rnd_position(seq=_seq, rates=rates, model=model, mean_logs=rates, sd_logs=sdlog, sample=sample) # For scancat, the random position sets will be empty
                            _rndmupos, _rndinpos, _rnddelpos = (rnd_position_perhomolo(dosage, _rndmupos, ploidy),
                                                    rnd_position_perhomolo(dosage, _rndinpos, ploidy),
                                                    rnd_position_perhomolo(dosage, _rnddelpos, ploidy))
                            rndmupos.append(_rndmupos)
                            rndinpos.append(_rndinpos)
                            rnddelpos.append(_rnddelpos)
                            
                        new_mupos = list()                        # mupos coordinates must be changed for MuGen.posmu()
                        new_delpos = list()                       # delpos coordinates must be changed for MuGen.posmu()
                        new_inpos = list()                        # inpos coordinates must be changed for MuGen.posmu()
                        if genomic:                               # If genomic coordinates are used for the input indel/mutation positions (default)
                            while (not end_of_mupos) and mupos[mupos_ind] < (_genomic_offset + contig_length):
                                if mupos[mupos_ind] >= _genomic_offset:    
                                    new_mupos.append(mupos[mupos_ind] - _genomic_offset)
                                else:
                                    pass
                                if mupos_ind < (len(mupos) - 1):
                                    mupos_ind += 1
                                else:
                                    end_of_mupos = True

                            while (not end_of_delpos) and delpos[delpos_ind] < (_genomic_offset + contig_length):
                                if delpos[delpos_ind] >= _genomic_offset:    
                                    new_delpos.append(delpos[delpos_ind] - _genomic_offset)
                                else:
                                    pass
                                if delpos_ind < (len(delpos) - 1):
                                    delpos_ind += 1
                                else:
                                    end_of_delpos = True

                            while (not end_of_inpos) and inpos[inpos_ind] < (_genomic_offset + contig_length):
                                if inpos[inpos_ind] >= _genomic_offset:
                                    new_inpos.append(inpos[inpos_ind] - _genomic_offset)
                                else:
                                    pass
                                if inpos_ind < (len(inpos) - 1):
                                    inpos_ind += 1
                                else:
                                    end_of_inpos = True

                            _genomic_offset = _genomic_offset + contig_length
                        else:                                     # If fasta coordinates are used for the input indel/mutation positions (default)
                            while (not end_of_mupos) and mupos[mupos_ind] < (offset + contig_length):
                                if mupos[mupos_ind] >= offset:    
                                    new_mupos.append(mupos[mupos_ind] - offset)
                                else:
                                    pass
                                if mupos_ind < (len(mupos) - 1):
                                    mupos_ind += 1
                                else:
                                    end_of_mupos = True

                            while (not end_of_delpos) and delpos[delpos_ind] < (offset + contig_length):
                                if delpos[delpos_ind] >= offset:    
                                    new_delpos.append(delpos[delpos_ind] - offset)
                                else:
                                    pass
                                if delpos_ind < (len(delpos) - 1):
                                    delpos_ind += 1
                                else:
                                    end_of_delpos = True

                            while (not end_of_inpos) and inpos[inpos_ind] < (offset + contig_length):
                                if inpos[inpos_ind] >= offset:    
                                    new_inpos.append(inpos[inpos_ind] - offset)
                                else:
                                    pass
                                if inpos_ind < (len(inpos) - 1):
                                    inpos_ind += 1
                                else:
                                    end_of_inpos = True
                        # n_contig is already correctly set from the loop - don't override it
                        new_mupos = sorted_union(rndmupos[n_contig - 1][homolo - 1], new_mupos)
                        new_inpos = sorted_union(rndinpos[n_contig - 1][homolo - 1], new_inpos)
                        new_delpos = sorted_union(rnddelpos[n_contig - 1][homolo - 1], new_delpos)
                        hapgen = MuGen(seq=_seq,                  # Make the new haplotype by changing the reference contig 
                            insertprob=new_insertprob, deleteprob=new_deleteprob, muprob=new_muprob, 
                            mualphabet=new_mualphabet, mupos=new_mupos, delpos=new_delpos, inpos=new_inpos, 
                            verbose=verbose)
                        if model == 3:                            # Apply the Scan-Catalyze model
                            hapgen.hapchanger()                   # Changes made to the specified positions as well as randomly for unspecified positions
                            if homolo == 1:
                                new_muprob.update(hapgen.get_muprob())              # The parameters should be updated so that all 
                                new_insertprob.update(hapgen.get_insertprob())      # contigs/haplotypes have the same working parameters
                                new_deleteprob.update(hapgen.get_deleteprob())
                                new_mualphabet.update(hapgen.get_mualphabet())
                        else:
                            hapgen.posmu()                        # Changes made to the specified and randomly generated positions
                            if homolo == 1:
                                new_mualphabet.update(hapgen.get_mualphabet())
                        if verbose:
                            print(f"Writing modified contig {_id} to tmp file.")
                        _new_seq = Seq(str(hapgen.get_hap()))       # Write the haplotype onto the output fasta file and make a new index file therefor
                        if ref_seqlen == 0:
                            new_offset = len(_descrip) + len(os.linesep) + 1
                        else:
                            new_offset = (ref_offset + ref_seqlen +
                                len(os.linesep) * int(ceil(float(ref_seqlen) / float(ref_linelen))) +
                                len(_descrip) + len(os.linesep) + 1)
                        ref_seqlen = len(_new_seq)                # update ref_seqlen
                        ref_offset = new_offset                   # update ref_offset
                        ref_linelen = linelen                     # update ref_seqlen
                        record.seq = _new_seq
                        out_handle.write('>' + _descrip + os.linesep)
                        out_handle.write((os.linesep).join(wrap(str(_new_seq), width=ref_linelen)) + os.linesep)
                        if verbose:
                            print(f"Indexing the modified contig {_id}.\n")
                        out_index_line = '\t'.join([str(contig_id), str(len(_new_seq)), str(new_offset),
                            str(linelen), str(linelen + len(os.linesep))])
                        out_index.write(out_index_line + os.linesep)  # write the index line for the haplotype onto the output index
                        _change_list = (hapgen.get_occuredmupos() +
                                        hapgen.get_occureddelpos() + hapgen.get_occuredinspos())
                        _ins_alleles = hapgen.get_ins_allele()
                        _mu_alleles = hapgen.get_mu_allele()
                        _new_hap_pos.extend(_change_list)             # Update the variation positions of the whole homologue
                        _new_hap_contig.extend([contig_id] * len(_change_list))  # Save the conitg name for the added posotions
                        _new_hap_ref.extend(list(hapgen.get_ref()[x] for x in _change_list))  # Save the reference alleles of the added variants
                        _new_hap_alt.extend(_mu_alleles + ['-'] * len(hapgen.get_occureddelpos()) + _ins_alleles)  # Save the alternative allele of the added variants
                        _change_list = []
                        _mu_alleles = []
                        _ins_alleles = []
                    _genomic_offset = None
                    _variant_homolos = zip(_new_hap_contig, _new_hap_pos,
                                            _new_hap_ref, _new_hap_alt)
                    _variant_homolos = sorted(_variant_homolos, key=lambda tup: tup[1])  # Sort the reference upon variant position
                    _tmphomolos = list([] for i in range(0, n_contig))
                    for ctgposal in _variant_homolos:            # Group the variants of the same contig together
                        _n = 0
                        _flag = True
                        while _n <= (n_contig - 1) and _flag:
                            if ctgposal[0] == Contiglist[_n]:    # Contigs come in the same order as in the fasta file 
                                _tmphomolos[_n].append(ctgposal)
                                _flag = False
                            else:
                                _n += 1
                    _n = None
                    for _tmpctg in _tmphomolos:
                        variant_homolos[homolo - 1] += _tmpctg
                    _tmpctg = None
                    _tmphomolos = None
                    _variant_homolos = None
                if verbose:
                    print(f"{n_contig} modified contigs ready to be written to {final[homolo]}.")
        if not _random:
            _rnddelpos, _rndmupos, _rndinpos = (None, None, None)
            _variant_homolos = []  # Check if any variants have been detected at all
            for x in variant_homolos:  # In case no change probability or positions are given, no variant will be present \
                _variant_homolos = _variant_homolos + x  # in the generated haplotypes. Then the _varianthaplos file will be empty.
            if _variant_homolos:  # This condition means that at least one homolgue is different from the reference
                _new_hap_contig = list(ploidy * 'x')  # list of the contig names for each variant haplotype
                _new_hap_pos = list(ploidy * 'x')  # list of the variation positions over the contigs for each variant haplotype
                _new_hap_ref = list(ploidy * 'x')  # list of the reference alleles at the variation positions of each variant haplotype
                _new_hap_alt = list(ploidy * 'x')  # list of the alternative alleles at the variation positions of each variant haplotype
                noNULL = None  # Temporary list to replace the empty list of a homologue which is the same as the reference
                _flag = True
                _n = 0
                while _flag and _n < len(variant_homolos):
                    if variant_homolos[_n]:
                        _flag = False
                        noNULL = variant_homolos[_n]
                    else:
                        _n += 1
                _n = None
                _flag = None
                for x in range(0, ploidy):
                    if variant_homolos[x]:
                        _new_hap_contig[x], _new_hap_pos[x], _new_hap_ref[x], _new_hap_alt[x] = zip(*variant_homolos[x])
                    else:
                        _new_hap_contig[x], _new_hap_pos[x], _new_hap_ref[x], _new_hap_alt[x] = zip(*noNULL)
                        _new_hap_alt[x] = copy.deepcopy(_new_hap_ref[x])  # Now you could go on with the non-modified homologue just like the other homolgues
                _polyploid_ref = set()
                for x in range(0, ploidy):  # The unified list of reference positions and alleles for all the variant haplotypes
                    _polyploid_ref = _polyploid_ref.union(set(zip(_new_hap_contig[x],
                                                                _new_hap_pos[x],
                                                                _new_hap_ref[x])))
                _polyploid_ref = list(_polyploid_ref)
                _polyploid_ref = sorted(_polyploid_ref, key=lambda tup: tup[1])  # Sort the reference upon variant position
                polyploid_ref = []
                _tmppolyref = list([] for i in range(0, n_contig))
                for ctgposal in _polyploid_ref:  # Group the variants of the same contig together
                    _n = 0
                    _flag = True
                    while _n <= (n_contig - 1) and _flag:
                        if ctgposal[0] == Contiglist[_n]:  # Contigs come in the same order as in the fasta file
                            _tmppolyref[_n].append(ctgposal)
                            _flag = False
                        else:
                            _n += 1
                _n = None
                for _tmpctg in _tmppolyref:
                    polyploid_ref += _tmpctg
                _tmpctg = None
                _tmppolyref = None
                _polyploid_ref = None
                # Compute variant density and average distance between variants
                try:
                    # map contig -> length
                    contig_len_map = dict(zip(Contiglist, ContigLengths)) if ContigLengths else {}
                    total_bases = sum(contig_len_map.values()) if contig_len_map else 0
                    variant_count = len(polyploid_ref)
                    density = (variant_count / total_bases) if total_bases > 0 else float('nan')
                    # collect inter-variant distances per contig
                    dists = []
                    from collections import defaultdict
                    pos_by_ctg = defaultdict(list)
                    for ctg, pos, ref_allele in polyploid_ref:
                        pos_by_ctg[ctg].append(pos)
                    for ctg, poss in pos_by_ctg.items():
                        poss_sorted = sorted(poss)
                        if len(poss_sorted) > 1:
                            for a, b in zip(poss_sorted, poss_sorted[1:]):
                                dists.append(b - a)
                    avg_dist = (sum(dists) / len(dists)) if dists else float('nan')
                    # print results — include ploidy and the primary mutation rate
                    try:
                        mut_rate_primary = None
                        if isinstance(rates, (list, tuple)) and len(rates) > 0:
                            mut_rate_primary = rates[0]
                        else:
                            mut_rate_primary = rates
                    except Exception:
                        mut_rate_primary = 'NA'

                    print(f"Ploidy: {ploidy}; Mutation rate (primary): {mut_rate_primary}")
                    print(f"Variant density: {variant_count} variants / {total_bases} bases = {density:.6g} variants per base")
                    print(f"Variant density: {density * 1e6:.6g} variants per Mb")
                    if not (avg_dist != avg_dist):
                        print(f"Average distance between variants: {avg_dist:.3f} bp")
                    else:
                        print("Average distance between variants: NA")
                except Exception:
                    # if anything goes wrong, don't crash the pipeline
                    pass
                ContigLen = [0] * n_contig  # Number of lines assigned to each contigeous segement in the haplotypes
                for i, ctg in enumerate(Contiglist):  # Needed to generate the ouput file in the correct format
                    for tup in polyploid_ref:
                        if tup[0] == ctg:
                            ContigLen[i] += 1
                    if i < (n_contig - 1):
                        ContigLen[i + 1] = ContigLen[i]
                ContigLen = RemoveDuplicate(ContigLen)
                Contiglist = None
                ref_code = 0  # The coding of the ref allele
                l = len(polyploid_ref)
                alt_code = list({} for x in range(0, l))  # The coding of the alternative allele. Multiple alleles get 2,3,...
                poly_alt = list('' for x in range(0, l))  # Tha actual alternative alleles to be added to each record of polyploid_ref
                _poly_haps = list([] for x in range(0, ploidy))  # final polyploid haplotypes containing all the variations, to be written to the output file
                l = None
                for x in range(0, ploidy):
                    _index = 0  # Index used to insert alleles in the variant haplotype
                    _n = 0  # Index to update the alternative allele of polyploid_ref
                    l = len(_new_hap_pos[x])
                    for ctg, pos, ref_allele in polyploid_ref:
                        if _index < l:
                            if ctg == _new_hap_contig[x][_index] and pos == _new_hap_pos[x][_index]:
                                if not (_new_hap_alt[x][_index] in set(poly_alt[_n].split(','))) and _new_hap_alt[x][_index] != ref_allele:
                                    poly_alt[_n] += (',' + str(_new_hap_alt[x][_index]))
                                    Code = len(alt_code[_n].keys()) + ref_code + 1
                                    alt_code[_n][_new_hap_alt[x][_index]] = Code  # The first alt_code is ref_code+1 and then the code is incremented by one for each new alt_allele
                                    _poly_haps[x].append((ctg, pos, ref_allele, alt_code[_n].get(_new_hap_alt[x][_index])))
                                elif _new_hap_alt[x][_index] != ref_allele:
                                    _poly_haps[x].append((ctg, pos, ref_allele, alt_code[_n].get(_new_hap_alt[x][_index])))
                                else:
                                    _poly_haps[x].append((ctg, pos, ref_allele, ref_code))
                                _index += 1
                            else:  # If the alternative is not present in this particular haplotype, i.e. \
                                _poly_haps[x].append((ctg, pos, ref_allele, ref_code))  # it has been present in another haplotype, the alternative allele will \
                        else:  # be the same as the reference allele.
                            _poly_haps[x].append((ctg, pos, ref_allele, ref_code))
                        _n += 1
                l = None
                Code = None
                _n = None
                if verbose:
                    print(f"Writing {ploidy} variant haplotyes to file...")
                with tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8') as varhaplo:  # Write the variant haplotypes to a file, first temporarily and after adjustment permanantly
                    varhaplo.write('var_id' + '\t' + 'contig' + '\t' + 'varpos' + '\t' + 'ref_allele' + '\t' + 'alt_allele' + '\t' +
                                '\t'.join(list('hap_' + str(x) for x in range(1, ploidy + 1))) + os.linesep)
                    for i, ctgposal in enumerate(polyploid_ref):
                        ctgposal = (ctgposal[0],) + (ctgposal[1] + 1,) + (ctgposal[2],)  # ADD 1 TO THE VARIANT POSITION SO THAT THE POSITION THAT APPEARS IN THE VARIANT HAPLO FILE CORRESPONDS WITH THE ONE IN VCF FILE
                        string = list(ctgposal) + [str(poly_alt[i][1:])] + (list(_poly_haps[x][i][3] for x in range(0, ploidy)))
                        varhaplo.write('\t'.join(str(x) for x in string) + os.linesep)
                    varhaplo.flush()
                string = None

            else:  # This condition means that no change has been made to the reference. All homologues will be the same as the reference.
                with tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8') as varhaplo:  # Write an empty file for variant haplotypes, first temporarily and after exception check permanantly
                    varhaplo.write('var_id' + '\t' + 'contig' + '\t' + 'varpos' + '\t' + 'ref_allele' + '\t' + 'alt_allele' + '\t' +
                                '\t'.join(list('hap_' + str(x) for x in range(1, ploidy + 1))) + os.linesep)
                ContigLen = [0]
            _variant_homolos = None
    except IOError as e:
        print(f"I/O error({e.errno}): {e.strerror} {e.filename}", file=sys.stderr)
        if out_index:
            out_index.close()
        for i in range(1, len(tmpindex)):
            os.remove(tmpindex[i])
        if out_handle:
            out_handle.close()
        for i in range(1, len(tmp)):
            os.remove(tmp[i])
        if varhaplo:
            varhaplo.close()
            os.remove(varhaplo.name)
        if handle:
            handle.close()
        if index:
            index.close()
        sys.exit(2)
    except CustomException as instance:
        print(instance, file=sys.stderr)
        if out_index:
            out_index.close()
        for i in range(1, len(tmpindex)):
            os.remove(tmpindex[i])
        if out_handle:
            out_handle.close()
        for i in range(1, len(tmp)):
            os.remove(tmp[i])
        if varhaplo:
            varhaplo.close()
            os.remove(varhaplo.name)
        if handle:
            handle.close()
        if index:
            index.close()
        sys.exit(2)
    except ValueError as e:
        print("ERROR: " + str(type(e)) + '\n' + '\n'.join(e.args), file=sys.stderr)
        if out_index:
            out_index.close()
        for i in range(1, len(tmpindex)):
            os.remove(tmpindex[i])
        if out_handle:
            out_handle.close()
        for i in range(1, len(tmp)):
            os.remove(tmp[i])
        if varhaplo:
            varhaplo.close()
            os.remove(varhaplo.name)
        if handle:
            handle.close()
        if index:
            index.close()
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
        traceback.print_exception(exc_type, exc_value, exc_traceback,
                                limit=2, file=sys.stdout)
        sys.exit(2)
    except Exception:
        print("Unexpected Error:", file=sys.stderr)
        if out_index:
            out_index.close()
        for i in range(1, len(tmpindex)):
            os.remove(tmpindex[i])
        if out_handle:
            out_handle.close()
        for i in range(1, len(tmp)):
            os.remove(tmp[i])
        if varhaplo:
            varhaplo.close()
            os.remove(varhaplo.name)
        if handle:
            handle.close()
        if index:
            index.close()
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback, limit=1, file=sys.stderr)
        traceback.print_exception(exc_type, exc_value, exc_traceback,
                                limit=2, file=sys.stderr)
        sys.exit(2)
    else:
        try:
            for i in range(1, len(final)):
                move(tmpindex[i], finalindex[i])
                move(tmp[i], final[i])
            if _random or varhaplo is None:
                pass
            else:
                _tmpname = varhaplo.name
                with open(_tmpname, 'r') as varhaplo_file, open(outputbase + "_varianthaplos.txt", 'w+') as varhaplo2:
                    _n = 0
                    header = ''
                    for i, line in enumerate(varhaplo_file):
                        if i in ContigLen[:-1]:
                            if i > 0:
                                varhaplo2.write(str(i) + '\t' + line +
                                                'Block length  ' + str(ContigLen[_n] - ContigLen[_n - 1]) + '\t' + header)
                                _n += 1
                            else:
                                header = line
                                varhaplo2.write('Block length  ' + str(ContigLen[1]) + '\t' + header)
                                _n += 2
                        else:
                            if i > 0:
                                varhaplo2.write(str(i) + '\t' + line)
                            else:
                                header = line
                                varhaplo2.write('Block length  ' + str(ContigLen[0]) + '\t' + header)
                                _n += 1
                os.remove(_tmpname)
                _tmpname = None
                _n = None
        except Exception:
            if varhaplo:
                os.remove(varhaplo.name)
                varhaplo.close()
            if varhaplo2:
                varhaplo2.close()
                os.remove(outputbase + "_varianthaplos.txt")
            raise
        else:
            if verbose:
                print(f"Modified contigs succesfully written to {', '.join(final[1:])}.")
                print(f"Indexes successfully written to {', '.join(finalindex[1:])}.")
                if not _random and varhaplo is not None:
                    print(f"Variant haplotypes successfully written to {outputbase + '_varianthaplos.txt'}")
            sys.exit(0)
