

#!/usr/bin/env python3

# -*- coding: utf-8 -*-
#
#Created on Jul 6 2017
#
#@author: Peter Edge, pedge@eng.ucsd.edu

# https://github.com/vibansal/HapCUT2/blob/master/utilities/LinkFragments.py

# Edited by Sina Majidian s_majidian@yahoo.com



from collections import defaultdict
from sys import argv

single_SNP_frags=False


class fragment:

    def __init__(self, seq, name, barcode, dtype = 0):
        self.seq = seq                         # list of (snp index, allele call, quality score) tuples in this version genomic index is removed 
        self.name = name                       # fragment ID / name
        self.barcode = barcode
        self.dtype = dtype
        self.used = False

    def __str__(self):
        fragstr = ''
        num_pairs = 0
        prev_snp_ix = -2
        qual = ' '
        for snp_ix, allele, q_char in self.seq:

            diff = snp_ix - prev_snp_ix

            if diff == 1:
                fragstr += allele
            else:
                num_pairs += 1
                barcode = self.barcode if self.barcode != None else 'NULL'
                fragstr += ' {} {}'.format(snp_ix+1, allele)

            prev_snp_ix = snp_ix
            qual += q_char

        fragstr += qual

        barcode = self.barcode if self.barcode != None else "NULL"
        #if self.dtype == 0:
            #prefix = '{} {} 0 -1 -1'.format(num_pairs,self.name)
        #elif self.dtype == 2:
            #prefix = '{} {} 2 {} -1'.format(num_pairs,self.name, self.barcode)
        prefix =  '{} {}'.format(num_pairs,self.name)
        fragstr = prefix + fragstr
        return fragstr

# read in a HapCUT2 format fragment file into a list of fragment objects
def read_fragment_matrix(frag_matrix):

    snp_ix = 0
    flist = []

    with open(frag_matrix,"r") as fm:
        for line in fm:
            if len(line) < 2: # 1 aaa 2 ACCGCTTCAGCAGATG-1 -1 6054 1 C
                continue

            el = line.strip().split()

            num_blks      = int(el[0])
            name = el[1]

            dtype = int(el[2])

            if not dtype == 2: # value 2 shows that it's a 10x data output of extracthairs
                print("Input to LinkFragments should be unlinked 10X fragments (datatype 2), obtained by running extractHAIRS with --10X 1 option. Current datatype is {}.".format(dtype))
                exit(1)

            barcode = el[3]

            if barcode == 'NULL':
                barcode = None

            call_list  = el[5:(5+2*num_blks)]              # extract base call part of line
            call_list  = zip(*[iter(call_list)]*2)             # list -> tuple list conversion: credit to http://stackoverflow.com/questions/23286254/convert-list-to-a-list-of-tuples-python
            call_list  = [(int(a)-1, b) for a,b in call_list]  # convert index to 0-based integer
            call_list2 = []

            for ix, blk in call_list:
                curr_ix = ix
                for a in blk:
                    call_list2.append((curr_ix, a))
                    curr_ix += 1

            qlist = el[-1]
            alist= [(a,b,c) for ((a,b),c) in zip(call_list2,qlist)] #  (a,vcf_dict[a],b,c)

            frag = fragment(alist,name,barcode,dtype=2)
            flist.append(frag)

    flist.sort(key=lambda x: x.seq[0][0])

    return flist # 





def link_fragments(flist):

    lines = []

    barcode_to_flist = defaultdict(list)

    for f in flist:

        read_id = f.name if f.name[-3:] != '_MP' else f.name[:-3]

        if f.barcode != None:
            barcode_to_flist[f.barcode].append(f)


    dup_snp_cover = 0


    # start is 0-indexed start point and stop is 1 past the last residue, 0-indexed
    #for chrom,start,stop,barcode in get_molecules(bam_file, curr_chrom, dist=dist):
    num_bad_snps = 0
    num_overlapped_match_snps = 0
    for barcode,barcode_flist in barcode_to_flist.items():

        seen_snps = defaultdict(int)
        bad_snps = set()
        new_fseq = []
        for f in barcode_flist:
            #if f.used:
            #    continue

            used = False
            for (snp_ix, allele_call, qual) in f.seq:
                #print(snp_ix)
                if snp_ix in bad_snps:
                    used = True
                    seen_snps[snp_ix] += 1
                    continue # this snp has a disagreement between fragments

                used = True

                if snp_ix in seen_snps:
                    seen_snps[snp_ix] += 1
                    # find location in the fragment we're building,
                    # that has this double-covered SNP
                    # to resolve the conflict
                    for i in range(len(new_fseq)):
                        if new_fseq[i][0] == snp_ix:

                            if new_fseq[i][1] == allele_call:
                                num_overlapped_match_snps += 1
                                q1 = ord(qual) - 33
                                q2 = ord(new_fseq[i][2]) - 33

                                Q = q1 + q2 # combined quality score
                                if Q > 93:
                                    Q = 93

                                Q_char = chr(33 + Q)

                                new_fseq[i] = (new_fseq[i][0], new_fseq[i][1], Q_char) 
                                #print(new_fseq)

                            else:
                                num_bad_snps += 1
                                del new_fseq[i]
                                bad_snps.add(snp_ix)

                            break
                else:
                    new_fseq.append((snp_ix, allele_call, qual))
                    seen_snps[snp_ix] += 1


            f.used = used #False #


        for k,v in seen_snps.items():
            if v > 1:
                dup_snp_cover += v - 1


        new_fseq.sort(key=lambda x: x[0])
        prev_snp_ix = -1
        for (snp_ix, allele_call, qual) in new_fseq:
            if snp_ix <= prev_snp_ix:
                #import pdb
                #pdb.set_trace()
                print("ERROR",file=sys.stderr)
                exit(1)
            prev_snp_ix = snp_ix


        new_id = str(barcode) #"{}:{}-{}:{}".format(chrom, start+1, stop+1, barcode)
        f = fragment(new_fseq,new_id,None,dtype=0)
        if (len(f.seq) == 0) or (not single_SNP_frags and len(f.seq) < 2):
            continue

        firstpos = f.seq[0][0]
        lines.append((firstpos, str(f)))

    linked_count = 0
    null_count = 0
    unlinked_count = 0
    for f in flist:

        if not f.used:
            if f.barcode == None:
                null_count += 1
            else:
                unlinked_count += 1
        else:
            if f.barcode == None:
                print("linked fragment with no barcode")
                exit(1)
            linked_count += 1

        if not f.used and (len(f.seq) >= 2 or (single_SNP_frags and len(f.seq) >= 1)):
            firstpos = f.seq[0][0]
            f.dtype = 0
            lines.append((firstpos, str(f)))
            f.used = True



    print("  {} fragments linked to larger molecules".format(linked_count))
    print("  {} unlinked fragments with barcodes".format(unlinked_count))
    print("  {} unlinked fragments without barcodes".format(null_count))
    print("  {} duplicate snp-cover. {} bad and {} matched-allele snps ".format(dup_snp_cover,num_bad_snps,num_overlapped_match_snps))


    del flist

    # write new fragments to file

    lines.sort()

    return  lines


# main function
# parse input and run function to call alleles
if __name__ == '__main__':

    hairs_file = argv[1]  #'frag.txt' #
    outfile= argv[2] #'frag_out.txt' 

    flist = read_fragment_matrix(hairs_file)
    #flist[0].barcode
    lines= link_fragments(flist)

    with open(outfile, 'w') as opf:
        for firstpos, line in lines:
            print(line, file=opf)







    
