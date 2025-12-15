"""Mutation Generator (MuGen) class
Written by Ehsan Motazedi, 14-08-2015, Wageningen UR
Updated for Python 3 and Biopython >=1.78 (no Bio.Alphabet)
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
import copy
import random
import sys

class CustomException(Exception):
    def __init__(self, value):
        self.parameter = value

    def __str__(self):
        return repr(self.parameter)

class MuGen(object):
    """Performs mutations and deletion/insertion with desired probability
    and desired structure. Gets a Seq object, a mutation or indel dictionary,
    and the probabilities for each item in those dictionaries.
    insertprob and deleteprob are base specific probabilities of length 4
    mualphabet is a dictionary specifying the possible mutations for each letter of
    the sequence alphabet.
    muprob gives the mutation probability for each letter of the sequence alphabet."""

    def __init__(self, seq, alphaproperty=None, insertprob=None,
                deleteprob=None, mualphabet=None,
                muprob=None, mupos=None, delpos=None, inpos=None,
                verbose=False):
        try:
            self.occureddel = []
            self.occuredmu = []
            self.occuredins = []
            self.inserted_allele = []
            self.alt_allele = []
            if not isinstance(verbose, bool):
                raise CustomException("ERROR: verbose must be set to either True or False. Default is False")
            else:
                self.verbose = verbose
            if isinstance(seq, str):
                self.seq = MutableSeq(seq)
            elif isinstance(seq, Seq):
                self.seq = MutableSeq(str(seq))
            elif isinstance(seq, MutableSeq):
                self.seq = copy.deepcopy(seq)
            else:
                raise CustomException("ERROR: Should provide a Seq or MutableSeq object, or a string sequence!")
            self.alphabet = set(str(self.seq))
            self.ref = str(self.seq)
            self.alphaproperty = None  # No longer used, but kept for compatibility

            self.delpos = list(delpos) if delpos and set(delpos).issubset(set(range(len(self.ref)))) else []
            if delpos and not set(delpos).issubset(set(range(len(self.ref)))):
                raise CustomException("ERROR: Deletion positions exceed the range of the reference or are not positive integers!")
            self.inpos = list(inpos) if inpos and set(inpos).issubset(set(range(len(self.ref)))) else []
            if inpos and not set(inpos).issubset(set(range(len(self.ref)))):
                raise CustomException("ERROR: Insertion positions exceed the range of the reference or are not positive integers!")
            self.mupos = list(mupos) if mupos and set(mupos).issubset(set(range(len(self.ref)))) else []
            if mupos and not set(mupos).issubset(set(range(len(self.ref)))):
                raise CustomException("ERROR: Mutation positions exceed the range of the reference or are not positive integers!")

            # Handle mutation alphabet
            if not mualphabet:
                if self.verbose:
                    print("WARNING: You have specified no mutation alphabet! Mutations are set to random letters!")
                self.mualphabet = {key: ''.join(self.alphabet - {key}) for key in self.alphabet}
            else:
                mualphabet = {str(k): str(v) for k, v in mualphabet.items()}
                for key, value in mualphabet.items():
                    if len(key) != 1:
                        raise CustomException("ERROR: the mutation alphabet deals with point mutations! Only single letters are allowed as keys!")
                    elif key in set(''.join(value)):
                        raise CustomException("ERROR: Wrong mutation values specified! A letter could just be substituted with a different letter for mutation!")
                if set(mualphabet.keys()) == self.alphabet and set(''.join(mualphabet.values())) <= self.alphabet:
                    self.mualphabet = copy.deepcopy(mualphabet)
                elif set(mualphabet.keys()) < self.alphabet and set(''.join(mualphabet.values())) < self.alphabet:
                    if self.verbose:
                        print("WARNING: Mutation is not specified for some letters! Those mutations are set to random letters!")
                    self.mualphabet = copy.deepcopy(mualphabet)
                    for key in self.alphabet - set(mualphabet.keys()):
                        self.mualphabet[key] = ''.join(self.alphabet - {key})
                else:
                    if self.verbose:
                        print("WARNING: Mutation alphabet is not compatible with sequence alphabet! Both alphabets are updated and unspecified mutations are set to random letters!")
                    new_mualphabet = dict()
                    for key, value in mualphabet.items():
                        self.alphabet.add(key)
                        self.alphabet |= (set(''.join(value)) - self.alphabet)
                        new_mualphabet.update({key: value})
                    for key in self.alphabet - set(new_mualphabet.keys()):
                        new_mualphabet[key] = ''.join(self.alphabet - {key})
                    self.mualphabet = copy.deepcopy(new_mualphabet)

            # Handle insertion probability
            if not insertprob:
                self.insertprob = {key: 0 for key in self.alphabet}
            else:
                if set(insertprob.keys()) != self.alphabet:
                    if self.verbose:
                        print("WARNING: Missing/Invalid letter(s) in insertion probability! Probabilities are set to zero for missing letters! Invalid letters are ignored!")
                new_insertprob = {}
                for key, value in insertprob.items():
                    if 0 <= value <= 1:
                        new_insertprob[key] = value
                    else:
                        raise CustomException("ERROR: Insertion probability must be >=0 and <=1!")
                for key in self.alphabet - set(new_insertprob.keys()):
                    new_insertprob[key] = 0
                self.insertprob = copy.deepcopy(new_insertprob)

            # Handle deletion probability
            if not deleteprob:
                self.deleteprob = {key: 0 for key in self.alphabet}
            else:
                if set(deleteprob.keys()) != self.alphabet:
                    if self.verbose:
                        print("WARNING: Missing/Invalid letter(s) in deletion probability! Probabilities are set to zero for missing letters! Invalid letters are ignored!")
                new_deleteprob = {}
                for key, value in deleteprob.items():
                    if 0 <= value <= 1:
                        new_deleteprob[key] = value
                    else:
                        raise CustomException("ERROR: Deletion probability must be >=0 and <=1!")
                for key in self.alphabet - set(new_deleteprob.keys()):
                    new_deleteprob[key] = 0
                self.deleteprob = copy.deepcopy(new_deleteprob)

            # Handle mutation probability
            if not muprob:
                self.muprob = {key: 0 for key in self.alphabet}
            else:
                if set(muprob.keys()) != self.alphabet:
                    if self.verbose:
                        print("WARNING: Missing/Invalid letter(s) in mutation probability! Probabilities are set to zero for missing letters! Invalid letters are ignored!")
                new_muprob = {}
                for key, value in muprob.items():
                    if 0 <= value <= 1:
                        new_muprob[key] = value
                    else:
                        raise CustomException("ERROR: Mutation probability must be >=0 and <=1!")
                for key in self.alphabet - set(new_muprob.keys()):
                    new_muprob[key] = 0
                self.muprob = copy.deepcopy(new_muprob)
        except CustomException as instance:
            print(instance)
            sys.exit(2)
        else:
            if self.verbose:
                print("MuGen object successfully created.\nWARNING: MuGen sequence is case sensitive!")

    def __repr__(self):
        return (f"Haplotype: {self.seq}, \n Reference sequence: {self.ref}, \n Mutation probability: {self.muprob}, \n"
                f"Mutations: {self.mualphabet}, \n Insertion probability: {self.insertprob}, \n Deletion Probability: {self.deleteprob}, \n"
                f"Insertion positions: {self.inpos}, \n Deletion positions: {self.delpos}, \n Mutation positions: {self.mupos} \n")

    def __str__(self):
        return repr(self)

    # Access Methods
    def get_hap(self):
        return self.seq

    def get_ref(self):
        return self.ref

    def get_insertprob(self):
        return self.insertprob

    def get_deleteprob(self):
        return self.deleteprob

    def get_muprob(self):
        return self.muprob

    def get_mualphabet(self):
        return self.mualphabet

    def get_mupos(self):
        return self.mupos

    def get_inpos(self):
        return self.inpos

    def get_delpos(self):
        return self.delpos

    def get_occureddelpos(self):
        return self.occureddel

    def get_occuredmupos(self):
        return self.occuredmu

    def get_occuredinspos(self):
        return self.occuredins

    def get_ins_allele(self):
        return self.inserted_allele

    def get_mu_allele(self):
        return self.alt_allele

    def set_ref(self, ref):
        """Changes the reference sequence of the MuGen object. Could become problematic if the new reference
        has a different length than the current reference, while indel and mutation positions are specified.
        A useful method if reference is a mutable seq entity which is constantly called and changed by other
        methods and classes."""
        try:
            if set(str(ref)).issubset(self.alphabet):
                if not set(self.mupos).issubset(set(range(len(str(ref))))):
                    raise CustomException("ERROR: Mutation positions exceed the range of the new reference!")
                elif not set(self.inpos).issubset(set(range(len(str(ref))))):
                    raise CustomException("ERROR: Insertion positions exceed the range of the new reference!")
                elif not set(self.delpos).issubset(set(range(len(str(ref))))):
                    raise CustomException("ERROR: Deletion positions exceed the range of the new reference!")
                else:
                    self.ref = str(ref)
            else:
                raise CustomException("ERROR: the new reference is not compatible with the current alphabet!")
        except CustomException as instance:
            print("Failed to update the reference!")
            print(instance)
        except Exception as e:
            print("Failed to update the reference!")
            raise
        else:
            if self.verbose:
                print("The reference sequence has been successfully updated!")

    def set_pos(self, inpos=None, delpos=None, mupos=None):
        """Changes the insertion, deletion and substitution sites of the MuGen object. A useful method if
        posmu and probmu methods are constantly called."""
        try:
            changedel = False
            changein = False
            changemu = False
            if delpos is not None:
                if set(delpos).issubset(set(range(len(self.ref)))):
                    changedel = True
                else:
                    raise CustomException("ERROR: New deletion positions exceed the range of the reference or are not positive integers!")
            if inpos is not None:
                if set(inpos).issubset(set(range(len(self.ref)))):
                    changein = True
                else:
                    raise CustomException("ERROR: New insertion positions exceed the range of the reference or are not positive integers!")
            if mupos is not None:
                if set(mupos).issubset(set(range(len(self.ref)))):
                    changemu = True
                else:
                    raise CustomException("ERROR: New mutation positions exceed the range of the reference or are not positive integers!")
            if changedel:
                self.delpos = list(delpos)
            if changein:
                self.inpos = list(inpos)
            if changemu:
                self.mupos = list(mupos)
        except CustomException as instance:
            print("Failed to update indel and mutation positions!")
            print(instance)
        except Exception as e:
            print("Failed to update indel and mutation positions!")
            raise
        else:
            if self.verbose:
                print("Indel and mutation positions successfully updated!")

    def set_prob(self, insertprob=None, deleteprob=None, muprob=None):
        """Changes the insertion, deletion and mutation probabilities of the MuGen object. A useful method if
        posmu and probmu methods are constantly called."""
        try:
            noinsert = -1
            nodel = -1
            nomu = -1
            if insertprob is None:
                noinsert = 0
            elif not insertprob:
                noinsert = 1
            elif set(insertprob.keys()) != self.alphabet:
                if self.verbose:
                    print("WARNING: Missing/Invalid letter(s) in insertion probability! Probabilities are set to zero for missing letters! Invalid letters are ignored!")
                new_insertprob = {}
                for key, value in insertprob.items():
                    if 0 <= value <= 1:
                        new_insertprob[key] = value
                    else:
                        raise CustomException("ERROR: Insertion probability must be >=0 and <=1!")
                for key in self.alphabet - set(new_insertprob.keys()):
                    new_insertprob[key] = 0
            else:
                new_insertprob = copy.deepcopy(insertprob)
            if deleteprob is None:
                nodel = 0
            elif not deleteprob:
                nodel = 1
            elif set(deleteprob.keys()) != self.alphabet:
                if self.verbose:
                    print("WARNING: Missing/Invalid letter(s) in deletion probability! Probabilities are set to zero for missing letters! Invalid letters are ignored!")
                new_deleteprob = {}
                for key, value in deleteprob.items():
                    if 0 <= value <= 1:
                        new_deleteprob[key] = value
                    else:
                        raise CustomException("ERROR: Deletion probability must be >=0 and <=1!")
                for key in self.alphabet - set(new_deleteprob.keys()):
                    new_deleteprob[key] = 0
            else:
                new_deleteprob = copy.deepcopy(deleteprob)
            if muprob is None:
                nomu = 0
            elif not muprob:
                nomu = 1
            elif set(muprob.keys()) != self.alphabet:
                if self.verbose:
                    print("WARNING: Missing/Invalid letter(s) in mutation probability! Probabilities are set to zero for missing letters! Invalid letters are ignored!")
                new_muprob = {}
                for key, value in muprob.items():
                    if 0 <= value <= 1:
                        new_muprob[key] = value
                    else:
                        raise CustomException("ERROR: Mutation probability must be >=0 and <=1!")
                for key in self.alphabet - set(new_muprob.keys()):
                    new_muprob[key] = 0
            else:
                new_muprob = copy.deepcopy(muprob)
            if nodel == 0:
                pass
            elif nodel == 1:
                self.deleteprob = {key: 0 for key in self.alphabet}
            else:
                self.deleteprob = copy.deepcopy(new_deleteprob)
            if nomu == 0:
                pass
            elif nomu == 1:
                self.muprob = {key: 0 for key in self.alphabet}
            else:
                self.muprob = copy.deepcopy(new_muprob)
            if noinsert == 0:
                pass
            elif noinsert == 1:
                self.insertprob = {key: 0 for key in self.alphabet}
            else:
                self.insertprob = copy.deepcopy(new_insertprob)
        except CustomException as instance:
            print(instance)
            print("Failed to update indel and mutation probabilities!")
        except Exception:
            print("Failed to update indel and mutation probabilities!")
            raise
        else:
            if self.verbose:
                print("Indel and mutation probabilities successfully updated!")

    def set_mualphabet(self, mualphabet=None):
        """Changes the mutation alphabet of the MuGen object. A useful method if posmu and probmu methods
        are constantly called."""
        try:
            if not mualphabet:
                if self.verbose:
                    print("WARNING: You have specified no mutation alphabet! Mutations are set to random letters!")
                self.mualphabet = {key: ''.join(self.alphabet - {key}) for key in self.alphabet}
            else:
                mualphabet = {str(k): str(v) for k, v in mualphabet.items()}
                for key, value in mualphabet.items():
                    if len(key) != 1:
                        raise CustomException("ERROR: the mutation alphabet deals with point mutations! Only single letters are allowed as keys!")
                    elif key in set(''.join(value)):
                        raise CustomException("ERROR: Wrong mutation values specified! A letter could just be substituted with a different letter for mutation!")
                if set(mualphabet.keys()) == self.alphabet and set(''.join(mualphabet.values())) <= self.alphabet:
                    self.mualphabet = copy.deepcopy(mualphabet)
                elif set(mualphabet.keys()) < self.alphabet and set(''.join(mualphabet.values())) < self.alphabet:
                    if self.verbose:
                        print("WARNING: Mutation is not specified for some letters! Those mutations are set to random letters!")
                    self.mualphabet = copy.deepcopy(mualphabet)
                    for key in self.alphabet - set(mualphabet.keys()):
                        self.mualphabet[key] = ''.join(self.alphabet - {key})
                else:
                    if self.verbose:
                        print("WARNING: Mutation alphabet is not compatible with sequence alphabet! Both alphabets are updated and unspecified mutations are set to random letters!")
                    new_mualphabet = dict()
                    for key, value in mualphabet.items():
                        self.alphabet.add(key)
                        self.alphabet |= (set(''.join(value)) - self.alphabet)
                        new_mualphabet.update({key: value})
                    for key in self.alphabet - set(new_mualphabet.keys()):
                        new_mualphabet[key] = ''.join(self.alphabet - {key})
                    self.mualphabet = copy.deepcopy(new_mualphabet)
        except CustomException as instance:
            print(instance)
            print("Mualphabet could not be updated!")
        except Exception:
            print("Mualphabet could not be updated!")
            raise
        else:
            if self.verbose:
                print("Mualphabet has been successfully updated!")

    def probmu(self):
        """Operates on a MuGen object, and returns a Seq object obtained by making random changes
        to the reference sequence of the MuGen object, using the probabilities given to MuGen."""
        self.occuredmu = []
        self.occureddel = []
        self.occuredins = []
        self.inserted_allele = []
        self.alt_allele = []
        seq_list = []
        for site, base in enumerate(self.ref):
            if site in set(self.mupos) | set(self.inpos) | set(self.delpos):
                seq_list.append(base)  # No change at indel/mutation positions
            else:
                prob = {
                    'ins': self.insertprob.get(base, 0),
                    'del': self.deleteprob.get(base, 0),
                    'sub': self.muprob.get(base, 0)
                }
                error = random.choice(['ins', 'del', 'sub', 'sub'])
                rnd = random.random()
                if rnd < prob.get(error, 0):
                    if error == 'sub':
                        alt = random.choice(self.mualphabet.get(base, base))
                        seq_list.append(alt)
                        self.occuredmu.append(site)
                        self.alt_allele.append(alt)
                    elif error == 'ins':
                        seq_list.append(base)
                        ins = random.choice(list(self.alphabet))
                        seq_list.append(ins)
                        self.occuredins.append(site)
                        self.inserted_allele.append(base + ins)
                    elif error == 'del':
                        self.occureddel.append(site)
                else:
                    seq_list.append(base)
        self.seq = MutableSeq(''.join(seq_list))
        # Sort and zip alleles
        if self.occuredins:
            ins_allele = sorted(zip(self.occuredins, self.inserted_allele), key=lambda tup: tup[0])
            self.occuredins, self.inserted_allele = zip(*ins_allele)
            self.occuredins = list(self.occuredins)
            self.inserted_allele = list(self.inserted_allele)
        else:
            self.inserted_allele = []
            self.occuredins = []
        if self.occuredmu:
            alt_allele = sorted(zip(self.occuredmu, self.alt_allele), key=lambda tup: tup[0])
            self.occuredmu, self.alt_allele = zip(*alt_allele)
            self.occuredmu = list(self.occuredmu)
            self.alt_allele = list(self.alt_allele)
        else:
            self.occuredmu = []
            self.alt_allele = []
        if self.occureddel:
            self.occureddel.sort()
        else:
            self.occureddel = []
        if self.verbose:
            print("WARNING: If indel/mutation positions are specified, MuGen.probmu() makes no change at those sites.\nUse MuGen.posmu() or MuGen.hapchanger() to apply changes at those sites!")
            print("Changes made to the haplotype!")

    def posmu(self):
        """Operates on a MuGen object, and returns a Seq object obtained by making specific changes
        at specific locations on the reference sequence of the MuGen object, using the
        indel and mutation positions already given to MuGen"""
        change = [None] * len(self.ref)
        self.occuredmu = []
        self.occureddel = []
        self.occuredins = []
        self.inserted_allele = []
        self.alt_allele = []
        for site in self.inpos:
            change[site] = 'ins'
        for site in self.delpos:
            change[site] = 'del'
        for site in self.mupos:
            change[site] = 'sub'
        seq_list = []
        for site, error in enumerate(change):
            base = self.ref[site]
            if error is None:
                seq_list.append(base)
            elif error == 'sub':
                alt = random.choice(self.mualphabet.get(base, base))
                seq_list.append(alt)
                self.occuredmu.append(site)
                self.alt_allele.append(alt)
            elif error == 'ins':
                seq_list.append(base)
                ins = random.choice(list(self.alphabet))
                seq_list.append(ins)
                self.occuredins.append(site)
                self.inserted_allele.append(base + ins)
            elif error == 'del':
                self.occureddel.append(site)
        self.seq = MutableSeq(''.join(seq_list))
        if self.occuredins:
            ins_allele = sorted(zip(self.occuredins, self.inserted_allele), key=lambda tup: tup[0])
            self.occuredins, self.inserted_allele = zip(*ins_allele)
            self.occuredins = list(self.occuredins)
            self.inserted_allele = list(self.inserted_allele)
        else:
            self.inserted_allele = []
            self.occuredins = []
        if self.occuredmu:
            alt_allele = sorted(zip(self.occuredmu, self.alt_allele), key=lambda tup: tup[0])
            self.occuredmu, self.alt_allele = zip(*alt_allele)
            self.occuredmu = list(self.occuredmu)
            self.alt_allele = list(self.alt_allele)
        else:
            self.occuredmu = []
            self.alt_allele = []
        if self.occureddel:
            self.occureddel.sort()
        else:
            self.occureddel = []
        if self.verbose:
            print("Changes made to the haplotype!")

    def hapchanger(self):
        """Operates on a MuGen object, and returns a Seq object obtained by making random and specified
        changes to the reference sequence of the MuGen object, using the probabilities as well as the
        positions given to MuGen."""
        self.occuredmu = []
        self.occureddel = []
        self.occuredins = []
        self.inserted_allele = []
        self.alt_allele = []
        seq_list = []
        for site, base in enumerate(self.ref):
            if site in set(self.mupos):
                alt = random.choice(self.mualphabet.get(base, base))
                seq_list.append(alt)
                self.occuredmu.append(site)
                self.alt_allele.append(alt)
            elif site in set(self.inpos):
                seq_list.append(base)
                ins = random.choice(list(self.alphabet))
                seq_list.append(ins)
                self.occuredins.append(site)
                self.inserted_allele.append(base + ins)
            elif site in set(self.delpos):
                self.occureddel.append(site)
            else:
                prob = {
                    'ins': self.insertprob.get(base, 0),
                    'del': self.deleteprob.get(base, 0),
                    'sub': self.muprob.get(base, 0)
                }
                error = random.choice(['ins', 'del', 'sub', 'sub'])
                rnd = random.random()
                if rnd < prob.get(error, 0):
                    if error == 'sub':
                        alt = random.choice(self.mualphabet.get(base, base))
                        seq_list.append(alt)
                        self.occuredmu.append(site)
                        self.alt_allele.append(alt)
                    elif error == 'ins':
                        seq_list.append(base)
                        ins = random.choice(list(self.alphabet))
                        seq_list.append(ins)
                        self.occuredins.append(site)
                        self.inserted_allele.append(base + ins)
                    elif error == 'del':
                        self.occureddel.append(site)
                else:
                    seq_list.append(base)
        self.seq = MutableSeq(''.join(seq_list))
        if self.occuredins:
            ins_allele = sorted(zip(self.occuredins, self.inserted_allele), key=lambda tup: tup[0])
            self.occuredins, self.inserted_allele = zip(*ins_allele)
            self.occuredins = list(self.occuredins)
            self.inserted_allele = list(self.inserted_allele)
        else:
            self.inserted_allele = []
            self.occuredins = []
        if self.occuredmu:
            alt_allele = sorted(zip(self.occuredmu, self.alt_allele), key=lambda tup: tup[0])
            self.occuredmu, self.alt_allele = zip(*alt_allele)
            self.occuredmu = list(self.occuredmu)
            self.alt_allele = list(self.alt_allele)
        else:
            self.occuredmu = []
            self.alt_allele = []
        if self.occureddel:
            self.occureddel.sort()
        else:
            self.occureddel = []
        if self.verbose:
            print("Changes made to the haplotype!")

