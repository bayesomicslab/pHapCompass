import vcf
import re
import pandas as pd


class InputHandler:

    def __init__(self, data_path, genotype_path, ploidy, alleles=None):
        # Initialize with necessary attributes
        self.data_path = data_path
        self.genotype_path = genotype_path
        self.ploidy = ploidy
        self.alleles = [int(a) for a in alleles] if alleles is not None else self.compute_alleles()
        self.genotype = self.parse_genotype()
        
    def compute_alleles(self):
        # Implement the logic to compute alleles based on ploidy
        al = set(list(range(self.ploidy)))
        return al
    
    def parse_genotype(self):
        with open(self.genotype_path, 'r') as f:
            genotype = f.readlines()
        genotype = genotype[0]
        return genotype
    
    def get_genotype_positions(self, positions):
        return ''.join([self.genotype[p] for p in positions])
        
    def parse_input(self):
        # Implement input parsing logic
        pass

    def validate_input(self):
        # Implement input validation logic
        pass

    def loadVCF(self, vcf_filename):
        '''
        credit: Derek Aguiar
        Loads in a fragment file created by extract-poly
        :param vcf_filename: The file path to the VCF file
        :return:
        '''
        data = []
        csq_fields = ["impact", "aa.pos", "aa.mut", "nuc.pos", "codon.change"]
        vcf_reader = vcf.Reader(open(vcf_filename))
    
        regex = re.compile("/|\|")
    
        # TODO: test this code on various VCF inputs (multiallelic, insertions, deletions, multiple samples)
        for rec in vcf_reader:
            t = {}
            t["POS"] = rec.POS
            t["REF"] = rec.REF
            t["ALT"] = rec.ALT[0]  # Since VCF has only 1 ALT per position
            for k, v in rec.INFO.items():
                if k == "CSQ":
                    for i, j in zip(v[0].split("|"), csq_fields):
                        if ".pos" in j:
                            i = int(i)
                        t["INFO.CSQ." + j] = i
                else:
                    t["INFO." + k] = v
            for s in rec.samples:
                for f in rec.FORMAT.split(":"):
                    if f == "GT":
                        # assumes unphased TODO: allow for phased inputs
                        for idx, gt in enumerate(regex.split(s[f])):
                            t["GT" + str(idx)] = int(gt)
                    else:
                        t["s." + f] = s[f]
        
            data.append(t)
    
        return pd.DataFrame(data)
