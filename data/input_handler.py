import vcf
import re
import pandas as pd
import subprocess
import os
from utils.utils import wsl_available
import shutil

class InputHandler:

    def __init__(self, args):
        # Initialize with necessary attributes
        self.data_path = args.data_path
        self.genotype_path = args.genotype_path
        self.ploidy = args.ploidy
        self.alleles = [int(a) for a in args.alleles] if args.alleles is not None else self.compute_alleles()
        self.genotype = self.parse_genotype()
        self.vcf_path = args.vcf_path
        self.vcf_df = self.load_vcf(self.vcf_path)
        self.root_dir = args.root_dir
        self.output_path = args.output_path
        self.bam_path = args.bam_path
        data_from_bam = self.bam2fragmentfile()
        
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

    def bam2fragmentfile(self):
        # self.vcf_df = self.load_vcf(self.vcf_path)
        # if self.data_path == None:
        test_path = self.convertBAM(self.bam_path, self.vcf_path, self.output_path, self.root_dir)
        print('aaaaaaaaaaaa')
        # self.G, self.fragments = loadFragments(self.data_path, self.vcf_df, self.ploidy)
        return test_path

    def validate_input(self):
        # Implement input validation logic
        pass

    def load_vcf(self, vcf_filename):
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


    def convertBAM(self, bam_filename, vcf_filename, output_dir, root_dir):
        '''
        Converts a BAM or SAM file into a fragment file
        :param bam_filename: input BAM or SAM file to be converted
        :return: the path to the fragment file
        '''
        # TODO: call the C++ code to extract fragment file
        prefix=""
        out_filename = os.path.join(output_dir, "bamfrags.txt")
        print('before')
        if wsl_available():
            prefix = "wsl"
            print('before2')
            root_dir = subprocess.check_output(["wsl", "wslpath", "-a", root_dir]).strip().decode()
            out_filename = subprocess.check_output(["wsl", "wslpath", "-a", out_filename]).strip().decode()
        print(root_dir)
        print(prefix)
        print(out_filename)
        print(bam_filename)
        print(vcf_filename)
        print(root_dir+"/extract-poly-src/build/extractHAIRS")
        command = [prefix, os.path.join(root_dir, "extract_poly/build/extractHAIRS"), "--bam", bam_filename,
                   "--vcf", vcf_filename, "--out", out_filename]
        print("Executing command:", ' '.join(command))
        subprocess.check_output(['ls -l', os.path.join(root_dir, "extract_poly/build")])
        print("Executing command:", ' '.join(command))
        subprocess.check_output(command)
        
        # subprocess.check_call(
        #     [prefix, root_dir + "/extract-poly/build/extractHAIRS", "--bam", bam_filename, "--vcf", vcf_filename,
        #      "--out", out_filename])
        # subprocess.check_call(
        #         [prefix, root_dir+"/../extract-poly-src/build/extractHAIRS", "--bam", bam_filename, "--vcf", vcf_filename,
        #          "--out", out_filename])
        os.makedirs(os.path.dirname(out_filename), exist_ok=True)

        try:
            subprocess.check_call(command)
        except subprocess.CalledProcessError as e:
            print("Error occurred:", e)
        except Exception as e:
            print("An unexpected error occurred:", e)
        
        
        
        
        
        subprocess.check_call(
                [prefix, root_dir+"/extract-poly/build/extractHAIRS", "--bam", bam_filename, "--vcf", vcf_filename,
                 "--out", out_filename])
        if wsl_available():
            out_filename = subprocess.check_output(["wsl", "wslpath", "-a", "-w", out_filename]).strip().decode()
    
        return out_filename

