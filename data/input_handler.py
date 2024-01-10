from cyvcf2 import VCF
import re
import pandas as pd
import subprocess
import os
from utils.utils import wsl_available
import shutil

class InputHandler:

    def __init__(self, args):
        # Initialize with necessary attributes
        self.data_path = args.data_path if args.data_path is not None else self.bam2fragmentfile()
        self.genotype_path = args.genotype_path
        self.ploidy = args.ploidy
        self.alleles = [int(a) for a in args.alleles] if args.alleles is not None else self.compute_alleles()
        self.genotype = self.parse_genotype()
        self.vcf_path = args.vcf_path
        self.vcf_df = self.load_vcf(self.vcf_path)
        self.root_dir = args.root_dir
        self.output_path = args.output_path
        self.bam_path = args.bam_path if args.bam_path is not None else None
        # self.data_path = self.convertBAM(self.bam_path, self.vcf_path, self.output_path, self.root_dir)
        # data_from_bam = self.bam2fragmentfile()
        # self.data_path = data_from_bam
        # self.data_path = args.data_path if args.data_path is not None else self.bam2fragmentfile()
    
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
        
        self.data_path = self.convertBAM(self.bam_path, self.vcf_path, self.output_path, self.root_dir)
        # self.G, self.fragments = loadFragments(self.data_path, self.vcf_df, self.ploidy)
        return self.data_path

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
            
        # TODO: test this code on various VCF inputs (multiallelic, insertions, deletions, multiple samples)
        for rec in VCF(vcf_filename):
            t = {}
            t["POS"] = rec.POS
            t["REF"] = rec.REF
            t["ALT"] = rec.ALT[0]  # Since VCF has only 1 ALT per position
            for k, v in rec.INFO:
                if k == "CSQ":
                    for i, j in zip(v[0].split("|"), csq_fields):
                        if ".pos" in j:
                            i = int(i)
                        t["INFO.CSQ." + j] = i
                else:
                    t["INFO." + k] = v
            for f in rec.FORMAT:
                if f == "GT":
                    for geno in rec.genotypes:
                        assert geno[-1]==False # assumes unphased TODO: allow for phased inputs
                        for idx, gt in enumerate(geno[:-1]):
                            t["GT" + str(idx)] = int(gt) 
                else:
                    if len(rec.format(f)) ==1:
                        t[f] = rec.format(f)[0]
                    else:
                        for idx,sample_val in enumerate(rec.format(f)):
                            t["s." + f + "."  + str(idx)] = sample_val
        
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
        print('root_dir:', root_dir)
        print('prefix:', prefix)
        print('out_filename:', out_filename)
        print('bam_filename:', bam_filename)
        print('vcf_filename:', vcf_filename)
        
        run_command = ['sh', os.path.join(root_dir, "extract_poly/build", 'test.sh')]
        print(' '.join(run_command))
        subprocess.run(run_command)
        
        subprocess.check_call(run_command)
        # print("ls command on output:")
        # subprocess.run(['ls', os.path.join(root_dir, "extract_poly/build/")])
        # print('ls done')
        # # print(root_dir+"/extract-poly-src/build/extractHAIRS")
        # # readlink_com = 'readlink -f' + ' ' + os.path.join(root_dir, "extract_poly/build/extractHAIRS")
        # readlink_com = 'readlink -f' + ' ' + vcf_filename
        #
        # print('readlink command:', readlink_com)
        # subprocess.run(['readlink', '-f', os.path.join(root_dir, "extract_poly/build/extractHAIRS")])
        # print('done')
        # subprocess.check_call(['readlink', '-f', os.path.join(root_dir, "extract_poly/build/extractHAIRS")])
        #
        # print("ls command on output:")
        # subprocess.run(['ls', os.path.join(root_dir, output_dir)])
        
        # print("mkdir on output:")
        # subprocess.run(['mkdir', os.path.join(root_dir, output_dir, 'test2')])

        command = [prefix, os.path.join(root_dir, "extract_poly/build/extractHAIRS"), "--bam", bam_filename,
                   "--vcf", vcf_filename, "--out", out_filename]
        print("Executing command:", ' '.join(command))
        subprocess.check_call(command)
        # subprocess.run(['ls', '/home/FCAM/mhosseini/HaplOrbit/output/'])
        # subprocess.run(['cat', '/home/FCAM/mhosseini/HaplOrbit/output/test.txt'])
        # subprocess.run(['mkdir', os.path.join(output_dir, 'test2')])

        # print("Executing command:", ' '.join(command))
        

        # subprocess.check_call(
        #     [prefix, root_dir + "/extract-poly/build/extractHAIRS", "--bam", bam_filename, "--vcf", vcf_filename,
        #      "--out", out_filename])
        # subprocess.check_call(
        #         [prefix, root_dir+"/../extract-poly-src/build/extractHAIRS", "--bam", bam_filename, "--vcf",
        #         vcf_filename, "--out", out_filename])
        os.makedirs(os.path.dirname(out_filename), exist_ok=True)

        try:
            subprocess.check_call(command)
        except subprocess.CalledProcessError as e:
            print("Error occurred:", e)
        except Exception as e:
            print("An unexpected error occurred:", e)
            
        # subprocess.check_call(
        #         [prefix, root_dir+"/extract-poly/build/extractHAIRS", "--bam", bam_filename, "--vcf", vcf_filename,
        #          "--out", out_filename])
        if wsl_available():
            out_filename = subprocess.check_output(["wsl", "wslpath", "-a", "-w", out_filename]).strip().decode()
    
        return out_filename

