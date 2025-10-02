import os
import pandas as pd
import random
import pysam
import sys
import numpy as np
# sys.path.append('/home/mok23003/BML/HaplOrbit')
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pickle
# from ..generate_simulated_graphs import generate_quotient_graph
from data.input_handler import InputHandler
from data.configuration import Configuration
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from multiprocessing import Pool
import matplotlib.pyplot as plt
import graph_tool.all as gt
from FFBS.test_FFBS_quotient_graph import *
import time
import re


class Simulator:
    def __init__(self, config):
        """
        Initialize the simulation class with configuration.
        :param config: A dictionary containing file paths and simulation parameters.
        """
        self.snp_df_path = config["snp_df_path"] # '/labs/Aguiar/pHapCompass/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv'
        self.input_vcf_path = config["input_vcf_path"] # '/labs/Aguiar/pHapCompass/simulated_data_NEW/hapref_chr21_filtered.vcf.bgz'
        self.contig_fasta = config["contig_fasta"] # '/labs/Aguiar/pHapCompass/references/EGA.GRCh37/chr21.fa'
        self.art_path = config["art_path"]
        self.main_path = config["main_path"] # '/labs/Aguiar/pHapCompass/simulated_data_NEW'
        self.extract_hairs_path = config["extract_hairs_path"]
        self.contig_lens = config.get("contig_lens", [100])
        self.ploidies = config.get("ploidies", [3])
        self.coverages = config.get("coverages", [10])
        self.read_length = config.get("read_length", 125)
        self.mil = config.get("mean_insert_length", 400)
        self.sil = config.get("std_insert_length", 50)
        self.n_samples = config.get("n_samples", 100)
        self.densify_snps = config.get("densify_snps", True)
        self.target_spacing = config.get("target_spacing", 350)
        self.main_sh = ''
        # '#!/bin/bash\n#BATCH --job-name=pyalb\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH 
        # --mail-type=END\n#SBATCH --mem=20G\n#SBATCH --mail-user=marjan.hosseini@uconn.edu\n#SBATCH -o script.out\n#SBATCH -e ont.err\n\necho `hostname`'


    def get_slurm_header(self, name):
        header = f"#!/bin/bash\n#BATCH --job-name={name}\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=20G\n#SBATCH --mail-user=marjan.hosseini@uconn.edu\n#SBATCH -o {name}.out\n#SBATCH -e {name}.err\n\necho `hostname`\n\n"
        return header


    def adjust_snp_spacing(self, snp_positions, max_length=None):
        """
        Adjust SNP positions to achieve approximately the desired spacing.
        :param snp_positions: List of original SNP positions (sorted).
        :param target_spacing: Desired average spacing between SNPs.
        :param max_length: Maximum position to consider (end of the sequence).
        :return: List of adjusted SNP positions.
        """
        adjusted_positions = []
        prev_position = snp_positions[0]
        adjusted_positions.append(prev_position)

        for pos in snp_positions[1:]:
            while prev_position + self.target_spacing < pos:
                prev_position += self.target_spacing
                adjusted_positions.append(prev_position)
            
            if pos != prev_position:  # Avoid duplicating SNPs
                adjusted_positions.append(pos)
                prev_position = pos

        # Add SNPs if needed up to max_length
        if max_length:
            while prev_position + self.target_spacing < max_length:
                prev_position += self.target_spacing
                adjusted_positions.append(prev_position)

        return adjusted_positions


    def old_generate_genomes_fasta(self):
        """
        Generate genomes and save to FASTA format for various configurations.
        """
        snp_df = pd.read_csv(self.snp_df_path)
        # snp_df['diff'] = snp_df['POS'].diff()
        snp_positions = list(snp_df['POS'].values)

        if self.densify_snps:
            adjusted_pos = self.adjust_snp_spacing(snp_positions, max_length=None)

        with open(self.contig_fasta) as f:
            contig_name = f.readline().strip()[1:]  # Get contig name without '>'
            contig_seq = f.read().replace("\n", "")

        for contig_len in self.contig_lens:
            contig_path = os.path.join(self.main_path, f'contig_{contig_len}')
            if not os.path.exists(contig_path):
                os.makedirs(contig_path, exist_ok=True)

            for ploidy in self.ploidies:
                ploidy_path = os.path.join(contig_path, f'ploidy_{ploidy}')
                if not os.path.exists(ploidy_path):
                    os.makedirs(ploidy_path, exist_ok=True)
                # stop
                snp_counter = 0
                
                # last_snp = snp_positions[:contig_len][-1]
                if self.densify_snps:
                    last_snp = adjusted_pos[:contig_len][-1]
                else:
                    last_snp = snp_positions[:contig_len][-1]

                haplotypes = [[] for _ in range(ploidy)]
                genomes = [list(contig_seq[:last_snp]) for _ in range(ploidy)]

                vcf_in = pysam.VariantFile(self.input_vcf_path)
                original_snps = {record.pos - 1: record for record in vcf_in.fetch()}  # Map positions to records

                # with pysam.VariantFile(self.input_vcf_path) as vcf_in:
                new_header = pysam.VariantHeader()
                for line in vcf_in.header.records:
                    new_header.add_line(str(line))
                for k in range(ploidy):
                    new_header.contigs.add(f"haplotype_{k + 1}", length=len(contig_seq))
                
                new_header.add_sample("NA12878")

                output_vcf = os.path.join(ploidy_path, f'maf0.01_hapref_chr21_filtered_NA12878_ploidy{ploidy}_contig{contig_len}.vcf.gz')
                vcf_out = pysam.VariantFile(output_vcf, "wz", header=new_header)
                
                if self.densify_snps:
                # Process adjusted SNP positions
                    for pos in adjusted_pos[:contig_len]:
                        pos = int(pos)
                        ref = contig_seq[pos]
                        if pos in original_snps:
                            # Use existing record
                            record = original_snps[pos]
                            alts = record.alts[:1]
                        else:
                            # Generate new SNP record
                            ref_upper = ref.upper()
                            bases = ["A", "C", "G", "T"]
                            bases.remove(ref_upper)
                            alts = [random.choice(bases)]
                        # print(snp_counter)
                        snp_counter += 1
                        g = random.randint(1, ploidy - 1)
                        selected_genomes = random.sample(range(ploidy), g)
                        genotype = tuple(1 if i in selected_genomes else 0 for i in range(ploidy))

                        for i in range(ploidy):
                            genomes[i][pos-1] = alts[0] if genotype[i] == 1 else ref
                            haplotypes[i].append(genotype[i])

                        for k in range(ploidy):
                            new_record = vcf_out.new_record(
                                contig=f"haplotype_{k + 1}",
                                start=pos,
                                stop=pos + 1,
                                alleles=(ref, alts[0]),
                                qual=75,
                                filter=["PASS"],
                                info={}
                            )
                            new_record.samples["NA12878"]["GT"] = genotype
                            new_record.samples["NA12878"].phased = True
                            vcf_out.write(new_record)

                    vcf_out.close()

                else:

                    for record in vcf_in.fetch():
                        pos = record.pos - 1
                        ref = record.ref
                        alts = record.alts[:1]
                        # if len(alts) > 0:
                        #     alts = record.alts[:1]
                        # else:
                        #     ref_upper = ref.upper()  # Convert ref to uppercase
                        #     bases = ["A", "C", "G", "T"]
                        #     bases.remove(ref_upper)  # Remove the ref base from the options
                        #     alts = random.choice(bases)
                        if len(alts) > 0 and pos + 1 in snp_positions and snp_counter < contig_len:
                        # if pos + 1 in adjusted_pos and snp_counter < contig_len:
                            snp_counter += 1
                            g = random.randint(1, ploidy - 1)
                            selected_genomes = random.sample(range(ploidy), g)
                            genotype = tuple(1 if i in selected_genomes else 0 for i in range(ploidy))

                            for i in range(ploidy):
                                genomes[i][pos-1] = alts[0] if genotype[i] == 1 else ref
                                haplotypes[i].append(genotype[i])
                            

                            for k in range(ploidy):
                                new_record = vcf_out.new_record(
                                    contig=f"haplotype_{k + 1}",
                                    start=record.start,
                                    stop=record.stop,
                                    alleles=(ref, alts[0]),
                                    qual=75,
                                    filter=record.filter.keys(),
                                    info=record.info,
                                )
                                new_record.samples["NA12878"]["GT"] = genotype
                                new_record.samples["NA12878"].phased = True
                                vcf_out.write(new_record)

                    vcf_out.close()

                output_fasta = os.path.join(ploidy_path, f'contig_{contig_len}_ploidy_{ploidy}.fa')
                with open(output_fasta, "w") as f_out:
                    for i in range(ploidy):
                        f_out.write(f">haplotype_{i + 1}\n")
                        f_out.write("".join(genomes[i]) + "\n")

                haplotype_df = pd.DataFrame(haplotypes).T
                haplotype_df.columns = [f"haplotype_{i + 1}" for i in range(ploidy)]
                haplotype_df.to_csv(os.path.join(ploidy_path, 'haplotypes.csv'), index=False)
                print(f"Generated genomes for contig length {contig_len} and ploidy {ploidy}.")


    def generate_genomes_fasta(self):
        """
        Generate genomes and save to FASTA format for various configurations.
        """
        # snp_df = pd.read_csv(self.snp_df_path)
        # # snp_df['diff'] = snp_df['POS'].diff()
        # snp_positions = list(snp_df['POS'].values)

        # if self.densify_snps:
        #     adjusted_pos = self.adjust_snp_spacing(snp_positions, max_length=None)
        vcf_in = pysam.VariantFile(self.input_vcf_path)
        snp_positions = [record.pos for record in vcf_in.fetch()]
        snp_positions = list(np.unique(snp_positions))
        # snp_differences = np.diff(snp_positions)

        # def detect_outliers_iqr(data):
        #     """
        #     Detects outliers in a NumPy array using the IQR method.

        #     Parameters:
        #         data (np.ndarray): The input data array.

        #     Returns:
        #         np.ndarray: A boolean array where True indicates an outlier.
        #     """
        #     # Ensure the data is a NumPy array
        #     data = np.asarray(data)

        #     # Calculate the first (Q1) and third (Q3) quartiles
        #     Q1 = np.percentile(data, 25)
        #     Q3 = np.percentile(data, 75)
        #     IQR = Q3 - Q1

        #     # Define the outlier boundaries
        #     lower_bound = Q1 - 1.5 * IQR
        #     upper_bound = Q3 + 1.5 * IQR

        #     # Identify outliers
        #     outliers = (data < lower_bound) | (data > upper_bound)

        #     return outliers


        # def detect_large_differences(data, target_diff=2000, tolerance=500):
        #     """
        #     Detects indices where the difference between consecutive numbers 
        #     is higher than the target_diff plus tolerance.

        #     Parameters:
        #         data (np.ndarray): The input data array.
        #         target_diff (float): The target difference between consecutive numbers.
        #         tolerance (float): The acceptable deviation for detection.

        #     Returns:
        #         np.ndarray: Indices where the difference is higher than target_diff + tolerance.
        #     """
        #     # Ensure the data is a NumPy array
        #     data = np.asarray(data)
            
        #     # Calculate differences between consecutive numbers
        #     differences = np.diff(data)
            
        #     # Identify indices where the difference is higher than the threshold
        #     high_difference_indices = np.where(differences > target_diff + tolerance)[0]
            
        #     return high_difference_indices

        # plt.plot(range(100), snp_differences[0:100])
        # plt.xlim(0, 100)
        # plt.show()

        # devis = []
        # for i in range(len(snp_positions) - 100):
        #     # np.mean(snp_differences[i: i+100])
        #     data = snp_positions[i: i+100]
        #     # outliers = detect_outliers_iqr(data)
        #     indices, differences = detect_large_differences(data, target_diff=3000, tolerance=500)
            
        #     deviation_values = [(data[idx], data[idx + 1]) for idx in indices]
        #     # print(i, len(deviation_values))
        #     devis.append(len(deviation_values))
        #     # print('----------------------------------')

        # np.min(devis)

        with open(self.contig_fasta) as f:
            contig_name = f.readline().strip()[1:]  # Get contig name without '>'
            contig_seq = f.read().replace("\n", "")

        for contig_len in self.contig_lens:
            contig_path = os.path.join(self.main_path, f'contig_{contig_len}')
            if not os.path.exists(contig_path):
                os.makedirs(contig_path, exist_ok=True)

            for ploidy in self.ploidies:
                ploidy_path = os.path.join(contig_path, f'ploidy_{ploidy}')
                if not os.path.exists(ploidy_path):
                    os.makedirs(ploidy_path, exist_ok=True)
                # stop
                
                last_snp = snp_positions[-1]
                # if self.densify_snps:
                #     last_snp = adjusted_pos[:contig_len][-1]
                # else:
                # last_snp = snp_positions[:contig_len][-1]

                haplotypes = [[] for _ in range(ploidy)]
                genomes = [list(contig_seq[:last_snp]) for _ in range(ploidy)]
                # genomes = [list(contig_seq) for _ in range(ploidy)]

                # vcf_in = pysam.VariantFile(self.input_vcf_path)
                # original_snps = {record.pos - 1: record for record in vcf_in.fetch()}  # Map positions to records

                # with pysam.VariantFile(self.input_vcf_path) as vcf_in:
                new_header = pysam.VariantHeader()
                for line in vcf_in.header.records:
                    new_header.add_line(str(line))
                # for k in range(ploidy):
                #     new_header.contigs.add(f"haplotype_{k + 1}", length=len(contig_seq))
                # new_header.contigs.add(f"contig_{contig_len}", length=len(contig_seq))
                # new_header.contigs.add('AHIQ01000001.1', length=len(contig_seq))
                
                new_header.add_sample("NA12878")

                output_vcf = os.path.join(ploidy_path, f'NA12878_ploidy{ploidy}_contig{contig_len}.vcf.gz')
                vcf_out = pysam.VariantFile(output_vcf, "wz", header=new_header)
                snp_counter = 0
                current_snp = 0

                for record in vcf_in.fetch():

                    pos = record.pos - 1
                    
                    ref = record.ref
                    alts = record.alts[:1]
                    # print(pos, record.start, record.stop, alts, ref, contig_seq[pos], len(alts[0]))
                    # print(ref, contig_seq[pos])
                    if len(ref) > 1 or 'N' in ref or 'n' in ref:
                        print('oomadam too 1')
                        ref = random.choice(['A', 'C', 'G', 'T'])
                    
                    if len(alts[0]) > 1 or len(alts) == 0 or 'N' in alts[0] or 'n' in alts[0] or alts[0] == ref:
                        print('oomadam too 2')
                        ref_upper = ref.upper()
                        alternatives = ['A', 'C', 'G', 'T']
                        alternatives.remove(ref_upper)
                        alts = [random.choice(alternatives)]
                    # if len(alts) > 0:
                    #     alts = record.alts[:1]
                    # else:
                    #     ref_upper = ref.upper()  # Convert ref to uppercase
                    #     bases = ["A", "C", "G", "T"]
                    #     bases.remove(ref_upper)  # Remove the ref base from the options
                    #     alts = random.choice(bases)
                    
                    # if len(alts) > 0 and snp_counter < contig_len and len(ref) == 1 and len(alts[0]) == 1 and ref != 'N' and pos != current_snp:
                    if len(alts) > 0 and len(ref) == 1 and len(alts[0]) == 1 and ref != 'N' and pos != current_snp:
                    # if pos + 1 in adjusted_pos and snp_counter < contig_len:
                        # print(pos, ref, genomes[0][pos], alts)
                        # print(pos, record.start, record.stop, alts, ref, contig_seq[pos], len(alts[0]))
                        print(pos, ref, contig_seq[pos], alts)
                        processed_info = {key: (value[0],) if isinstance(value, tuple) and len(value) > 1 else value for key, value in record.info.items()}
                        snp_counter += 1
                        current_snp = pos
                        g = random.randint(1, ploidy - 1)
                        selected_genomes = random.sample(range(ploidy), g)
                        genotype = tuple(1 if i in selected_genomes else 0 for i in range(ploidy))

                        for i in range(ploidy):
                            genomes[i][pos] = alts[0] if genotype[i] == 1 else ref
                            haplotypes[i].append(genotype[i])
                        
                        new_record = vcf_out.new_record(
                            # pos=record.pos,
                            # contig=f"contig_{contig_len}",
                            contig='chr21',
                            start=record.start,
                            stop=record.stop,
                            alleles=(ref, alts[0]),
                            qual=75,
                            filter=record.filter.keys(),
                            # info=record.info,
                            info=processed_info)
                        new_record.samples["NA12878"]["GT"] = genotype
                        new_record.samples["NA12878"].phased = True
                        vcf_out.write(new_record)

                    # else:
                    #     break

                vcf_out.close()

                output_fasta = os.path.join(ploidy_path, f'contig_{contig_len}_ploidy_{ploidy}.fa')
                with open(output_fasta, "w") as f_out:
                    for i in range(ploidy):
                        f_out.write(f">haplotype_{i + 1}\n")
                        f_out.write("".join(genomes[i]) + "\n")

                haplotype_df = pd.DataFrame(haplotypes).T
                haplotype_df.columns = [f"haplotype_{i + 1}" for i in range(ploidy)]
                haplotype_df.to_csv(os.path.join(ploidy_path, 'haplotypes.csv'), index=False)
                print(f"Generated genomes for contig length {contig_len} and ploidy {ploidy}.")


    def simulate_fastq_art(self):
        """
        Simulate FASTQ files using ART.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                fasta_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'contig_{contig_len}_ploidy_{ploidy}.fa')
                this_sh_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'01_simulate_{contig_len}_{ploidy}.sh')

                # to_print = self.get_slurm_header('fastq')
                to_print = ''
                for coverage in self.coverages:
                    cov_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'cov_{coverage}')
                    if not os.path.exists(cov_path):
                        os.makedirs(cov_path, exist_ok=True)
                    fastq_path = os.path.join(cov_path, 'fastq')
                    if not os.path.exists(fastq_path):
                        os.makedirs(fastq_path, exist_ok=True)

                    for rd in range(self.n_samples):
                        rs = random.randint(1, 2**32)
                        command = f'{self.art_path} -ss HS25 -rs {rs} -i {fasta_path} -p -na -l {self.read_length} -f {coverage} -m {self.mil} -s {self.sil} -o {fastq_path}/{str(rd).zfill(2)}\n'
                        to_print += command

                with open(this_sh_path, 'w') as f:
                    f.write(to_print)
                self.main_sh += f'sh {this_sh_path}\n'


    def align_fastq_files(self):
        """
        Align simulated FASTQ files using BWA.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                # fasta_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'contig_{}_ploidy_{}.fa'.format(contig_len, ploidy))
                # fa_index = 'bwa index {}\n\n'.format(fasta_path)
                this_sh_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '02_align_{}_{}.sh'.format(contig_len, ploidy))
                # to_print = self.get_slurm_header('align')
                # to_print += 'module load bwa-mem2/2.1\nmodule load bwa/0.7.17\n' + fa_index
                # to_print = fa_index
                to_print = ''
                for coverage in self.coverages:
                    this_cov_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                    fastq_path = os.path.join(this_cov_path, 'fastq')
                    bam_path = os.path.join(this_cov_path, 'bam')
                    if not os.path.exists(bam_path):
                        os.makedirs(bam_path)
                    for rd in range(self.n_samples):
                        command = 'bwa mem {} {}/{}1.fq {}/{}2.fq > {}/{}.sam\n'.format(self.contig_fasta, fastq_path, str(rd).zfill(2), fastq_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
                        sort_com = 'samtools view -Sb {}/{}.sam | samtools sort -o {}/{}.bam\n'.format(bam_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
                        index_com = 'samtools index {}/{}.bam\n'.format(bam_path, str(rd).zfill(2))
                        rmv_com = 'rm {}/{}.sam\n\n'.format(bam_path, str(rd).zfill(2))
                        to_print += command
                        to_print += sort_com
                        to_print += index_com
                        to_print += rmv_com
                with open(this_sh_path, 'w') as f:
                    f.write(to_print)
                self.main_sh += f'sh {this_sh_path}\n'


    def extract_hairs(self):
        """
        Extract haplotypes using ExtractHAIRS.
        """
        vcf_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/updated_NA12878_extracted.vcf'
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                this_sh_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '03_extract_hairs_{}_{}.sh'.format(contig_len, ploidy))
                # to_print = self.get_slurm_header('exttract_hairs')
                # to_print += 'module load bcftools/1.20\nmodule load htslib/1.20\n\n'
                to_print = ''
                # this_vcf_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'AWRI_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
                # this_vcf_sorted_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'sorted_AWRI_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
                # unzipped_vcf = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'sorted_AWRI_ploidy{}_contig{}.vcf'.format(ploidy, contig_len))
                # sort_cmd = 'bcftools sort {} -Oz -o {}\n'.format(this_vcf_path, this_vcf_sorted_path)
                # idx_cmd = 'tabix -p vcf {}\n'.format(this_vcf_sorted_path)
                # unzip_cmd = 'gunzip {}\n\n'.format(this_vcf_sorted_path)
                # to_print += sort_cmd
                # to_print += idx_cmd
                # to_print += unzip_cmd
                for coverage in self.coverages:
                    this_cov_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                    bam_path = os.path.join(this_cov_path, 'bam')
                    frag_path = os.path.join(this_cov_path, 'frag')
                    if not os.path.exists(frag_path):
                        os.makedirs(frag_path)
                    for rd in range(self.n_samples): 
                        command = '{} --bam {}/{}.bam --vcf {} --out {}/{}.frag\n'.format(self.extract_hairs_path, bam_path, str(rd).zfill(2), vcf_path, frag_path, str(rd).zfill(2))
                        to_print += command
                with open(this_sh_path, 'w') as f:
                    f.write(to_print)
                self.main_sh += f'sh {this_sh_path}\n'


    def simulate(self):
        """
        Run the entire simulation pipeline.
        """
        self.generate_genomes_fasta()
        self.simulate_fastq_art()
        self.align_fastq_files()
        self.extract_hairs()
        with open(os.path.join(self.main_path, 'main.sh'), 'w') as f:
            f.write(self.main_sh)


class SimulatorAWRI:
    def __init__(self, config):
        """
        Initialize the simulation class with configuration.
        :param config: A dictionary containing file paths and simulation parameters.
        """
        # self.snp_df_path = config["snp_df_path"] # '/labs/Aguiar/pHapCompass/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv'
        self.input_vcf_path = config["input_vcf_path"] # '/labs/Aguiar/pHapCompass/simulated_data_NEW/hapref_chr21_filtered.vcf.bgz'
        self.contig_fasta = config["contig_fasta"] # '/labs/Aguiar/pHapCompass/references/EGA.GRCh37/chr21.fa'
        self.art_path = config["art_path"]
        self.main_path = config["main_path"] # '/labs/Aguiar/pHapCompass/simulated_data_NEW'
        self.extract_hairs_path = config["extract_hairs_path"]
        self.contig_lens = config.get("contig_lens", [100])
        self.ploidies = config.get("ploidies", [3])
        self.coverages = config.get("coverages", [10])
        self.read_length = config.get("read_length", 125)
        self.mil = config.get("mean_insert_length", 400)
        self.sil = config.get("std_insert_length", 50)
        self.n_samples = config.get("n_samples", 100)
        self.densify_snps = config.get("densify_snps", True)
        self.target_spacing = config.get("target_spacing", 350)
        self.sbatch_str = config.get("sbatch_str", '')
        self.main_sh = self.sbatch_str + '\n\n' + 'echo `hostname`\n\nmodule load bwa/0.7.17\nmodule load bwa-mem2/2.0 bwa-mem2/2.1\n\n\n'
        # '#!/bin/bash\n#BATCH --job-name=pyalb\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH 
        # --mail-type=END\n#SBATCH --mem=20G\n#SBATCH --mail-user=marjan.hosseini@uconn.edu\n#SBATCH -o script.out\n#SBATCH -e ont.err\n\necho `hostname`'

    def generate_genomes_fasta(self):
        """
        Generate genomes and save to FASTA format for various configurations.
        """
        # snp_df = pd.read_csv(self.snp_df_path)
        # # snp_df['diff'] = snp_df['POS'].diff()
        # snp_positions = list(snp_df['POS'].values)

        # if self.densify_snps:
        #     adjusted_pos = self.adjust_snp_spacing(snp_positions, max_length=None)
        vcf_in = pysam.VariantFile(self.input_vcf_path)
        snp_positions = [record.pos for record in vcf_in.fetch()]
        # snp_positions = np.unique(snp_positions)

        with open(self.contig_fasta) as f:
            contig_name = f.readline().strip()[1:]  # Get contig name without '>'
            contig_seq = f.read().replace("\n", "")

        for contig_len in self.contig_lens:
            contig_path = os.path.join(self.main_path, f'contig_{contig_len}')
            if not os.path.exists(contig_path):
                os.makedirs(contig_path, exist_ok=True)

            for ploidy in self.ploidies:
                ploidy_path = os.path.join(contig_path, f'ploidy_{ploidy}')
                if not os.path.exists(ploidy_path):
                    os.makedirs(ploidy_path, exist_ok=True)
                # stop
                
                last_snp = snp_positions[:contig_len][-1]
                # if self.densify_snps:
                #     last_snp = adjusted_pos[:contig_len][-1]
                # else:
                # last_snp = snp_positions[:contig_len][-1]

                haplotypes = [[] for _ in range(ploidy)]
                genomes = [list(contig_seq[:last_snp]) for _ in range(ploidy)]
                # genomes = [list(contig_seq) for _ in range(ploidy)]

                # vcf_in = pysam.VariantFile(self.input_vcf_path)
                # original_snps = {record.pos - 1: record for record in vcf_in.fetch()}  # Map positions to records

                # with pysam.VariantFile(self.input_vcf_path) as vcf_in:
                new_header = pysam.VariantHeader()
                for line in vcf_in.header.records:
                    new_header.add_line(str(line))
                # for k in range(ploidy):
                #     new_header.contigs.add(f"haplotype_{k + 1}", length=len(contig_seq))
                # new_header.contigs.add(f"contig_{contig_len}", length=len(contig_seq))
                # new_header.contigs.add('AHIQ01000001.1', length=len(contig_seq))
                
                new_header.add_sample("AWRI")

                output_vcf = os.path.join(ploidy_path, f'AWRI_ploidy{ploidy}_contig{contig_len}.vcf.gz')
                vcf_out = pysam.VariantFile(output_vcf, "wz", header=new_header)
                snp_counter = 0
                # current_snp = 0
                for record in vcf_in.fetch():

                    pos = record.pos - 1
                    
                    ref = record.ref
                    alts = record.alts[:1]
                    # print(pos, record.start, record.stop, alts, ref, contig_seq[pos], len(alts[0]))
                    # print(ref, contig_seq[pos])
                    if len(ref) > 1 or 'N' in ref or 'n' in ref:
                        ref = random.choice(['A', 'C', 'G', 'T'])
                    
                    if len(alts[0]) > 1 or len(alts) == 0 or 'N' in alts[0] or 'n' in alts[0] or alts[0] == ref:
                        ref_upper = ref.upper()
                        alternatives = ['A', 'C', 'G', 'T']
                        alternatives.remove(ref_upper)
                        alts = [random.choice(alternatives)]
                    # if len(alts) > 0:
                    #     alts = record.alts[:1]
                    # else:
                    #     ref_upper = ref.upper()  # Convert ref to uppercase
                    #     bases = ["A", "C", "G", "T"]
                    #     bases.remove(ref_upper)  # Remove the ref base from the options
                    #     alts = random.choice(bases)
                    
                    if len(alts) > 0 and snp_counter < contig_len and len(ref) == 1 and len(alts[0]) == 1 and ref != 'N':
                    # if pos + 1 in adjusted_pos and snp_counter < contig_len:
                        # print(pos, ref, genomes[0][pos], alts)
                        # print(pos, record.start, record.stop, alts, ref, contig_seq[pos], len(alts[0]))
                        print(pos, ref, contig_seq[pos], alts)
                        processed_info = {key: (value[0],) if isinstance(value, tuple) and len(value) > 1 else value for key, value in record.info.items()}
                        snp_counter += 1
                        # current_snp = pos
                        g = random.randint(1, ploidy - 1)
                        selected_genomes = random.sample(range(ploidy), g)
                        genotype = tuple(1 if i in selected_genomes else 0 for i in range(ploidy))

                        for i in range(ploidy):
                            genomes[i][pos] = alts[0] if genotype[i] == 1 else ref
                            haplotypes[i].append(genotype[i])
                        
                        # for k in range(ploidy):
                        #     new_record = vcf_out.new_record(
                        #         # pos=record.pos,
                        #         contig=f"haplotype_{k + 1}",
                        #         start=record.start,
                        #         stop=record.stop,
                        #         alleles=(ref, alts[0]),
                        #         qual=75,
                        #         filter=record.filter.keys(),
                        #         # info=record.info,
                        #         info=processed_info
                        #     )
                        #     new_record.samples["AWRI"]["GT"] = genotype
                        #     new_record.samples["AWRI"].phased = True
                        #     vcf_out.write(new_record)

                        new_record = vcf_out.new_record(
                            # pos=record.pos,
                            # contig=f"contig_{contig_len}",
                            contig='AHIQ01000001.1',
                            start=record.start,
                            stop=record.stop,
                            alleles=(ref, alts[0]),
                            qual=75,
                            filter=record.filter.keys(),
                            # info=record.info,
                            info=processed_info)
                        new_record.samples["AWRI"]["GT"] = genotype
                        new_record.samples["AWRI"].phased = True
                        vcf_out.write(new_record)

                vcf_out.close()

                output_fasta = os.path.join(ploidy_path, f'contig_{contig_len}_ploidy_{ploidy}.fa')
                with open(output_fasta, "w") as f_out:
                    for i in range(ploidy):
                        f_out.write(f">haplotype_{i + 1}\n")
                        f_out.write("".join(genomes[i]) + "\n")

                haplotype_df = pd.DataFrame(haplotypes).T
                haplotype_df.columns = [f"haplotype_{i + 1}" for i in range(ploidy)]
                haplotype_df.to_csv(os.path.join(ploidy_path, 'haplotypes.csv'), index=False)
                print(f"Generated genomes for contig length {contig_len} and ploidy {ploidy}.")


    def simulate_fastq_art(self):
        """
        Simulate FASTQ files using ART.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                fasta_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'contig_{contig_len}_ploidy_{ploidy}.fa')
                this_sh_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'01_simulate_{contig_len}_{ploidy}.sh')

                # to_print = self.get_slurm_header('fastq')
                to_print = ''
                for coverage in self.coverages:
                    cov_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'cov_{coverage}')
                    if not os.path.exists(cov_path):
                        os.makedirs(cov_path, exist_ok=True)
                    fastq_path = os.path.join(cov_path, 'fastq')
                    if not os.path.exists(fastq_path):
                        os.makedirs(fastq_path, exist_ok=True)

                    for rd in range(self.n_samples):
                        rs = random.randint(1, 2**32)
                        command = f'{self.art_path} -ss HS25 -rs {rs} -i {fasta_path} -p -na -l {self.read_length} -f {coverage} -m {self.mil} -s {self.sil} -o {fastq_path}/{str(rd).zfill(2)}\n'
                        to_print += command

                with open(this_sh_path, 'w') as f:
                    f.write(to_print)
                self.main_sh += f'sh {this_sh_path}\n'


    def align_fastq_files(self):
        """
        Align simulated FASTQ files using BWA.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                # fasta_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'contig_{}_ploidy_{}.fa'.format(contig_len, ploidy))
                # fa_index = 'bwa index {}\n\n'.format(fasta_path)
                this_sh_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '02_align_{}_{}.sh'.format(contig_len, ploidy))
                # to_print = self.get_slurm_header('align')
                # to_print += 'module load bwa-mem2/2.1\nmodule load bwa/0.7.17\n' + fa_index
                # to_print = fa_index
                to_print = ''
                for coverage in self.coverages:
                    this_cov_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                    fastq_path = os.path.join(this_cov_path, 'fastq')
                    bam_path = os.path.join(this_cov_path, 'bam')
                    if not os.path.exists(bam_path):
                        os.makedirs(bam_path)
                    for rd in range(self.n_samples):
                        command = 'bwa mem {} {}/{}1.fq {}/{}2.fq > {}/{}.sam\n'.format(self.contig_fasta, fastq_path, str(rd).zfill(2), fastq_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
                        sort_com = 'samtools view -Sb {}/{}.sam | samtools sort -o {}/{}.bam\n'.format(bam_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
                        index_com = 'samtools index {}/{}.bam\n'.format(bam_path, str(rd).zfill(2))
                        rmv_com = 'rm {}/{}.sam\n\n'.format(bam_path, str(rd).zfill(2))
                        to_print += command
                        to_print += sort_com
                        to_print += index_com
                        to_print += rmv_com
                with open(this_sh_path, 'w') as f:
                    f.write(to_print)
                self.main_sh += f'sh {this_sh_path}\n'


    def extract_hairs(self):
        """
        Extract haplotypes using ExtractHAIRS.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                this_sh_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '03_extract_hairs_{}_{}.sh'.format(contig_len, ploidy))
                # to_print = self.get_slurm_header('exttract_hairs')
                # to_print += 'module load bcftools/1.20\nmodule load htslib/1.20\n\n'
                to_print = ''
                # this_vcf_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'AWRI_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
                # this_vcf_sorted_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'sorted_AWRI_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
                # unzipped_vcf = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'sorted_AWRI_ploidy{}_contig{}.vcf'.format(ploidy, contig_len))
                # sort_cmd = 'bcftools sort {} -Oz -o {}\n'.format(this_vcf_path, this_vcf_sorted_path)
                # idx_cmd = 'tabix -p vcf {}\n'.format(this_vcf_sorted_path)
                # unzip_cmd = 'gunzip {}\n\n'.format(this_vcf_sorted_path)
                # to_print += sort_cmd
                # to_print += idx_cmd
                # to_print += unzip_cmd
                for coverage in self.coverages:
                    this_cov_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                    bam_path = os.path.join(this_cov_path, 'bam')
                    frag_path = os.path.join(this_cov_path, 'frag')
                    if not os.path.exists(frag_path):
                        os.makedirs(frag_path)
                    for rd in range(self.n_samples): 
                        command = '{} --bam {}/{}.bam --vcf {} --out {}/{}.frag\n'.format(self.extract_hairs_path, bam_path, str(rd).zfill(2), self.input_vcf_path, frag_path, str(rd).zfill(2))
                        to_print += command
                with open(this_sh_path, 'w') as f:
                    f.write(to_print)
                self.main_sh += f'sh {this_sh_path}\n'


    def simulate(self):
        """
        Run the entire simulation pipeline.
        """
        self.generate_genomes_fasta()
        self.simulate_fastq_art()
        self.align_fastq_files()
        self.extract_hairs()
        print('save main.sh')
        with open(os.path.join(self.main_path, 'main.sh'), 'w') as f:
            f.write(self.main_sh)


class SimulatorNA12878:
    def __init__(self, config):
        """
        Initialize the simulation class with configuration.
        :param config: A dictionary containing file paths and simulation parameters.
        """
        # self.snp_df_path = config["snp_df_path"] # '/labs/Aguiar/pHapCompass/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv'
        self.input_vcf_path = config["input_vcf_path"] # '/labs/Aguiar/pHapCompass/simulated_data_NEW/hapref_chr21_filtered.vcf.bgz'
        self.contig_fasta = config["contig_fasta"] # '/labs/Aguiar/pHapCompass/references/EGA.GRCh37/chr21.fa'
        self.art_path = config["art_path"]
        self.main_path = config["main_path"] # '/labs/Aguiar/pHapCompass/simulated_data_NEW'
        self.extract_hairs_path = config["extract_hairs_path"]
        self.contig_lens = config.get("contig_lens", [100])
        self.ploidies = config.get("ploidies", [3])
        self.coverages = config.get("coverages", [10])
        self.read_length = config.get("read_length", 125)
        self.mil = config.get("mean_insert_length", 400)
        self.sil = config.get("std_insert_length", 50)
        self.n_samples = config.get("n_samples", 100)
        self.densify_snps = config.get("densify_snps", True)
        self.target_spacing = config.get("target_spacing", 350)
        self.sbatch_str = config.get("sbatch_str", '')
        self.main_sh = self.sbatch_str + '\n\n' + 'echo `hostname`\n\nmodule load bwa/0.7.17\nmodule load bwa-mem2/2.0 bwa-mem2/2.1\n\n\n'
        # '#!/bin/bash\n#BATCH --job-name=pyalb\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH 
        # --mail-type=END\n#SBATCH --mem=20G\n#SBATCH --mail-user=marjan.hosseini@uconn.edu\n#SBATCH -o script.out\n#SBATCH -e ont.err\n\necho `hostname`'


    def generate_genomes_fasta(self):
        """
        Generate genomes and save to FASTA format for various configurations.
        """
        # snp_df = pd.read_csv(self.snp_df_path)
        # # snp_df['diff'] = snp_df['POS'].diff()
        # snp_positions = list(snp_df['POS'].values)

        # if self.densify_snps:
        #     adjusted_pos = self.adjust_snp_spacing(snp_positions, max_length=None)
        vcf_in = pysam.VariantFile(self.input_vcf_path)
        snp_positions = [record.pos for record in vcf_in.fetch()]
        # snp_positions = np.unique(snp_positions)

        with open(self.contig_fasta) as f:
            contig_name = f.readline().strip()[1:]  # Get contig name without '>'
            contig_seq = f.read().replace("\n", "")

        for contig_len in self.contig_lens:
            contig_path = os.path.join(self.main_path, f'contig_{contig_len}')
            if not os.path.exists(contig_path):
                os.makedirs(contig_path, exist_ok=True)

            for ploidy in self.ploidies:
                ploidy_path = os.path.join(contig_path, f'ploidy_{ploidy}')
                if not os.path.exists(ploidy_path):
                    os.makedirs(ploidy_path, exist_ok=True)
                # stop
                
                last_snp = snp_positions[:contig_len][-1]

                genomes = [list(contig_seq[:last_snp]) for _ in range(ploidy)]
                # genomes = [list(contig_seq) for _ in range(ploidy)]
                new_fasta = list(contig_seq[:last_snp])
                haplotype_df = pd.DataFrame(columns=[f"haplotype_{i + 1}" for i in range(ploidy)], index=range(contig_len))
                # vcf_in = pysam.VariantFile(self.input_vcf_path)
                # original_snps = {record.pos - 1: record for record in vcf_in.fetch()}  # Map positions to records

                # with pysam.VariantFile(self.input_vcf_path) as vcf_in:
                new_header = pysam.VariantHeader()
                for line in vcf_in.header.records:
                    # print(line)
                    if not 'contig' in str(line):
                        # print(line)
                        new_header.add_line(str(line))

                new_header.add_line('##contig=<ID=chr21_{},length={}>'.format(contig_len, len(contig_seq[:last_snp])))
                new_header.add_sample("NA12878")

                output_vcf = os.path.join(ploidy_path, f'NA12878_ploidy{ploidy}_contig{contig_len}.vcf.gz')
                vcf_out = pysam.VariantFile(output_vcf, "w", header=new_header)
                snp_counter = 0
                # current_snp = 0
                for record in vcf_in.fetch():
                    
                    pos = record.pos - 1
                    
                    ref = record.ref
                    alts = record.alts[:1]
                    # print(pos, record.start, record.stop, alts, ref, contig_seq[pos], len(alts[0]))
                    # print(ref, contig_seq[pos])
                    # if len(ref) > 1 or 'N' in ref or 'n' in ref:
                    if len(ref) > 1 or ref not in ['A', 'C', 'G', 'T']:
                        ref = random.choice(['A', 'C', 'G', 'T'])
                    
                    # if len(alts[0]) > 1 or len(alts) == 0 or 'N' in alts[0] or 'n' in alts[0] or alts[0] == ref:
                    if len(alts[0]) > 1 or len(alts) == 0 or alts[0] not in ['A', 'C', 'G', 'T'] or alts[0] == ref:
                        ref_upper = ref.upper()
                        alternatives = ['A', 'C', 'G', 'T']
                        alternatives.remove(ref_upper)
                        alts = [random.choice(alternatives)]

                    if len(alts) > 0 and snp_counter < contig_len and len(ref) == 1 and len(alts[0]) == 1 and ref != 'N':
                    # if pos + 1 in adjusted_pos and snp_counter < contig_len:
                        # print(pos, ref, genomes[0][pos], alts)
                        # print(pos, record.start, record.stop, alts, ref, contig_seq[pos], len(alts[0]))
                        # print(pos, ref, contig_seq[pos], alts)
                        processed_info = {key: (value[0],) if isinstance(value, tuple) and len(value) > 1 else value for key, value in record.info.items()}
                        
                        # current_snp = pos
                        g = random.randint(1, ploidy - 1)
                        selected_genomes = random.sample(range(ploidy), g)
                        genotype = tuple(1 if i in selected_genomes else 0 for i in range(ploidy))

                        new_fasta[pos] = ref
                        haplotype_df.loc[snp_counter, :] = genotype
                        
                        for i in range(ploidy):
                            genomes[i][pos] = alts[0] if genotype[i] == 1 else ref
                            # haplotypes[i].append(genotype[i])
                        
                        new_record = vcf_out.new_record(
                            # pos=record.pos,
                            # contig=f"contig_{contig_len}",
                            contig='chr21_{}'.format(contig_len),
                            # pos=record.pos,
                            start=record.start,
                            stop=record.stop,
                            alleles=(ref, alts[0]),
                            qual=75,
                            filter=record.filter.keys(),
                            # info=record.info,
                            info=processed_info)
                        new_record.samples["NA12878"]["GT"] = genotype
                        new_record.samples["NA12878"].phased = False
                        vcf_out.write(new_record)
                        snp_counter += 1
                vcf_out.close()

                output_fasta = os.path.join(ploidy_path, f'contig_{contig_len}_ploidy_{ploidy}.fa')
                with open(output_fasta, "w") as f_out:
                    for i in range(ploidy):
                        f_out.write(f">haplotype_{i + 1}\n")
                        f_out.write("".join(genomes[i]) + "\n")

                ref_fasta = os.path.join(ploidy_path, f'contig_{contig_len}.fa')
                with open(ref_fasta, "w") as f_out:
                    f_out.write('>chr21_{}\n'.format(contig_len))
                    f_out.write("".join(new_fasta) + "\n")

                # haplotype_df = haplotype_df.T
                haplotype_df.columns = [f"haplotype_{i + 1}" for i in range(ploidy)]
                haplotype_df.to_csv(os.path.join(ploidy_path, 'haplotypes.csv'), index=False)

                print(f"Generated genomes for contig length {contig_len} and ploidy {ploidy}.")


    def simulate_fastq_art(self):
        """
        Simulate FASTQ files using ART.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                fasta_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'contig_{contig_len}_ploidy_{ploidy}.fa')
                this_sh_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'01_simulate_{contig_len}_{ploidy}.sh')

                # to_print = self.get_slurm_header('fastq')
                to_print = ''
                for coverage in self.coverages:
                    cov_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'cov_{coverage}')
                    if not os.path.exists(cov_path):
                        os.makedirs(cov_path, exist_ok=True)
                    fastq_path = os.path.join(cov_path, 'fastq')
                    if not os.path.exists(fastq_path):
                        os.makedirs(fastq_path, exist_ok=True)

                    for rd in range(self.n_samples):
                        rs = random.randint(1, 2**32)
                        command = f'{self.art_path} -ss HS25 -rs {rs} -i {fasta_path} -p -na -l {self.read_length} -f {coverage} -m {self.mil} -s {self.sil} -o {fastq_path}/{str(rd).zfill(2)}\n'
                        to_print += command

                with open(this_sh_path, 'w') as f:
                    f.write(to_print)
                self.main_sh += f'sh {this_sh_path}\n'


    def align_fastq_files(self):
        """
        Align simulated FASTQ files using BWA.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                ref_fasta = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'contig_{}.fa'.format(contig_len))
                # fasta_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'contig_{}_ploidy_{}.fa'.format(contig_len, ploidy))
                # fa_index = 'bwa index {}\n\n'.format(fasta_path)
                this_sh_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '02_align_{}_{}.sh'.format(contig_len, ploidy))
                # to_print = self.get_slurm_header('align')
                # to_print += 'module load bwa-mem2/2.1\nmodule load bwa/0.7.17\n' + fa_index
                # to_print = fa_index
                to_print = 'bwa index {}\n\n'.format(ref_fasta)
                for coverage in self.coverages:
                    this_cov_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                    fastq_path = os.path.join(this_cov_path, 'fastq')
                    bam_path = os.path.join(this_cov_path, 'bam')
                    if not os.path.exists(bam_path):
                        os.makedirs(bam_path)
                    for rd in range(self.n_samples):
                        command = 'bwa mem {} {}/{}1.fq {}/{}2.fq > {}/{}.sam; '.format(ref_fasta, fastq_path, str(rd).zfill(2), fastq_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
                        sort_com = 'samtools view -Sb {}/{}.sam | samtools sort -o {}/{}.bam; '.format(bam_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
                        index_com = 'samtools index {}/{}.bam; '.format(bam_path, str(rd).zfill(2))
                        rmv_com = 'rm {}/{}.sam; \n'.format(bam_path, str(rd).zfill(2))
                        to_print += command
                        to_print += sort_com
                        to_print += index_com
                        to_print += rmv_com
                with open(this_sh_path, 'w') as f:
                    f.write(to_print)
                self.main_sh += f'sh {this_sh_path}\n'


    def extract_hairs(self):
        """
        Extract haplotypes using ExtractHAIRS.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                this_sh_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '03_extract_hairs_{}_{}.sh'.format(contig_len, ploidy))

                vcf_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'NA12878_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
                to_print = 'gunzip {}\n\n'.format(vcf_path)
                vcf_path_unzipped = vcf_path.split('.gz')[0]
                for coverage in self.coverages:
                    this_cov_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                    bam_path = os.path.join(this_cov_path, 'bam')
                    frag_path = os.path.join(this_cov_path, 'frag')
                    if not os.path.exists(frag_path):
                        os.makedirs(frag_path)
                    for rd in range(self.n_samples): 
                        command = '{} --bam {}/{}.bam --vcf {} --out {}/{}.frag\n'.format(self.extract_hairs_path, bam_path, str(rd).zfill(2), vcf_path_unzipped, frag_path, str(rd).zfill(2))
                        to_print += command
                with open(this_sh_path, 'w') as f:
                    f.write(to_print)
                self.main_sh += f'sh {this_sh_path}\n'


    def simulate(self):
        """
        Run the entire simulation pipeline.
        """
        self.generate_genomes_fasta()
        self.simulate_fastq_art()
        self.align_fastq_files()
        self.extract_hairs()
        print('save main.sh')
        with open(os.path.join(self.main_path, 'main.sh'), 'w') as f:
            f.write(self.main_sh)


def remove_duplicate_positions():
    input_vcf = "/mnt/research/aguiarlab/proj/HaplOrbit/reference/chr21_NA12878/updated_NA12878_extracted.vcf"  # Change to your actual input VCF file
    output_vcf = "/mnt/research/aguiarlab/proj/HaplOrbit/reference/chr21_NA12878/updated_NA12878_extracted_no_duplicate.vcf"

    # Open input VCF file
    vcf_in = pysam.VariantFile(input_vcf, "r")

    # Create output VCF file with the same header
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    # Dictionary to store seen (chromosome, position) pairs
    seen_positions = set()

    # Iterate through VCF records
    for record in vcf_in:
        key = (record.chrom, record.pos)  # Unique identifier (Chromosome, Position)
        
        if key not in seen_positions:
            vcf_out.write(record)  # Write the first occurrence
            seen_positions.add(key)  # Mark this (chrom, pos) as seen

    # Close VCF files
    vcf_in.close()
    vcf_out.close()


def extract_positions(vcf_path):
    positions = []
    vcf_in = pysam.VariantFile(vcf_path)  # Open the VCF file
    for record in vcf_in.fetch():
        positions.append(record.pos)  # Extract the position (1-based)
    return positions


def plot_positions_distances(positions):
    # vcf_path = '/labs/Aguiar/pHapCompass/datasets/SRR10489264/variants_freebayes.vcf'
    vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/datasets/gold_standard_phasing/vcf-diploid/HG00514.chr1.vcf.gz'
    vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878/contig_100/ploidy_3/NA12878_ploidy3_contig100.vcf'
    positions = extract_positions(vcf_path)
    positions = sorted(positions)
    differences = [positions[i] - positions[i - 1] for i in range(1, len(positions))]
    plt.figure(figsize=(12, 6))
    plt.hist(differences, bins=50, color='tab:blue', alpha=0.7, edgecolor='black')
    plt.xlabel('Distances (bp)', fontdict={'size': 20})
    plt.ylabel('Frequency', fontdict={'size': 20})
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    # plt.title('Histogram of Differences Between Consecutive VCF Positions')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xscale('log')
    # plt.xlim(0, 100000)
    plt.show()
    plt.savefig('/mnt/research/aguiarlab/proj/HaplOrbit/results/position_differences.png')


def extract_column_NA12878(file_path):
    # file_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.1/windows/50000/sample_NA12878.txt'
    df = pd.read_csv(file_path, skiprows=48, sep=' ')
    hetero_df = df[df['NA12878'].isin(['1|0', '0|1'])].reset_index(drop=True)
    hetero_df = hetero_df.sort_values(by='POS').reset_index(drop=True)
    hetero_df['sim_pos'] = hetero_df.index
    return hetero_df


def extract_column_NA12878_v2(file_path):

    # file_path = "/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/windows/10/hapref_chr21_filtered.vcf.bgz_sample0_len10.bcf"
    
    file_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/hapref_chr21_filtered.vcf.bgz'
    bcf_file = pysam.VariantFile(file_path)


    output_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.vcf.gz'
    output_vcf = pysam.VariantFile(output_file, 'w', header=bcf_file.header)

    # Specify the sample you want to extract
    target_sample = "NA12878"

    # Create an empty list to store the records
    records = []

    # Iterate through the records in the BCF file
    for record in bcf_file.fetch():
        if target_sample in record.samples:
            # Extract genotype information
            genotype = record.samples[target_sample]["GT"]
            # Add record to DataFrame if it's relevant
            records.append({
                "CHROM": record.chrom,
                "POS": record.pos,
                "REF": record.ref,
                "ALT": ",".join(record.alts) if record.alts else None,
                "GENOTYPE": genotype
            })
        # Write the record to the new VCF file
        output_vcf.write(record)

    # Close the output VCF file
    output_vcf.close()
    
    # Convert the list of records into a DataFrame
    df = pd.DataFrame(records)

    # Filter for heterozygous genotypes (1|0 or 0|1)
    df["GENOTYPE"] = df["GENOTYPE"].apply(lambda x: "|".join(map(str, x)))  # Convert tuples to string
    # hetero_df = df[df["GENOTYPE"].isin(["1|0", "0|1"])].reset_index(drop=True)

    # # Sort by position and add simulated positions
    # hetero_df = hetero_df.sort_values(by="POS").reset_index(drop=True)
    # hetero_df["sim_pos"] = hetero_df.index
    # hetero_df.to_csv('/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv', index=False)
    df = df.sort_values(by="POS").reset_index(drop=True)
    # df['difference'] = df['POS'].diff()
    df.to_csv('/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv', index=False)


def extract_sample_from_vcf(file_path, target_sample):
    # target_sample = 'HG00514'
    file_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/hapref_chr21_filtered.vcf.bgz'
    bcf_file = pysam.VariantFile(file_path)

    output_file = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/maf0.01_hapref_chr21_filtered_{}.vcf.gz'.format(target_sample)
    output_vcf = pysam.VariantFile(output_file, 'w', header=bcf_file.header)

    # Specify the sample you want to extract
    # target_sample = "NA12878"

    # Create an empty list to store the records
    records = []

    # Iterate through the records in the BCF file
    for record in bcf_file.fetch():
        if target_sample in record.samples:
            # Extract genotype information
            genotype = record.samples[target_sample]["GT"]
            # Add record to DataFrame if it's relevant
            records.append({
                "CHROM": record.chrom,
                "POS": record.pos,
                "REF": record.ref,
                "ALT": ",".join(record.alts) if record.alts else None,
                "GENOTYPE": genotype
            })
        # Write the record to the new VCF file
        output_vcf.write(record)

    # Close the output VCF file
    output_vcf.close()
    
    # Convert the list of records into a DataFrame
    df = pd.DataFrame(records)

    # Filter for heterozygous genotypes (1|0 or 0|1)
    df["GENOTYPE"] = df["GENOTYPE"].apply(lambda x: "|".join(map(str, x)))  # Convert tuples to string
    # hetero_df = df[df["GENOTYPE"].isin(["1|0", "0|1"])].reset_index(drop=True)

    # # Sort by position and add simulated positions
    # hetero_df = hetero_df.sort_values(by="POS").reset_index(drop=True)
    # hetero_df["sim_pos"] = hetero_df.index
    # hetero_df.to_csv('/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv', index=False)
    df = df.sort_values(by="POS").reset_index(drop=True)
    # df['difference'] = df['POS'].diff()
    df_name = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/maf0.01_hapref_chr21_filtered_{}.csv'.format(target_sample)
    df.to_csv(df_name, index=False)


def check_compatibility_vcf_fasta():
    # vcf_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/_EGAZ00001239288_HRC.r1-1.EGA.GRCh37.chr21.haplotypes.vcf.gz'
    vcf_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/hapref_chr21_filtered.vcf.gz'

    vcf_df = pd.read_csv(vcf_path, skiprows=48, nrows=10, sep='\t')

    hetero_df = extract_column_NA12878(vcf_path)
    filtered_positions = hetero_df['POS'].values

    # Path to the FASTA file
    # chr21_fasta_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/chr21.fna'
    chr21_fasta_path = '/mnt/research/aguiarlab/data/hg19/chr21.fa'
    chr21_legend_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/_EGAZ00001239288_HRC.r1-1.EGA.GRCh37.chr21.legend.gz'
    # Read the FASTA file, skipping the first line (chromosome name line)
    with open(chr21_fasta_path, 'r') as f:
        fasta_file = f.readlines()

    # Combine the sequence lines into a single string, skipping the first line
    fasta_sequence = "".join(line.strip() for line in fasta_file if not line.startswith(">"))

    legend_df = pd.read_csv(chr21_legend_path, sep=' ', compression='gzip')

    for i in range(10):
        pos = legend_df.loc[i, 'position']
        ref_pos = legend_df.loc[i, 'a0']
        # alt_pos = legend_df.loc[i, 'a1']    
        print(pos, ref_pos, fasta_sequence[pos - 1], legend_df.loc[i, 'a1'])
        # print('----------------------------------')

    legend_df = legend_df[legend_df['position'].isin(list(filtered_positions))].reset_index(drop=True)
    for i in range(10):
        pos = legend_df.loc[i, 'position']
        ref_pos = legend_df.loc[i, 'a0']
        # alt_pos = legend_df.loc[i, 'a1']
        print(pos, ref_pos, fasta_sequence[pos - 1], legend_df.loc[i, 'a1'])


def check_integrity_vcf_fasta_sim_fasta_ref_():
    fasta_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/contig_10_ploidy_6.fa'
    vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/AWRI_ploidy6_contig10.vcf'
    genotypes_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/haplotypes.csv'

    vcf_in = pysam.VariantFile(vcf_path)
    # original_snps = {record.pos - 1: record for record in vcf_in.fetch()}  # Map positions to records
    
    def read_fasta_file(file_path):
        haplotype_dict = {}
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for i in range(0, len(lines), 2):
                key = lines[i].strip()[1:]  # Remove '>' and strip any whitespace
                value = lines[i + 1].strip()  # Strip any whitespace
                haplotype_dict[key] = value
        return haplotype_dict

    hap_dict = read_fasta_file(fasta_path)
    # [len(hap_dict[k]) for k in hap_dict.keys()]
    
    gen_df = pd.read_csv(genotypes_path)
    gen_df = gen_df.loc[gen_df.index.repeat(3)].reset_index(drop=True)
    snp_counter = 0
    for record in vcf_in.fetch():
        pos = record.pos - 1
        ref = record.ref
        alts = record.alts[:1] 
        print(snp_counter, pos, 'REF:', ref, 'ALT:', alts, 'HAPs:',[hap_dict[key][pos] for key in sorted(hap_dict.keys())], 'GT:', gen_df.iloc[snp_counter].values)
        snp_counter += 1

    
def simulate_NA12878():
    
    beagle_config_NA12878 = {
        # "input_vcf_path": '/mnt/research/aguiarlab/proj/HaplOrbit/reference/chr21_NA12878/updated_NA12878_extracted.vcf',
        "input_vcf_path": "/mnt/research/aguiarlab/proj/HaplOrbit/reference/chr21_NA12878/updated_NA12878_extracted_no_duplicate.vcf",
        "contig_fasta": '/mnt/research/aguiarlab/proj/HaplOrbit/reference/chr21_NA12878/GRCh37_chr21.fna',
        "main_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878_2',
        "art_path": 'art_illumina',
        "extract_hairs_path": 'extractHAIRS',
        "n_samples": 100, 
        "target_spacing": 100,
        "densify_snps": False, 
        "contig_lens": [100000], 
        "ploidies": [3, 4, 6, 8],
        "coverages": [10, 30, 50, 70, 100],
        "read_length": 150,
        "mean_insert_length": 500,
        "std_insert_length": 50, 
        "sbatch_str": ''
        }

    simulator = SimulatorNA12878(beagle_config_NA12878)
    simulator.simulate()


def extract_chromosome(fasta_path, chrom, out_path=None):
    """
    Extract one chromosome/contig from a plain-text FASTA by matching the
    first token after '>' (e.g., '>Chr8 something' -> 'Chr8').

    If out_path is given, writes the FASTA block there.
    Otherwise returns the FASTA-formatted string.
    """
    out_lines = []
    found = False
    write_block = False

    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                # End current block if we were writing and we hit a new header
                if found and not write_block:
                    break
                # Decide whether this header matches the requested chrom
                key = line[1:].split()[0] if len(line) > 1 else ""
                write_block = (key == chrom)
                if write_block:
                    # out_lines.append(line)
                    out_lines.append('>'+key+'\n')
                    found = True
            else:
                if write_block:
                    out_lines.append(line)

    if not found:
        raise ValueError(f"Contig '{chrom}' not found in {fasta_path}")

    if out_path:
        with open(out_path, "w", encoding="utf-8") as out:
            out.writelines(out_lines)
    else:
        return "".join(out_lines)


def tsv_chrom_to_vcf(tsv_path, chromosome, out_vcf_path, sample_name="SAMPLE"):
    """
    Convert a whitespace-delimited variants table to a VCF for one chromosome.
    Expected leading columns:
      var_id contig varpos ref_allele alt_allele hap_1 hap_2 [hap_3 hap_4 ...]
    Extra columns are ignored. Lines with too few columns are skipped.

    Writes a phased GT built from all hap_* columns present.
    """
    strip_prefix = re.compile(r'^\s*Block\s+length\s+\d+\s+')

    with open(out_vcf_path, "w", encoding="utf-8") as out, \
         open(tsv_path, "r", encoding="utf-8") as fh:

        # Minimal header
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##contig=<ID={chromosome}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype from hap_* columns (0=REF,1=ALT)">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        first_line = True
        for line in fh:
            if first_line:
                line = strip_prefix.sub("", line, count=1)  # remove "Block length <num> "
                first_line = False

            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()  # whitespace-split

            if parts[0].lower() == "var_id":
                continue  # skip header row
            if len(parts) < 7:           # <-- changed from 9 to 7 (supports diploid+)
                continue

            var_id, contig, pos, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            if contig != chromosome:
                continue

            haps = parts[5:]             # <-- changed: take all hap_* columns present
            gt_fields = []
            for x in haps:
                if x in {"0", "1"}:
                    gt_fields.append(x)
                elif x == ".":
                    gt_fields.append(".")
                else:
                    gt_fields.append(".")  # unknown -> missing
            gt = "|".join(gt_fields)

            out.write(f"{contig}\t{pos}\t{var_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")



def print_fasta_lengths(fasta_path, header_prefix=None):
    """
    Print 'header<TAB>length' for each FASTA record in fasta_path.
    If header_prefix is given (e.g., 'haplotype_1'), only print records whose
    header (first token after '>') starts with that prefix.
    """
    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                print(line, end="")           # print header line itself
            else:
                print(len(line.rstrip("\r\n")))  # length of sequence line (no newline)





def write_concat_chromosome_simple(fasta_paths, chrom, out_path):
    """
    Given a list of FASTA file paths, extract the block whose first header token == `chrom`
    from each file and write a single multi-FASTA to `out_path`, with headers named by
    hap number in the filename (e.g., '*hap1*' -> '>haplotype_1'). Sequences are one line.
    """
    with open(out_path, "w", encoding="utf-8") as out:
        for raw_fp in fasta_paths:
            fp = os.path.expanduser(raw_fp.strip())  # handle spaces, ~
            fname = os.path.basename(fp)

            # hap number from filename (e.g., hap1, hap_2, haplotype-3)
            m = re.search(r'hap(?:lotype[_-]?)?(\d+)', fname, re.IGNORECASE)
            if not m:
                raise ValueError(f"Cannot determine hap number from filename: {fname}")
            hap_no = int(m.group(1))

            # extract sequence for `chrom`
            seq_chunks, writing, found = [], False, False
            with open(fp, "r", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith(">"):
                        key = line[1:].split()[0] if len(line) > 1 else ""
                        writing = (key == chrom)
                        if writing:
                            found = True
                        continue
                    if writing:
                        seq_chunks.append("".join(line.split()))  # strip whitespace

            if not found:
                raise ValueError(f"Contig '{chrom}' not found in file: {fp}")

            seq = "".join(seq_chunks).upper()
            out.write(f">haplotype_{hap_no}\n{seq}\n")



def simulate_long_reads():
    pbsim_location = '/home/mah19006/pHapCompass_v2/pbsim3/src/pbsim'
    model_path = '/home/mah19006/pHapCompass_v2/pbsim3/data/QSHMM-ONT-HQ.model' # QSHMM-ONT-HQ.model
    # model_path = '/home/mah19006/pHapCompass_v2/pbsim3/data/QSHMM-RSII.model'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_auto'
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    coverages = [2, 10, 30, 50, 70, 100]
    ploidies = [2, 3, 4, 6, 8]
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    chrom = 'Chr1'
    for ploidy in ploidies:
        for mr in mut_rates:
            for sample in samples:
                
                haplotype_path = os.path.join(hap_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2))
                haplotype_names = [os.path.join(haplotype_path, str(sample).zfill(2) + '_' + 'hap' + str(i + 1) + '.fa.gz') for i in range(ploidy)]
                # haplotype_names = ['test_hap1.fa', 'test_hap2.fa', 'test_hap3.fa', 'test_hap4.fa']
                # out_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/potato_sim_fasta'
                variation_txt_path = os.path.join(haplotype_path, str(sample).zfill(2) + '_varianthaplos.txt')
                out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                fasta_file_name = os.path.join(haplotype_path, chrom + '_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(2) +'.fa.gz')
                write_concat_chromosome_simple(haplotype_names, chrom, fasta_file_name)
                tsv_chrom_to_vcf(variation_txt_path, chrom, out_vcf_path, sample_name="SAMPLE")

                for coverage in coverages:
                    cov_path = os.path.join(dataset_path, ploidy, coverage, str(sample).zfill(2))
                    if not os.path.exists(cov_path):
                        os.makedirs(cov_path, exist_ok=True)
                    # pbsim_command = 'pbsim --strategy wgs --method qshmm --qshmm data/QSHMM-ONT.model --depth 30 --genome ref.fa       --prefix ont_r10x_30x'
                    # pbsim_command = '{} --strategy wgs --method qshmm --qshmm {} --depth {} --genome {} --pass-num 10 --prefix sim'.format(pbsim_location, model_path, coverage, fasta_file_name)
                    pbsim_command = '{} --strategy wgs --method qshmm --qshmm {} --depth {} --genome {}  --accuracy-mean 0.999 0.0005 --prefix sim'.format(pbsim_location, model_path, coverage, fasta_file_name)
                    # pbsim_command = 

# def simulate_na12878():

#     beagle_config_human = {
#         "snp_df_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv',
#         # "input_vcf_path": '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/hapref_chr21_filtered.vcf.bgz',
#         "input_vcf_path": '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/updated_NA12878_extracted.vcf',
#         "contig_fasta": '/mnt/research/aguiarlab/data/hg19/chr21.fa',
#         "main_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878',
#         "art_path": 'art_illumina',
#         "extract_hairs_path": 'extractHAIRS',
#         "n_samples": 2, 
#         "target_spacing": 100,
#         "densify_snps": False, 
#         "contig_lens": [100], 
#         "ploidies": [3],
#         "coverages": [10],
#         "read_length": 150,
#         "mean_insert_length": 2000,
#         "std_insert_length": 300
#         }


#     simulator = Simulator(beagle_config_human)
#     # simulator.generate_genomes_fasta()

#     # simulator.simulate()
    
#     inputs = make_inputs_for_generate_qoutient_graph(simulator)

#     pool = Pool(2)
#     pool.map(generate_quotient_graph, inputs)

#     # for inp in inputs:
#     #     print(inp[4])
#     #     generate_quotient_graph(inp)


def simulate_awri():

    server3_config_AWRI = {
        "snp_df_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv',
        "input_vcf_path": '/mnt/research/aguiarlab/proj/HaplOrbit/SRR942191/vcf_files/SRR942191_AHIQ01000001.1.vcf',
        # "contig_fasta": '/mnt/research/aguiarlab/proj/HaplOrbit/reference/AWRI1499/contigs_noamb/contig_AHIQ01000001.1.fa',
        "contig_fasta": '/mnt/research/aguiarlab/proj/HaplOrbit/reference/AWRI1499/contigs/contig_AHIQ01000001.1.fa',
        "main_path": '/home/mok23003/BML/HaplOrbit/simulated_data_awri',
        "art_path": 'art_illumina',
        "extract_hairs_path": 'extractHAIRS',
        "n_samples": 100, 
        "target_spacing": 100,
        "densify_snps": False, 
        "contig_lens": [10, 100, 1000], 
        "ploidies": [3, 4, 6, 8, 10],
        "coverages": [10, 20, 30, 40, 50],
        "read_length": 75,
        "mean_insert_length": 300,
        "std_insert_length": 30, 
        "sbatch_str": ''
        }

    beagle_config_AWRI = {
        "snp_df_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv',
        "input_vcf_path": '/mnt/research/aguiarlab/proj/HaplOrbit/SRR942191/vcf_files/SRR942191_AHIQ01000001.1.vcf',
        # "contig_fasta": '/mnt/research/aguiarlab/proj/HaplOrbit/reference/AWRI1499/contigs_noamb/contig_AHIQ01000001.1.fa',
        "contig_fasta": '/mnt/research/aguiarlab/proj/HaplOrbit/reference/AWRI1499/contigs/contig_AHIQ01000001.1.fa',
        "main_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri',
        "art_path": 'art_illumina',
        "extract_hairs_path": 'extractHAIRS',
        "n_samples": 100, 
        "target_spacing": 100,
        "densify_snps": False, 
        "contig_lens": [10, 100, 1000], 
        "ploidies": [3, 4, 6, 8, 10],
        "coverages": [10, 20, 30, 40, 50],
        "read_length": 75,
        "mean_insert_length": 300,
        "std_insert_length": 30, 
        "sbatch_str": ''
        }

    xanadu_config_AWRI = {
        "snp_df_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv',
        "input_vcf_path": '/labs/Aguiar/pHapCompass/datasets/SRR942191/vcf_files/SRR942191_AHIQ01000001.1.vcf',
        # "contig_fasta": '/mnt/research/aguiarlab/proj/HaplOrbit/reference/AWRI1499/contigs_noamb/contig_AHIQ01000001.1.fa',
        "contig_fasta": '/labs/Aguiar/pHapCompass/references/AWRI1499/contigs/contig_AHIQ01000001.1.fa',
        "main_path": '/labs/Aguiar/pHapCompass/datasets/simulated/simulated_data_awri',
        "art_path": 'art_illumina',
        "extract_hairs_path": 'extractHAIRS',
        "n_samples": 1, 
        "target_spacing": 100,
        "densify_snps": False, 
        "contig_lens": [10, 100, 1000], 
        "ploidies": [3, 4, 6, 8, 10],
        "coverages": [10, 20, 30, 40, 50],
        "read_length": 75,
        "mean_insert_length": 300,
        "std_insert_length": 30,
        "sbatch_str": "#!/bin/bash\n#SBATCH --job-name={}\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=50G\n#SBATCH --mail-user=marjan.hosseini@uconn.edu\n#SBATCH -o {}.out\n#SBATCH -e {}.err".format('job', 'sim', 'sim')
        }

    # simulator = SimulatorAWRI(xanadu_config_AWRI)
    simulator = SimulatorAWRI(beagle_config_AWRI)
    # simulator = SimulatorAWRI(server3_config_AWRI)
    # simulator.generate_genomes_fasta()
    # # simulator.simulate()
    
    # next_inputs = make_inputs_for_running_FFBS(simulator)
    # print('number of inputs:', len(next_inputs))
    # pool = Pool(30)
    # pool.map(run_FFBS_quotient, next_inputs)

    # next_inputs = make_inputs_for_run_count(simulator)
    # print('number of inputs:', len(next_inputs))
    # # pool = Pool(30)
    # # pool.map(run_FFBS_quotient_count, next_inputs)
    # for inp in next_inputs:
    #     print(inp[4])
    #     run_FFBS_quotient_count(inp)

    # next_inputs = make_inputs_for_run_likelihood(simulator)
    # print('number of inputs:', len(next_inputs))
    # # pool = Pool(10)
    # # pool.map(run_FFBS_quotient_likelihood, next_inputs)
    # for inp in next_inputs:
    #     print(inp[4])
    #     run_FFBS_quotient_likelihood(inp)

