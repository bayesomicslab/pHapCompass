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
from FFBS.FFBS_quotient_graph import *
import time


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
                        print(pos, ref, contig_seq[pos], alts)
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

                haplotype_df = haplotype_df.T
                # haplotype_df.columns = [f"haplotype_{i + 1}" for i in range(ploidy)]
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
                        command = 'bwa mem {} {}/{}1.fq {}/{}2.fq > {}/{}.sam\n'.format(ref_fasta, fastq_path, str(rd).zfill(2), fastq_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
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

                # haplotypes = [[] for _ in range(ploidy)]
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
                    if not 'AHIQ01000' in str(line):
                        # print(line)
                        new_header.add_line(str(line))
                # for k in range(ploidy):
                #     new_header.contigs.add(f"haplotype_{k + 1}", length=len(contig_seq))
                # new_header.contigs.add(f"contig_{contig_len}", length=len(contig_seq))
                # new_header.contigs.add('AHIQ01000001.1', length=len(contig_seq))
                new_header.add_line('##contig=<ID=AHIQ01000001.1_{},length={}>'.format(contig_len, len(contig_seq[:last_snp])))
                new_header.add_sample("AWRI")

                output_vcf = os.path.join(ploidy_path, f'AWRI_ploidy{ploidy}_contig{contig_len}.vcf.gz')
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
                        print(pos, ref, contig_seq[pos], alts)
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
                            contig='AHIQ01000001.1_{}'.format(contig_len),
                            # pos=record.pos,
                            start=record.start,
                            stop=record.stop,
                            alleles=(ref, alts[0]),
                            qual=75,
                            filter=record.filter.keys(),
                            # info=record.info,
                            info=processed_info)
                        new_record.samples["AWRI"]["GT"] = genotype
                        new_record.samples["AWRI"].phased = False
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
                    f_out.write('>AHIQ01000001.1_{}\n'.format(contig_len))
                    f_out.write("".join(new_fasta) + "\n")

                haplotype_df = haplotype_df.T
                # haplotype_df.columns = [f"haplotype_{i + 1}" for i in range(ploidy)]
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
                        command = 'bwa mem {} {}/{}1.fq {}/{}2.fq > {}/{}.sam\n'.format(ref_fasta, fastq_path, str(rd).zfill(2), fastq_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
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
                
                # this_vcf_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'AWRI_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
                # this_vcf_sorted_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'sorted_AWRI_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
                # unzipped_vcf = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'sorted_AWRI_ploidy{}_contig{}.vcf'.format(ploidy, contig_len))
                # sort_cmd = 'bcftools sort {} -Oz -o {}\n'.format(this_vcf_path, this_vcf_sorted_path)
                # idx_cmd = 'tabix -p vcf {}\n'.format(this_vcf_sorted_path)
                # unzip_cmd = 'gunzip {}\n\n'.format(this_vcf_sorted_path)
                # to_print += sort_cmd
                # to_print += idx_cmd
                # to_print += unzip_cmd
                vcf_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'AWRI_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
                to_print = 'gunzip {}\n\n'.format(vcf_path)
                
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
        print('save main.sh')
        with open(os.path.join(self.main_path, 'main.sh'), 'w') as f:
            f.write(self.main_sh)


def check_integrity_vcf_fastas_genotypes():
    fasta_ref_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test/contig_10/ploidy_6/contig_10.fa'
    fasta_sim_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test/contig_10/ploidy_6/contig_10_ploidy_6.fa'
    vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test/contig_10/ploidy_6/AWRI_ploidy6_contig10.vcf.gz'
    genotypes_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test/contig_10/ploidy_6/haplotypes.csv'

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

    hap_dict = read_fasta_file(fasta_sim_path)
    ref_fasta = read_fasta_file(fasta_ref_path)
    # [len(hap_dict[k]) for k in hap_dict.keys()]
    
    gen_df = pd.read_csv(genotypes_path)

    # gen_df = gen_df.loc[gen_df.index.repeat(3)].reset_index(drop=True)
    snp_counter = 0
    for record in vcf_in.fetch():
        pos = record.pos - 1
        ref = record.ref
        alts = record.alts[:1]
        print(snp_counter, pos, 'REF:', ref, 'ALT:', alts, 'HAPs:',[hap_dict[key][pos] for key in sorted(hap_dict.keys())], 'GT:', gen_df.loc[:, str(snp_counter)], 'REF_fasta:', ref_fasta['AHIQ01000001.1_10'][pos])
        print('----------------------------')
        snp_counter += 1


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
        "main_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test',
        "art_path": 'art_illumina',
        "extract_hairs_path": 'extractHAIRS',
        "n_samples": 10, 
        "target_spacing": 100,
        "densify_snps": False, 
        "contig_lens": [10, 100, 1000], 
        "ploidies": [3, 4, 6],
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


    simulator = SimulatorAWRI(beagle_config_AWRI)

    # simulator.simulate()
    bam_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test/contig_10/ploidy_6/cov_10/bam/00.bam'
    fragments = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue  # Skip unmapped or non-primary alignments

            read_name = read.query_name

            chrom = bam.get_reference_name(read.reference_id)  # Chromosome for this read
            
            aligned_positions = read.get_reference_positions(full_length=True)
            read_sequence = read.query_sequence
            if read_name == 'haplotype_1-18':
                print(read_name)
                # stop
            # Process only positions in the current chromosome's variants
            fragment_positions = []
            fragment_bases = []
            for pos in variants[chrom]:
                if pos - 1 in aligned_positions:  # BAM positions are 0-based
                    base_index = aligned_positions.index(pos - 1)
                    base = read_sequence[base_index]
                    if base == variants[chrom][pos]["ref"]:
                        fragment_bases.append("0")
                    elif base in variants[chrom][pos]["alt"]:
                        fragment_bases.append("1")
                    else:
                        fragment_bases.append("-")
                    fragment_positions.append(pos)  # Keep as 1-based for output


def simulate_NA12878():
    
    beagle_config_NA12878 = {
        # "snp_df_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv',
        "input_vcf_path": '/mnt/research/aguiarlab/proj/HaplOrbit/reference/chr21_NA12878/updated_NA12878_extracted.vcf',
        "contig_fasta": '/mnt/research/aguiarlab/proj/HaplOrbit/reference/chr21_NA12878/GRCh37_chr21.fna',
        "main_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test',
        "art_path": 'art_illumina',
        "extract_hairs_path": 'extractHAIRS',
        "n_samples": 10, 
        "target_spacing": 100,
        "densify_snps": False, 
        "contig_lens": [100], 
        "ploidies": [6],
        "coverages": [10, 50, 100],
        "read_length": 150,
        "mean_insert_length": 500,
        "std_insert_length": 50, 
        "sbatch_str": ''
        }

    simulator = SimulatorNA12878(beagle_config_NA12878)
    simulator.simulate()


if __name__ == '__main__':
    simulate_NA12878()