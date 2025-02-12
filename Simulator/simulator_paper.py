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


# def get_block_info(quotient_g, predicted_haplotypes, true_haplotypes, fragment_model):
#     component_labels, _ = gt.label_components(quotient_g.graph)
#     components = {}
#     for v in quotient_g.graph.vertices():
#         comp_id = component_labels[v]  # Get component ID of the vertex
#         if comp_id not in components:
#             components[comp_id] = {'blocks': []}
#         components[comp_id]['blocks'].append(quotient_g.graph.vertex_properties["v_label"][v])

#     for key in components.keys():
#         block = components[key]['blocks']
#         positions = sorted(set(int(num) for r in block for num in r.split('-')))
#         positions_ind = [p-1 for p in positions]
#         block_pred_haplotype = predicted_haplotypes[positions_ind].to_numpy()
#         block_true_haplotype = true_haplotypes[positions_ind].to_numpy()
#         block_vector_error_rate, block_vector_error, _, _ = compute_vector_error_rate(block_pred_haplotype, block_true_haplotype)
#         block_accuracy, _ = calculate_accuracy(block_pred_haplotype, block_true_haplotype)
#         block_mismatch_error, _ = calculate_mismatch_error(block_pred_haplotype, block_true_haplotype)
#         block_mec_ = mec(block_pred_haplotype, fragment_model.fragment_list)
#         components[key]['evaluation'] = {'vector_error_rate': block_vector_error_rate, 'vector_error': block_vector_error, 
#                                          'accuracy': block_accuracy, 'mismatch_error': block_mismatch_error, 'mec': block_mec_}
#         components[key]['block_size'] = len(positions_ind)

#     block_info = {'vector_error_rate': np.mean([components[key]['evaluation']['vector_error_rate'] for key in components.keys()]), 
#                   'vector_error': np.mean([components[key]['evaluation']['vector_error'] for key in components.keys()]),
#                   'accuracy': np.mean([components[key]['evaluation']['accuracy'] for key in components.keys()]),
#                   'mismatch_error': np.mean([components[key]['evaluation']['mismatch_error'] for key in components.keys()]),
#                   'mec': np.mean([components[key]['evaluation']['mec'] for key in components.keys()]), 
#                   'average_block_size': np.mean([components[key]['block_size'] for key in components.keys()])}

#     return block_info, components


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


# def emissions_v2(ploidy, quotient_g, quotient_g_v_label_reversed, error_rate):
#     """
#     The only difference in version 2 is that we are using quotient graph directly:
#     quotient_g.graph ==> quotient_g
#     """
#     emission_dict = {}
#     # Calculate emissions for each state and populate the emission probability matrix
#     for state in quotient_g_v_label_reversed.keys():
#         emission_dict[state] = {}

#         node_elements = state.split('-')  # For example, '1-2' -> ['1', '2']
#         node_length = len(node_elements)  # Determine the number of elements in the node

#         # Generate all possible combinations of 0 and 1 based on the node length
#         possible_emissions = generate_binary_combinations(node_length)
#         v = quotient_g_v_label_reversed[state]
#         phasings = quotient_g.vertex_properties["v_weights"][v]['weight'].keys()
#         for phasing in phasings:
#             emission_dict[state][phasing] = {}
#             phasing_np = str_2_phas_1(phasing, ploidy)  # Compute phasing for the state key
#             for emission in possible_emissions:
#                 likelihood = compute_likelihood(np.array(emission), phasing_np, error_rate)
#                 emission_dict[state][phasing][''.join([str(e) for e in emission])] = likelihood

#     return emission_dict


# def transition_matrices_v2(quotient_g, edges_map_quotient, ploidy, config, fragment_model):
#     """
#     The only difference in version 2 is that we are using quotient graph directly:
#     quotient_g.graph ==> quotient_g
#     """
#     transitions_dict = {}
#     transitions_dict_extra = {}
#     for edge in edges_map_quotient.keys():
#         transitions_dict_extra[edge] = {}
#         source = edges_map_quotient[edge][0]
#         target = edges_map_quotient[edge][1]
#         source_weights = quotient_g.vertex_properties["v_weights"][source]['weight']
#         target_weights = quotient_g.vertex_properties["v_weights"][target]['weight']
#         source_label = quotient_g.vertex_properties["v_label"][source]
#         target_label = quotient_g.vertex_properties["v_label"][target]
#         common_ff, common_sf = find_common_element_and_index(source_label, target_label)
#         source_phasings = list(source_weights.keys())
#         target_phasings = list(target_weights.keys())
#         # transitions_dict = {'source': source_phasings, 'target': target_phasings}
#         transitions_mtx = np.zeros((len(source_phasings), len(target_phasings)))

#         for i, ffstr in enumerate(source_phasings):
#             for j, sfstr in enumerate(target_phasings):
#                 transitions_dict_extra[edge][str(i) + '-' + str(j)] = {}
#                 transitions_dict_extra[edge][str(i) + '-' + str(j)]['source_phasing'] = ffstr
#                 transitions_dict_extra[edge][str(i) + '-' + str(j)]['target_phasing'] = sfstr
#                 transitions_dict_extra[edge][str(i) + '-' + str(j)]['matched_phasings'] = {}
#                 matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, ploidy), str_2_phas_1(sfstr, ploidy), common_ff, common_sf, source_label, target_label)
#                 sorted_phasings = []
#                 for mtx in matched_phasings:
#                     sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
#                     sorted_phasings.append(sorted_matrix)
                
#                 matched_phasings_str = list(set([phas_2_str(pm) for pm in sorted_phasings]))
#                 # print(i, ffstr, j, sfstr)
#                 # print('matched phasings:', matched_phasings_str, len(matched_phasings_str))
#                 # if len(matched_phasings_str) > 1:
#                 #     print('More than one matching phasing')
#                 #     # stop
#                 poss = sorted(list(set([int(ss) for ss in source_label.split('-')] + [int(tt) for tt in target_label.split('-')])))
#                 match_reads = get_matching_reads_for_positions([int(i) for i in poss], fragment_model.fragment_list)
#                 wei = 0
#                 for phas in matched_phasings_str:
#                     this_phas_weight = 0
#                     for indc, this_po, obs in match_reads:
#                         this_phas_read_weight = compute_likelihood_generalized_plus(np.array(obs), str_2_phas_1(phas, ploidy), indc, list(range(len(indc))), 
#                                                                    config.error_rate)
#                         wei += this_phas_read_weight
#                         this_phas_weight += this_phas_read_weight
#                     transitions_dict_extra[edge][str(i) + '-' + str(j)]['matched_phasings'][phas] = this_phas_weight
#                 transitions_mtx[i, j] = wei

#         transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)
#         transitions_dict[edge] = transitions_mtx
#     return transitions_dict, transitions_dict_extra


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
    plt.hist(differences, bins=50, color='blue', alpha=0.7, edgecolor='black')
    plt.xlabel('Distances (bp)')
    plt.ylabel('Frequency')
    plt.title('Histogram of Differences Between Consecutive VCF Positions')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xscale('log')
    # plt.xlim(0, 100000)
    # plt.show()
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

    hetero_df = extract_column_NA12878()
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


# def make_inputs_for_generate_qoutient_graph(simulator):
#     inputs = []
#     simulator.contig_lens = [100]
#     simulator.ploidies = [3, 4, 6, 8]
#     for contig_len in simulator.contig_lens:
#         for ploidy in simulator.ploidies:
#             # stop
#             genotype_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
#             # genotype_df = pd.read_csv(genotype_path)
#             for coverage in simulator.coverages:
#                 this_cov_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
#                 frag_path = os.path.join(this_cov_path, 'frag')
#                 frag_graph_path = os.path.join(this_cov_path, 'fgraph')
#                 quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
#                 qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')

#                 if not os.path.exists(frag_graph_path):
#                     os.makedirs(frag_graph_path)
#                 if not os.path.exists(quotient_graph_path):
#                     os.makedirs(quotient_graph_path)
#                 if not os.path.exists(qgraph_reverse_maps_path):
#                     os.makedirs(qgraph_reverse_maps_path)

#                 existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
#                 existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
#                 existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
#                 existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
#                 existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
#                 existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]

#                 for rd in range(simulator.n_samples):
#                     if 'qg_e_label_' + str(rd).zfill(2) + '.pkl' not in existing_files_qg_e or 'qg_v_label_' + str(rd).zfill(2) + '.pkl' not in existing_files_qg_v or 'fg_e_label_' + str(rd).zfill(2) + '.pkl' not in existing_files_fg_e or 'fg_v_label_' + str(rd).zfill(2) + '.pkl' not in existing_files_fg_v or str(rd).zfill(2) + '.gt.gz' not in existing_fg or str(rd).zfill(2) + '.gt.gz' not in existing_qg:
#                         inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path]
#                         inputs.append(inp)
#     return inputs


# def make_inputs_for_run_count(simulator):
#     simulator.contig_lens = [100]
#     inputs = []
#     for contig_len in simulator.contig_lens:
#         for ploidy in simulator.ploidies:
#             # stop
#             genotype_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
#             # genotype_df = pd.read_csv(genotype_path)
#             for coverage in simulator.coverages:
#                 this_cov_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
#                 frag_path = os.path.join(this_cov_path, 'frag')
#                 frag_graph_path = os.path.join(this_cov_path, 'fgraph')
#                 quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
#                 qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')
#                 results_path = os.path.join(this_cov_path, 'results_counts')

#                 if not os.path.exists(frag_graph_path):
#                     os.makedirs(frag_graph_path)
#                 if not os.path.exists(quotient_graph_path):
#                     os.makedirs(quotient_graph_path)
#                 if not os.path.exists(qgraph_reverse_maps_path):
#                     os.makedirs(qgraph_reverse_maps_path)
#                 if not os.path.exists(results_path):
#                     os.makedirs(results_path)

#                 # existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
#                 # existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
#                 # existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
#                 # existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
#                 # existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
#                 # existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]
#                 existing_results = [ff for ff in os.listdir(results_path) if 'FFBS' in ff]
#                 # existing_results = []
#                 for rd in range(simulator.n_samples):
#                     if 'FFBS_{}.pkl'.format(str(rd).zfill(2)) not in existing_results and \
#                         '{}.gt.gz'.format(str(rd).zfill(2)) in os.listdir(quotient_graph_path) and \
#                         'qg_e_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path) and \
#                         'qg_v_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path):
#                         inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path, results_path]
#                         # inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path]

#                         inputs.append(inp)
#         inputs = sorted(inputs, key=lambda x: x[0])
#     return inputs


# def make_inputs_for_run_likelihood(simulator):
#     simulator.contig_lens = [100]
#     simulator.ploidies = [3, 4, 6]
#     inputs = []
#     for contig_len in simulator.contig_lens:
#         for ploidy in simulator.ploidies:
#             # stop
#             genotype_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
#             # genotype_df = pd.read_csv(genotype_path)
#             for coverage in simulator.coverages:
#                 this_cov_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
#                 frag_path = os.path.join(this_cov_path, 'frag')
#                 frag_graph_path = os.path.join(this_cov_path, 'fgraph')
#                 quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
#                 qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')
#                 results_path = os.path.join(this_cov_path, 'results_likelihood')

#                 if not os.path.exists(frag_graph_path):
#                     os.makedirs(frag_graph_path)
#                 if not os.path.exists(quotient_graph_path):
#                     os.makedirs(quotient_graph_path)
#                 if not os.path.exists(qgraph_reverse_maps_path):
#                     os.makedirs(qgraph_reverse_maps_path)
#                 if not os.path.exists(results_path):
#                     os.makedirs(results_path)

#                 # existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
#                 # existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
#                 # existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
#                 # existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
#                 # existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
#                 # existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]
#                 # existing_results = [ff for ff in os.listdir(results_path) if 'FFBS' in ff]
#                 existing_results = []
#                 # for rd in range(90, 100):
#                 for rd in range(simulator.n_samples):
#                     if 'FFBS_{}.pkl'.format(str(rd).zfill(2)) not in existing_results and \
#                         '{}.gt.gz'.format(str(rd).zfill(2)) in os.listdir(quotient_graph_path) and \
#                         'qg_e_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path) and \
#                         'qg_v_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path):
#                         print(frag_path)
#                         inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path, results_path]
#                         # inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path]

#                         inputs.append(inp)
#         inputs = sorted(inputs, key=lambda x: x[0], reverse=True)
#     return inputs


# def make_inputs_for_running_FFBS(simulator):
#     simulator.contig_lens = [10]
#     inputs = []
#     for contig_len in simulator.contig_lens:
#         for ploidy in simulator.ploidies:
#             # stop
#             genotype_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
#             # genotype_df = pd.read_csv(genotype_path)
#             for coverage in simulator.coverages:
#                 this_cov_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
#                 frag_path = os.path.join(this_cov_path, 'frag')
#                 frag_graph_path = os.path.join(this_cov_path, 'fgraph')
#                 quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
#                 qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')
#                 results_path = os.path.join(this_cov_path, 'results_algorithm_v2')

#                 if not os.path.exists(frag_graph_path):
#                     os.makedirs(frag_graph_path)
#                 if not os.path.exists(quotient_graph_path):
#                     os.makedirs(quotient_graph_path)
#                 if not os.path.exists(qgraph_reverse_maps_path):
#                     os.makedirs(qgraph_reverse_maps_path)
#                 if not os.path.exists(results_path):
#                     os.makedirs(results_path)

#                 # existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
#                 # existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
#                 # existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
#                 # existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
#                 # existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
#                 # existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]
#                 # existing_results = [ff for ff in os.listdir(results_path) if 'FFBS' in ff]
#                 existing_results = []
#                 for rd in range(simulator.n_samples):
#                     if 'FFBS_{}.pkl'.format(str(rd).zfill(2)) not in existing_results and \
#                         '{}.gt.gz'.format(str(rd).zfill(2)) in os.listdir(quotient_graph_path) and \
#                         'qg_e_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path) and \
#                         'qg_v_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path):
#                         inp = [frag_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path, results_path]
#                         inputs.append(inp)
#         inputs = sorted(inputs, key=lambda x: x[0])
#     return inputs


# def make_inputs_for_chordal_contraction(simulator):
#     simulator.contig_lens = [100]
#     k = 10
#     inputs = []
#     for contig_len in simulator.contig_lens:
#         for ploidy in simulator.ploidies:
#             # stop
#             genotype_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
#             # genotype_df = pd.read_csv(genotype_path)
#             for coverage in simulator.coverages:
#                 this_cov_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
#                 frag_path = os.path.join(this_cov_path, 'frag')
#                 frag_graph_path = os.path.join(this_cov_path, 'fgraph')
#                 quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
#                 qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')
#                 chordal_graph_path = os.path.join(this_cov_path, 'chgraph')

#                 if not os.path.exists(frag_graph_path):
#                     os.makedirs(frag_graph_path)
#                 if not os.path.exists(quotient_graph_path):
#                     os.makedirs(quotient_graph_path)
#                 if not os.path.exists(qgraph_reverse_maps_path):
#                     os.makedirs(qgraph_reverse_maps_path)
#                 if not os.path.exists(chordal_graph_path):
#                     os.makedirs(chordal_graph_path)

#                 # existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
#                 # existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
#                 # existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
#                 # existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
#                 # existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
#                 # existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]
#                 existing_results = [ff for ff in os.listdir(chordal_graph_path) if '.gt.gz' in ff]

#                 for rd in range(simulator.n_samples):
#                     if '{}.gt.gz'.format(str(rd).zfill(2)) not in existing_results and \
#                         '{}.gt.gz'.format(str(rd).zfill(2)) in os.listdir(quotient_graph_path):
#                         inp = [frag_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), chordal_graph_path, genotype_path, ploidy, k]
                        
#                         inputs.append(inp)
#         inputs = sorted(inputs, key=lambda x: x[0])
#     return inputs


# def chordal_contraction_graph_tool_top_k(inp):
#     this_frag_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, chordal_graph_path, genotype_path, ploidy, k = inp

#     # save_path, subg_id, subg, config, fragment_model, k = inp
#     # this_path = os.path.join(save_path, 'chordal_sub_' + str(subg_id) + '.gt.gz')

#     print('Working on:', os.path.join(this_frag_path, frag_file))
#     chordal_v_label_revered_path = os.path.join(this_reverse_maps_path, 'ch_v_label_' + frag_file.split('.')[0] + '.pkl')
#     chordal_e_label_revered_path = os.path.join(this_reverse_maps_path, 'ch_e_label_' + frag_file.split('.')[0] + '.pkl')


#     class Args:
#         def __init__(self):
#             self.vcf_path = 'example/62_ID0.vcf'
#             self.data_path = os.path.join(this_frag_path, frag_file)
#             # self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_90.frag.txt'
#             self.bam_path = 'example/example.bam'
#             self.genotype_path = genotype_path
#             self.ploidy = ploidy
#             self.error_rate = 0.001
#             self.epsilon = 0.0001
#             self.output_path = 'output'
#             self.root_dir = 'D:/UCONN/HaplOrbit'
#             self.alleles = [0, 1]

#     # Create the mock args object
#     args = Args()

#     # Initialize classes with parsed arguments
#     input_handler = InputHandler(args)

#     config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

#     fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
#     fragment_model.construct(input_handler, config)

#     quotient_g_path = os.path.join(this_quotient_coverage_path, frag_file.split('.')[0] + '.gt.gz')
#     quotient_g = gt.load_graph(quotient_g_path)

#     new_graph = quotient_g.copy()
#     e_weights = new_graph.edge_properties["e_weights"]
#     # new_graph.clear_filters()
#     e_entropy = new_graph.new_edge_property("double")

#     # Loop over edges and assign entropy from the e_weights property
#     for e in new_graph.edges():
#         e_entropy[e] = e_weights[e]['entropy']

#     new_graph.ep['e_entropy'] = e_entropy

#     chordless_cycles = get_chordless_cycles(new_graph)

#     to_be_removed_nodes = []

#     for cyc_id, cyc in enumerate(chordless_cycles):
#         # print(cyc_id)        
#         edges = [new_graph.edge(cyc[-1], cyc[0])]
#         for i in range(len(cyc) - 1):
#             edges += [new_graph.edge(cyc[i], cyc[i+1])]
#         edges = [x for x in edges if x is not None]
#         while len(edges) > 3:
#             min_edge = min(edges, key=lambda e: new_graph.ep['e_entropy'][e])
#             source_label = new_graph.vp['v_label'][min_edge.source()]
#             target_label = new_graph.vp['v_label'][min_edge.target()]
#             # new node positions
#             poss = sorted(set([int(nn) for nn in source_label.split('-')] + [int(nn) for nn in target_label.split('-')]))
#             # new vertex properties:
#             new_vertex_name = '-'.join([str(nnn) for nnn in poss])
#             vertex_weights = new_graph.ep['e_weights'][min_edge]
#             vertex_weights_appr = get_top_k_weights(vertex_weights, k)

#             new_graph.vertex_properties["v_weights"][min_edge.source()] = vertex_weights_appr
#             new_graph.vertex_properties["v_label"][min_edge.source()] = new_vertex_name

#             source_nbrs = [n for n in min_edge.source().all_neighbors() if n != min_edge.target()]
#             target_nbrs = [n for n in min_edge.target().all_neighbors() if n != min_edge.source()]
#             common_nbrs = set(source_nbrs).intersection(set(target_nbrs))

#             for n in common_nbrs:
                
#                 v_label = new_graph.vertex_properties["v_label"][n]
#                 # e_poss = sorted(set([int(nn) for nn in v_label.split('-')] + poss))
#                 # print(len(e_poss))
#                 # new_edge_name = '-'.join([str(nnn) for nnn in e_poss])
#                 sorted_labels = sort_nodes([new_vertex_name, v_label])
#                 new_edge_name = '--'.join(sorted_labels)
#                 (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

#                 first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
#                 second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
#                 final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)
#                 final_weights_appr = get_top_k_weights(final_weight, k)

#                 e1 = new_graph.edge(min_edge.source(), n)
#                 e2 = new_graph.edge(min_edge.target(), n)
                
#                 new_graph.edge_properties["e_weights"][e1] = final_weights_appr
#                 new_graph.edge_properties["e_label"][e1] = new_edge_name
#                 new_graph.edge_properties['e_entropy'][e1] = final_weights_appr['entropy']
#                 new_graph.remove_edge(e2)

#             for n in set(source_nbrs)-common_nbrs:
                
#                 v_label = new_graph.vertex_properties["v_label"][n]

#                 sorted_labels = sort_nodes([new_vertex_name, v_label])
#                 new_edge_name = '--'.join(sorted_labels)
#                 (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

#                 first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
#                 second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
#                 final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)
#                 final_weights_appr = get_top_k_weights(final_weight, k)


#                 e1 = new_graph.edge(min_edge.source(), n)
#                 # e2 = new_graph.edge(min_edge.target(), n)
#                 new_graph.edge_properties["e_weights"][e1] = final_weights_appr
#                 new_graph.edge_properties["e_label"][e1] = new_edge_name
#                 new_graph.edge_properties['e_entropy'][e1] = final_weights_appr['entropy']
#                 # new_graph.edge_properties["e_weights"][e2]

            
#             for n in set(target_nbrs)-common_nbrs:
                
#                 v_label = new_graph.vertex_properties["v_label"][n]
#                 sorted_labels = sort_nodes([new_vertex_name, v_label])
#                 new_edge_name = '--'.join(sorted_labels)

#                 (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

#                 first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
#                 second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
#                 final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)
#                 final_weights_appr = get_top_k_weights(final_weight, k)

#                 e2 = new_graph.edge(min_edge.target(), n)
#                 new_graph.remove_edge(e2)
#                 e1 = new_graph.add_edge(min_edge.source(), n)
#                 new_graph.edge_properties["e_weights"][e1] = final_weights_appr
#                 new_graph.edge_properties["e_label"][e1] = new_edge_name
#                 new_graph.edge_properties['e_entropy'][e1] = final_weights_appr['entropy']
            
#             # to_be_removed_nodes += [min_edge.target()]
#             to_be_removed_nodes.append(min_edge.target())
#             new_graph.remove_edge(min_edge)
#             edges.remove(min_edge)

#     new_graph.remove_vertex(to_be_removed_nodes)

#     e_labels_ch = new_graph.edge_properties["e_label"]
#     v_labels_ch = new_graph.vertex_properties["v_label"]
        
#     v_label_reversed = {}
#     for v in new_graph.vertices():
#         v_label = v_labels_ch[v]
#         v_label_reversed[v_label] = int(v)
    
#     e_label_reversed = {}
#     for e in new_graph.edges():
#         e_label = e_labels_ch[e]
#         v1_label, v2_label = e_label.split('--')
#         e_label_reversed[e_label] = [v_label_reversed[v1_label], v_label_reversed[v2_label]]
        

#     this_path = os.path.join(chordal_graph_path, frag_file.split('.')[0] + '.gt.gz')
#     new_graph.save(this_path)

#     with open(chordal_v_label_revered_path, 'wb') as f:
#         pickle.dump(v_label_reversed, f)

#     with open(chordal_e_label_revered_path, 'wb') as f:
#         pickle.dump(e_label_reversed, f)

#     print('[Done]', this_path)


# def run_FFBS_quotient(inp):
#     this_frag_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
#     print('Working on:', os.path.join(this_frag_path, frag_file))
#     # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/test.frag'
#     # frag_path = '/labs/Aguiar/pHapCompass/test/test2.frag'
#     # ploidy= 3
#     # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/haplotypes.csv'
#     # genotype_path = '/labs/Aguiar/pHapCompass/test/haplotypes.csv'

#     class Args:
#         def __init__(self):
#             self.vcf_path = 'example/62_ID0.vcf'
#             self.data_path = os.path.join(this_frag_path, frag_file)
#             # self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_90.frag.txt'
#             self.bam_path = 'example/example.bam'
#             self.genotype_path = genotype_path
#             self.ploidy = ploidy
#             self.error_rate = 0.001
#             self.epsilon = 0.0001
#             self.output_path = 'output'
#             self.root_dir = 'D:/UCONN/HaplOrbit'
#             self.alleles = [0, 1]

#     # Create the mock args object
#     args = Args()

#     # Initialize classes with parsed arguments
#     input_handler = InputHandler(args)

#     config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

#     # frag_graph_path = os.path.join(this_frag_graph_path, frag_file.split('.')[0] + '.gt.gz')
#     # frag_graph = gt.load_graph(frag_graph_path)

#     fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
#     fragment_model.construct(input_handler, config)

#     quotient_g_path = os.path.join(this_quotient_coverage_path, frag_file.split('.')[0] + '.gt.gz')
#     quotient_g = gt.load_graph(quotient_g_path)

#     edge_map_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + frag_file.split('.')[0] + '.pkl')
#     with open(edge_map_path, 'rb') as f:
#         edges_map_quotient = pickle.load(f)
        
#     quotient_g_v_label_reversed_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + frag_file.split('.')[0] + '.pkl')
#     with open(quotient_g_v_label_reversed_path, 'rb') as f:
#         quotient_g_v_label_reversed = pickle.load(f)

#     start_time = time.time()

#     transitions_dict, transitions_dict_extra = transition_matrices_v2(quotient_g, edges_map_quotient, ploidy, config, fragment_model)
#     emission_dict = emissions_v2(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

#     nodes = list(emission_dict.keys())
#     edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

#     slices, _ =  assign_slices_and_interfaces(nodes, edges)

#     assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, args.data_path)

#     forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, args.data_path)

#     # backward_messages = compute_backward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, args.data_path)

#     samples = sample_states_book(slices, edges, forward_messages, transitions_dict)
#     # samples = sample_states(slices, edges, forward_messages, transitions_dict)
#     # for k in samples.keys():
#     #     kedges = samples[k].keys()
#     #     for e in kedges:
#     #         print(e, samples[k][e])
    
#     predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config)
#     end_time = time.time()
#     elapsed_time = round(end_time - start_time, 2)


#     # print('Predicted Haplotypes:\n', predicted_haplotypes)
#     # print('\nTrue Haplotypes:\n', pd.read_csv(genotype_path).T)
#     sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]

#     predicted_haplotypes_np = predicted_haplotypes[sampled_positions].to_numpy()
#     true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()[:, sampled_positions]

#     vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes)
#     accuracy, _ = calculate_accuracy(predicted_haplotypes_np, true_haplotypes)
#     mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes)
#     mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)
#     results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
#     results = {}
#     results['evaluation'] = {'vector_error_rate': vector_error_rate, 'vector_error': vector_error, 'backtracking_steps': backtracking_steps, 
#                              'dp_table': dp_table, 'accuracy': accuracy, 'mismatch_error': mismatch_error, 'mec': mec_}
#     results['predicted_haplotypes'] = predicted_haplotypes_np
#     results['true_haplotypes'] = true_haplotypes
#     results['forward_messages'] = forward_messages
#     # results['backward_messages'] = backward_messages
#     results['transitions_dict'] = transitions_dict
#     results['transitions_dict_extra'] = transitions_dict_extra
#     results['emission_dict'] = emission_dict
#     results['assignment_dict'] = assignment_dict
#     results['samples'] = samples
#     results['slices'] = slices
#     results['best_permutation'] = best_permutation
#     results['fragment_list'] = fragment_model.fragment_list
#     results['time'] = elapsed_time

#     with open(os.path.join(results_path, results_name), 'wb') as f:
#         pickle.dump(results, f)

#     print('Saved results in {}.'.format(os.path.join(results_path, results_name)), 'vector_error_rate', vector_error_rate, 'accuracy', accuracy, 'mismatch_error', mismatch_error, 'mec', mec_)
    

# def run_FFBS_quotient_count(inp):
#     this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
#     file_name = frag_file.split('.')[0]

#     fragment_v_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_v_label_' + file_name + '.pkl')
#     fragment_e_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_e_label_' + file_name + '.pkl')
#     quotient_v_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + file_name + '.pkl')
#     quotient_e_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + file_name + '.pkl')
#     print('Working on:', os.path.join(this_frag_path, frag_file))

#     class Args:
#         def __init__(self):
#             self.vcf_path = 'example/62_ID0.vcf'
#             self.data_path = os.path.join(this_frag_path, frag_file)
#             # self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_90.frag.txt'
#             self.bam_path = 'example/example.bam'
#             self.genotype_path = genotype_path
#             self.ploidy = ploidy
#             self.error_rate = 0.001
#             self.epsilon = 0.0001
#             self.output_path = 'output'
#             self.root_dir = 'D:/UCONN/HaplOrbit'
#             self.alleles = [0, 1]

#     # Create the mock args object
#     args = Args()

#     start_time = time.time()

#     # Initialize classes with parsed arguments
#     input_handler = InputHandler(args)

#     config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

#     fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
#     fragment_model.construct(input_handler, config)

#     frag_graph_path = os.path.join(this_fragment_coverage_path, file_name + '.gt.gz')
#     fragment_model.graph.save(frag_graph_path)

#     with open(fragment_v_label_revered_path, "wb") as f:
#         pickle.dump(fragment_model.v_label_reversed, f)

#     edges_map_fragment = {}
#     for k in fragment_model.e_label_reversed.keys():
#         edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

#     with open(fragment_e_label_revered_path, "wb") as f:
#         pickle.dump(edges_map_fragment, f)

#     # create quotient graph
#     quotient_g = QuotientGraph(fragment_model)
#     quotient_g.construct(input_handler, config)

#     # save quotient graph
#     quot_graph_path = os.path.join(this_quotient_coverage_path, file_name + '.gt.gz')
#     quotient_g.graph.save(quot_graph_path)

#     with open(quotient_v_label_revered_path, "wb") as f:
#         pickle.dump(quotient_g.v_label_reversed, f)

#     edges_map_quotient = {}
#     for k in quotient_g.e_label_reversed.keys():
#         edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

#     with open(quotient_e_label_revered_path, "wb") as f:
#         pickle.dump(edges_map_quotient, f)

#     quotient_g_v_label_reversed = quotient_g.v_label_reversed

#     edges_map_quotient = {}
#     for k in quotient_g.e_label_reversed.keys():
#         edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

#     transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config)
#     emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

#     nodes = list(emission_dict.keys())
#     edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

#     slices, interfaces =  assign_slices_and_interfaces(nodes, edges)

#     assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, input_handler.data_path)

#     forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)

#     samples = sample_states_book(slices, edges, forward_messages, transitions_dict)

#     predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="counts")

#     end_time = time.time()

#     elapsed_time = round(end_time - start_time, 2)

#     sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]

#     predicted_haplotypes_np = predicted_haplotypes[sampled_positions].to_numpy()
#     true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()[:, sampled_positions]

#     vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes)
#     accuracy, _ = calculate_accuracy(predicted_haplotypes_np, true_haplotypes)
#     mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes)
#     mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)
#     results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
#     results = {}
#     results['evaluation'] = {'vector_error_rate': vector_error_rate, 'vector_error': vector_error, 'backtracking_steps': backtracking_steps, 
#                              'dp_table': dp_table, 'accuracy': accuracy, 'mismatch_error': mismatch_error, 'mec': mec_}
#     results['predicted_haplotypes'] = predicted_haplotypes_np
#     results['true_haplotypes'] = true_haplotypes
#     results['forward_messages'] = forward_messages
#     results['transitions_dict'] = transitions_dict
#     results['transitions_dict_extra'] = transitions_dict_extra
#     results['emission_dict'] = emission_dict
#     results['assignment_dict'] = assignment_dict
#     results['samples'] = samples
#     results['slices'] = slices
#     results['best_permutation'] = best_permutation
#     results['fragment_list'] = fragment_model.fragment_list
#     results['time'] = elapsed_time
#     # print('Results:', results['evaluation'])

#     with open(os.path.join(results_path, results_name), 'wb') as f:
#         pickle.dump(results, f)

#     print('Saved results in {}.'.format(os.path.join(results_path, results_name)), 'vector_error_rate', vector_error_rate, 'accuracy', accuracy, 'mismatch_error', mismatch_error, 'mec', mec_)
    

# def run_FFBS_quotient_likelihood(inp):
#     this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
#     file_name = frag_file.split('.')[0]

#     fragment_v_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_v_label_' + file_name + '.pkl')
#     fragment_e_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_e_label_' + file_name + '.pkl')
#     quotient_v_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + file_name + '.pkl')
#     quotient_e_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + file_name + '.pkl')
#     print('Working on:', os.path.join(this_frag_path, frag_file))

#     class Args:
#         def __init__(self):
#             self.vcf_path = 'example/62_ID0.vcf'
#             self.data_path = os.path.join(this_frag_path, frag_file)
#             # self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_90.frag.txt'
#             self.bam_path = 'example/example.bam'
#             self.genotype_path = genotype_path
#             self.ploidy = ploidy
#             self.error_rate = 0.001
#             self.epsilon = 0.0001
#             self.output_path = 'output'
#             self.root_dir = 'D:/UCONN/HaplOrbit'
#             self.alleles = [0, 1]

#     # Create the mock args object
#     args = Args()

#     start_time = time.time()

#     # Initialize classes with parsed arguments
#     input_handler = InputHandler(args)

#     config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

#     fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
#     fragment_model.construct(input_handler, config)

#     frag_graph_path = os.path.join(this_fragment_coverage_path, file_name + '.gt.gz')
#     fragment_model.graph.save(frag_graph_path)

#     frag_graph_plot_path = os.path.join(this_fragment_coverage_path, file_name + '.png')

#     e_labels = fragment_model.graph.edge_properties["e_label"]
#     v_labels = fragment_model.graph.vertex_properties["v_label"]
#     gt.graph_draw(fragment_model.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
#     edge_font_size=10, output=frag_graph_plot_path)

#     with open(fragment_v_label_revered_path, "wb") as f:
#         pickle.dump(fragment_model.v_label_reversed, f)

#     edges_map_fragment = {}
#     for k in fragment_model.e_label_reversed.keys():
#         edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

#     with open(fragment_e_label_revered_path, "wb") as f:
#         pickle.dump(edges_map_fragment, f)

#     # create quotient graph
#     quotient_g = QuotientGraph(fragment_model)
#     quotient_g.construct(input_handler, config)

#     # save quotient graph
#     quot_graph_path = os.path.join(this_quotient_coverage_path, file_name + '.gt.gz')
#     quotient_g.graph.save(quot_graph_path)


#     quotient_graph_plot_path = os.path.join(this_quotient_coverage_path, file_name + '.png')
#     e_labels = quotient_g.graph.edge_properties["e_label"]
#     v_labels = quotient_g.graph.vertex_properties["v_label"]
#     gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
#     edge_font_size=10, output=quotient_graph_plot_path)

#     # gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
#     # edge_font_size=10)

#     with open(quotient_v_label_revered_path, "wb") as f:
#         pickle.dump(quotient_g.v_label_reversed, f)

#     edges_map_quotient = {}
#     for k in quotient_g.e_label_reversed.keys():
#         edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

#     with open(quotient_e_label_revered_path, "wb") as f:
#         pickle.dump(edges_map_quotient, f)

#     quotient_g_v_label_reversed = quotient_g.v_label_reversed

#     edges_map_quotient = {}
#     for k in quotient_g.e_label_reversed.keys():
#         edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

#     transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config)
#     emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

#     nodes = list(emission_dict.keys())
#     edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

#     slices, interfaces =  assign_slices_and_interfaces(nodes, edges)

#     assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, input_handler.data_path)

#     forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)

#     # backward_messages = compute_backward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)   


#     # samples = sample_states_book(slices, edges, forward_messages, transitions_dict)
#     samples = sample_states_book_multiple_times(slices, edges, forward_messages, transitions_dict, n=100)
#     # samples = sample_states_ground_truth(slices, nodes, genotype_path)
#     # fragment_list = fragment_model.fragment_list
#     # reads_dict = calculate_pair_counts(fragment_list)

#     ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, samples, ploidy)
#     # print('FFBS Accuracy:', ffbs_acc)
#     predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")

#     # for _ in range(10):
#     #     samples = sample_states_book(slices, edges, forward_messages, transitions_dict)
#     #     predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")
#     #     ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, samples, ploidy)
#     #     this_phasing_likelihood = compute_global_phasing_likelihood(predicted_haplotypes, fragment_model, config)
#     #     print(ffbs_acc, this_phasing_likelihood)
#     #     compute_vector_error_rate(predicted_haplotypes.to_numpy(), true_haplotypes.to_numpy())

#     end_time = time.time()

#     elapsed_time = round(end_time - start_time, 2)

#     true_haplotypes = pd.read_csv(genotype_path).T

#     block_info, components = get_block_info(quotient_g, predicted_haplotypes, true_haplotypes, fragment_model)

#     sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]

#     predicted_haplotypes_np = predicted_haplotypes[sampled_positions].to_numpy()
#     # true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()[:, sampled_positions]
    
#     true_haplotypes_np = true_haplotypes.to_numpy()[:, sampled_positions]

#     vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes_np)
#     accuracy, _ = calculate_accuracy(predicted_haplotypes_np, true_haplotypes_np)
#     mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes_np)
#     mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)
#     results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
#     results = {}
#     results['block_evaluation'] = block_info
#     results['components'] = components
#     results['n_blocks'] = len(components.keys())
#     results['average_block_size'] = block_info['average_block_size']
#     results['length_phased'] = len(sampled_positions)
#     results['evaluation'] = {'vector_error_rate': vector_error_rate, 'vector_error': vector_error, 'backtracking_steps': backtracking_steps, 
#                              'dp_table': dp_table, 'accuracy': accuracy, 'mismatch_error': mismatch_error, 'mec': mec_, 'ffbs_acc': ffbs_acc}
#     results['predicted_haplotypes'] = predicted_haplotypes
#     results['true_haplotypes'] = pd.read_csv(genotype_path).T
#     results['forward_messages'] = forward_messages
#     results['transitions_dict'] = transitions_dict
#     results['transitions_dict_extra'] = transitions_dict_extra
#     results['emission_dict'] = emission_dict
#     results['assignment_dict'] = assignment_dict
#     results['samples'] = samples
#     results['slices'] = slices
#     results['best_permutation'] = best_permutation
#     results['fragment_list'] = fragment_model.fragment_list
#     results['time'] = elapsed_time
#     # print('Results:', results['evaluation'])

#     with open(os.path.join(results_path, results_name), 'wb') as f:
#         pickle.dump(results, f)

#     print('Saved results in {}.'.format(os.path.join(results_path, results_name)), 'vector_error_rate', vector_error_rate, 'vector_error', vector_error, 'mismatch_error', mismatch_error, 'mec', mec_, 'ffbs_acc', ffbs_acc)
    

def simulate_na12878():

    beagle_config_human = {
        "snp_df_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv',
        # "input_vcf_path": '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/hapref_chr21_filtered.vcf.bgz',
        "input_vcf_path": '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/updated_NA12878_extracted.vcf',
        "contig_fasta": '/mnt/research/aguiarlab/data/hg19/chr21.fa',
        "main_path": '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878',
        "art_path": 'art_illumina',
        "extract_hairs_path": 'extractHAIRS',
        "n_samples": 2, 
        "target_spacing": 100,
        "densify_snps": False, 
        "contig_lens": [100], 
        "ploidies": [3],
        "coverages": [10],
        "read_length": 150,
        "mean_insert_length": 2000,
        "std_insert_length": 300
        }


    simulator = Simulator(beagle_config_human)
    # simulator.generate_genomes_fasta()

    # simulator.simulate()
    
    inputs = make_inputs_for_generate_qoutient_graph(simulator)

    pool = Pool(2)
    pool.map(generate_quotient_graph, inputs)

    # for inp in inputs:
    #     print(inp[4])
    #     generate_quotient_graph(inp)


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

    next_inputs = make_inputs_for_run_likelihood(simulator)
    print('number of inputs:', len(next_inputs))
    # pool = Pool(10)
    # pool.map(run_FFBS_quotient_likelihood, next_inputs)
    for inp in next_inputs:
        print(inp[4])
        run_FFBS_quotient_likelihood(inp)


# def save_inputs(inputs, output_dir):
#     """
#     Save each input as a separate pickle file in the specified output directory.
#     """
#     output_dir = '/mnt/research/aguiarlab/proj/HaplOrbit/inputs100_2'
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     for i, inp in enumerate(inputs):
#         input_file = os.path.join(output_dir, f"input_{i}.pkl")
#         with open(input_file, "wb") as f:
#             pickle.dump(inp, f)
#     print(f"Saved {len(inputs)} inputs to {output_dir}")


# def run_FFBS_quotient_likelihood_from_input(input_file):
#     with open(input_file, "rb") as f:
#         inp = pickle.load(f)

#     run_FFBS_quotient_likelihood(inp)


# if __name__ == '__main__':

#     # simulate_na12878()
#     # simulate_awri()

#     if len(sys.argv) != 2:
#         print("Usage: python3 simulator_paper.py <input_file>")
#         sys.exit(1)
    
#     input_file = sys.argv[1]
#     run_FFBS_quotient_likelihood_from_input(input_file)