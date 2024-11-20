import os
import pandas as pd
import random
import pysam



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


    def generate_genomes_fasta(self):
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
                            genomes[i][pos] = alts[0] if genotype[i] == 1 else ref
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
                                genomes[i][pos] = alts[0] if genotype[i] == 1 else ref
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


    def simulate_fastq_art(self):
        """
        Simulate FASTQ files using ART.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                fasta_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'contig_{contig_len}_ploidy_{ploidy}.fa')
                this_sh_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'01_simulate_{contig_len}_{ploidy}.sh')

                to_print = self.get_slurm_header('fastq')
                for coverage in self.coverages:
                    cov_path = os.path.join(self.main_path, f'contig_{contig_len}', f'ploidy_{ploidy}', f'cov_{coverage}')
                    if not os.path.exists(cov_path):
                        os.makedirs(cov_path, exist_ok=True)
                    fastq_path = os.path.join(cov_path, 'fastq')
                    if not os.path.exists(fastq_path):
                        os.makedirs(fastq_path, exist_ok=True)

                    for rd in range(self.n_samples):
                        command = f'{self.art_path} -ss HS25 -i {fasta_path} -p -na -l {self.read_length} -f {coverage} -m {self.mil} -s {self.sil} -o {fastq_path}/{str(rd).zfill(2)}\n'
                        to_print += command

                with open(this_sh_path, 'w') as f:
                    f.write(to_print)


    def align_fastq_files(self):
        """
        Align simulated FASTQ files using BWA.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                fasta_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'contig_{}_ploidy_{}.fa'.format(contig_len, ploidy))
                fa_index = 'bwa index {}\n\n'.format(fasta_path)
                this_sh_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '02_align_{}_{}.sh'.format(contig_len, ploidy))
                to_print = self.get_slurm_header('align')
                to_print += 'module load bwa-mem2/2.1\nmodule load bwa/0.7.17\n' + fa_index
                for coverage in self.coverages:
                    this_cov_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                    fastq_path = os.path.join(this_cov_path, 'fastq')
                    bam_path = os.path.join(this_cov_path, 'bam')
                    if not os.path.exists(bam_path):
                        os.makedirs(bam_path)
                    for rd in range(self.n_samples):
                        command = 'bwa mem {} {}/{}1.fq {}/{}2.fq > {}/{}.sam\n'.format(fasta_path, fastq_path, str(rd).zfill(2), fastq_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
                        sort_com = 'samtools view -Sb {}/{}.sam | samtools sort -o {}/{}.bam\n'.format(bam_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
                        index_com = 'samtools index {}/{}.bam\n'.format(bam_path, str(rd).zfill(2))
                        rmv_com = 'rm {}/{}.sam\n\n'.format(bam_path, str(rd).zfill(2))
                        to_print += command
                        to_print += sort_com
                        to_print += index_com
                        to_print += rmv_com
                with open(this_sh_path, 'w') as f:
                    f.write(to_print)


    def extract_hairs(self):
        """
        Extract haplotypes using ExtractHAIRS.
        """
        for contig_len in self.contig_lens:
            for ploidy in self.ploidies:
                this_sh_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '03_extract_hairs_{}_{}.sh'.format(contig_len, ploidy))
                to_print = self.get_slurm_header('exttract_hairs')
                to_print += 'module load bcftools/1.20\nmodule load htslib/1.20\n\n'
                this_vcf_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'maf0.01_hapref_chr21_filtered_NA12878_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
                this_vcf_sorted_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'sorted_maf0.01_hapref_chr21_filtered_NA12878_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
                unzipped_vcf = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'sorted_maf0.01_hapref_chr21_filtered_NA12878_ploidy{}_contig{}.vcf'.format(ploidy, contig_len))
                sort_cmd = 'bcftools sort {} -Oz -o {}\n'.format(this_vcf_path, this_vcf_sorted_path)
                idx_cmd = 'tabix -p vcf {}\n'.format(this_vcf_sorted_path)
                unzip_cmd = 'gunzip {}\n\n'.format(this_vcf_sorted_path)
                to_print += sort_cmd
                to_print += idx_cmd
                to_print += unzip_cmd
                for coverage in self.coverages:
                    this_cov_path = os.path.join(self.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                    bam_path = os.path.join(this_cov_path, 'bam')
                    frag_path = os.path.join(this_cov_path, 'frag')
                    if not os.path.exists(frag_path):
                        os.makedirs(frag_path)
                    for rd in range(self.n_samples): 
                        command = '{} --bam {}/{}.bam --noquality --vcf {} --out {}/{}.frag\n'.format(self.extract_hairs_path, bam_path, str(rd).zfill(2), unzipped_vcf, frag_path, str(rd).zfill(2))
                        to_print += command
                with open(this_sh_path, 'w') as f:
                    f.write(to_print)


    def simulate(self):
        """
        Run the entire simulation pipeline.
        """
        self.generate_genomes_fasta()
        self.simulate_fastq_art()
        self.align_fastq_files()
        self.extract_hairs()



if __name__ == '__main__':

    config = {
        "snp_df_path": '/labs/Aguiar/pHapCompass/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv',
        "input_vcf_path": '/labs/Aguiar/pHapCompass/simulated_data_NEW/hapref_chr21_filtered.vcf.bgz',
        "contig_fasta": '/labs/Aguiar/pHapCompass/references/EGA.GRCh37/chr21.fa',
        "art_path": '/labs/Aguiar/pHapCompass/ART/art_bin_MountRainier/art_illumina',
        "main_path": '/labs/Aguiar/pHapCompass/simulated_data_long',
        "extract_hairs_path": '/home/FCAM/mhosseini/HaplOrbit/extract_poly/build/extractHAIRS',
        "n_samples": 3, 
        "target_spacing": 100,
        "densify_snps": True, 
        "contig_lens": [100], 
        "ploidies": [3],
        "coverages": [10],
        "read_length": 150,
        "mean_insert_length": 800,
        "std_insert_length": 150
        }

    simulator = Simulator(config)
    simulator.simulate()