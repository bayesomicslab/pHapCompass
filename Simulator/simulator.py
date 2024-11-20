import numpy as np
import pandas as pd
import random
import os
from functools import reduce
from multiprocessing import Pool


def find_NA12878_column():
    file_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.1/windows/50000/hapref_chr21_filtered.vcf.bgz_sample9_len50000.bcf'
    df = pd.read_csv(file_path, skiprows=48, nrows=2, sep='\t')
    columns = list(df.columns.values)
    index = columns.index('NA12878')

    # Father: NA12891
    # Mother: NA12892
    # child: NA12878

def extract_column_NA12878():
    file_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.1/windows/50000/sample_NA12878.txt'
    df = pd.read_csv(file_path, skiprows=48, sep=' ')
    hetero_df = df[df['NA12878'].isin(['1|0', '0|1'])].reset_index(drop=True)
    hetero_df = hetero_df.sort_values(by='POS').reset_index(drop=True)
    hetero_df['sim_pos'] = hetero_df.index
    return hetero_df


def swap_fragments_based_on_start(s):
    # Split the string into parts
    parts = s.split()

    # Extract the two numbers and their associated sections
    num1, zeros1, num2, zeros2 = parts

    # Compare the two numbers as integers
    if int(num1) > int(num2):
        # Swap the numbers along with their associated 0s/1s sections
        return f"{num2} {zeros2} {num1} {zeros1} "
    else:
        # Return the original string if no swap is needed
        return s


def flip_with_error_rate(data, error_rate=0.001):
    # Iterate over the list and flip each value with the given error rate
    return [1 - x if random.random() < error_rate else x for x in data]


def gen_coverages(expected_avg, n, a=3, b=9):
    while True:
        l = [random.randint(a, b) for i in range(n)]
        avg = reduce(lambda x, y: x + y, l) / len(l)

        if avg == expected_avg:
            return l
        

def generate_reads(inp):
    file_path, df, cov, min_read_length, max_read_length = inp
    
    to_print = ''  # To store all generated reads
    position_coverage = np.zeros(len(df))  # Track how many reads cover each position
    # position_coverage_limit = gen_coverages(cov, len(df), a=3, b=9)
    read_id = 0
    # Step 2: Loop until every position is covered by at least `min_reads_per_pos` reads
    # while (position_coverage < position_coverage_limit).any():

    while reduce(lambda x, y: x + y, position_coverage) / len(df) <= cov:
        # print(reduce(lambda x, y: x + y, position_coverage) / len(df))
        n_frag = np.random.choice([1, 2], p=[0.95, 0.05])
        fragment_line = '{} '.format(n_frag)
        read_info = ''
        allels = ''
        for _ in range(n_frag):
            # print(nf)
            # Randomly pick a start position for the new read
            start = np.random.randint(0, len(df) - min_read_length + 1)

            # Randomly decide the read length (between min_read_length and max_read_length)
            read_length = np.random.randint(min_read_length, max_read_length + 1)

            # Ensure the read does not exceed the dataframe's length
            end = min(start + read_length, len(df))

            # Store the read (subset of the dataframe)
            read = df.iloc[start:end]
            # reads.append(read)
            
            haplotype_no = np.random.choice([0, 2])
            hap_dict = {0: 1, 2: 2}
            read_major = list(read['REF'].values)
            read_minor = list(read['ALT'].values)
            
            haplotypes = [int(elem[haplotype_no]) for elem in list(read['NA12878'].values)]
            random_hap = flip_with_error_rate(haplotypes, error_rate=0.001)
            hap = ''.join([str(i) for i in random_hap])
            allels += ''.join([read_major[elem_id] if elem == 0 else read_minor[elem_id] for elem_id, elem in enumerate(random_hap)])
            read_info += '{} {} '.format(start, hap)
            # print(read_info)
            # Update the coverage for each position in this read
            position_coverage[start:end] += 1
        # print(len(read_info.split(' ')), read_info.split(' '))
        if len(read_info.split(' ')) == 5:
            read_info = swap_fragments_based_on_start(read_info)
            # print(read_info)

        fragment_line += 'h{}_NA12878_{} {}{}'.format(hap_dict[haplotype_no], f"{(read_id):06d}", read_info, allels)
        to_print += '{}\n'.format(fragment_line)
        read_id += 1
    with open(file_path, 'w') as f:
        f.write(to_print)



def check_compatibility_vcf_fasta():
    vcf_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/_EGAZ00001239288_HRC.r1-1.EGA.GRCh37.chr21.haplotypes.vcf.gz'
    vcf_df = pd.read_csv(vcf_path, skiprows=48, nrows=10, sep='\t')
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


def generate_short_reads():
    inputs = []
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data'
    min_read_length=2
    max_read_length=3
    contigs = ['Contig1_k2']
    contigs_dict = {'Contig1_k2': {'start': 0, 'end': 1000}}
    coverages = [6]
    df = extract_column_NA12878()
    df.to_csv(os.path.join(main_path, 'complete_data.csv'), index=False)

    for cont in contigs:
        contig_path = os.path.join(main_path, cont)
        if not os.path.exists(contig_path):
            os.makedirs(contig_path)

        this_df = df.iloc[contigs_dict[cont]['start']: contigs_dict[cont]['end']].reset_index()
        this_df.to_csv(os.path.join(contig_path, 'contig_vcf_data.csv'), index=False)
        real_hap = this_df[['NA12878']]
        real_hap[['h1', 'h2']] = real_hap['NA12878'].str.split('|', expand=True)
        real_hap = real_hap[['h1', 'h2']]
        real_hap.to_csv(os.path.join(contig_path, 'real_haps_{}.txt'.format(cont.lower())), sep='\t', index=False, header=False)

        for cov in coverages:
            this_save_path = os.path.join(main_path, cont, 'c' + str(cov))
            if not os.path.exists(this_save_path):
                os.makedirs(this_save_path)

            for sample in range(100):
                sample_path = os.path.join(this_save_path, 'SIM_' + str(sample) + '.frag.txt')
                inputs.append([sample_path, this_df, cov, min_read_length, max_read_length])
    
    
    pool = Pool(20)
    pool.map(generate_reads, inputs)




def generate_genomes_fasta():
     # snp_df_NA12878_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv'
    snp_df_NA12878_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv'
    snp_df_NA12878 = pd.read_csv(snp_df_NA12878_path)
    snp_df_NA12878['diff'] = snp_df_NA12878['POS'].diff()
    for i in range(len(snp_df_NA12878)):
        print(snp_df_NA12878.loc[i, 'POS'], snp_df_NA12878.loc[i, 'diff'])
    
    snp_positions = list(snp_df_NA12878['POS'].values)
    # input_vcf_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/hapref_chr21_filtered.vcf.bgz'
    input_vcf_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW/hapref_chr21_filtered.vcf.bgz'
    # contig_fasta = '/mnt/research/aguiarlab/data/hg19/chr21.fa'
    contig_fasta = '/labs/Aguiar/pHapCompass/references/EGA.GRCh37/chr21.fa'
    # output_fasta_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW'
    output_fasta_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW'
    is_phased = True
    # output_vcf = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.vcf.gz'
    with open(contig_fasta) as f:
        contig_name = f.readline().strip()[1:]  # Get contig name without '>'
        contig_seq = f.read().replace("\n", "")  # Get sequence without newlines

    contig_lens = [100] # [10, 100, 300, 1000] 
    ploidies = [3] #, 4, 6, 8, 10]
    for contig_len in contig_lens:
        contig_path = os.path.join(output_fasta_path, 'contig_{}'.format(contig_len))
        if not os.path.exists(contig_path):
            os.makedirs(contig_path)
        for ploidy in ploidies:
            contig_ploidy_path = os.path.join(contig_path, 'ploidy_{}'.format(ploidy))
            if not os.path.exists(contig_ploidy_path):
                os.makedirs(contig_ploidy_path)
            



            snp_counter = 0
            last_snp = snp_positions[0:contig_len][-1]
            output_vcf = os.path.join(contig_ploidy_path, 'maf0.01_hapref_chr21_filtered_NA12878_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
            haplotypes = [[] for _ in range(ploidy)]  # Initialize ploidy haplotype sequences

            # Initialize ploidy genomes with the reference sequence
            genomes = [list(contig_seq[:last_snp]) for _ in range(ploidy)]  # List of lists for mutable strings

            vcf_in = pysam.VariantFile(input_vcf_path)

            # Create a new header for the output VCF file with only NA12878
            new_header = pysam.VariantHeader()
            for line in vcf_in.header.records:
                new_header.add_line(str(line))  # Copy metadata lines

            # Add custom contigs for haplotypes
            for k in range(ploidy):
                new_header.contigs.add(f"haplotype_{k + 1}", length=len(contig_seq))

            # Add the sample
            new_header.add_sample("NA12878")  # Add only NA12878

            # Open the output VCF file with the new header
            vcf_out = pysam.VariantFile(output_vcf, "wz", header=new_header)

            for record in vcf_in.fetch():
                # print(record)
                pos = record.pos - 1  # Convert 1-based VCF position to 0-based for indexing
                ref = record.ref
                alts = record.alts[:1]  # Only use the first alternative allele

                # Generate a random heterozygous genotype for ploidy = 3
                if len(alts) > 0 and pos + 1 in snp_positions and snp_counter < contig_len:
                    # print(pos)
                    snp_counter += 1
                    print(snp_counter)
                    # stop
                    g = random.randint(1, ploidy - 1)
                    selected_genomes = random.sample(range(ploidy), g)
                    genotype = tuple(1 if i in selected_genomes else 0 for i in range(ploidy))

                    # Assign alleles in the genomes based on the selected genotype
                    for i in range(ploidy):
                        genomes[i][pos] = alts[0] if genotype[i] == 1 else ref
                        haplotypes[i].append(genotype[i])

                    # Format the genotype as VCF-style (e.g., 0|1|0)
                    # genotype_str = "|".join(map(str, genotype))
                    qual = 75 # round(random.uniform(3, 75), 4)
                    
                    for k in range(ploidy):
                        # Create a new VCF record with only NA12878
                        new_record = vcf_out.new_record(
                            # contig=record.chrom,
                            contig='haplotype_{}'.format(str(k + 1)),
                            start=record.start,
                            stop=record.stop,
                            alleles=(ref, alts[0]),
                            id=record.id,
                            qual=qual,
                            filter=record.filter.keys(),
                            info=record.info
                        )
                        new_record.samples["NA12878"]["GT"] = genotype
                        new_record.samples["NA12878"].phased = is_phased
                        vcf_out.write(new_record)

            vcf_out.close()

            output_fasta = os.path.join(contig_ploidy_path, f'contig_{contig_len}_ploidy_{ploidy}.fa')
            with open(output_fasta, "w") as f_out:
                for i in range(ploidy):
                    f_out.write(f">haplotype_{i + 1}\n")
                    f_out.write("".join(genomes[i]) + "\n")

            # write haplotypes to a DataFrame
            haplotype_df = pd.DataFrame(haplotypes).T  # Transpose so each column represents a haplotype

            # Rename columns for clarity
            haplotype_df.columns = [f"haplotype_{i + 1}" for i in range(len(haplotypes))]

            # Save the DataFrame as a CSV if needed
            haplotype_df.to_csv(os.path.join(contig_ploidy_path, 'haplotypes.csv'), index=False)


def simulate_fastq_art():
    # main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW'
    main_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW'
    # art_path = '/home/mah19006/downloads/art_bin_MountRainier/art_illumina'
    art_path = '/labs/Aguiar/pHapCompass/ART/art_bin_MountRainier/art_illumina'
    # sh_path = os.path.join('/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation', '01_simulate_illumina.sh')
    sh_path = os.path.join('/labs/Aguiar/pHapCompass/scripts/simulation', '01_simulate_illumina.sh')
    
    # to_print = ''
    contig_lens = [100] #[10, 100, 300, 1000]
    ploidies = [3] #[3, 4, 6, 8, 10]
    coverages = [10] #, 20, 30, 40, 50]
    mil = 400
    sil = 50
    read_length = 150
    for contig_len in contig_lens:
        for ploidy in ploidies:
            fasta_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'contig_{}_ploidy_{}.fa'.format(contig_len, ploidy))
            this_sh_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '01_simulate_{}_{}.sh'.format(contig_len, ploidy))
            to_print = ''
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                if not os.path.exists(this_cov_path):
                    os.makedirs(this_cov_path)
                fastq_path = os.path.join(this_cov_path, 'fastq')
                if not os.path.exists(fastq_path):
                    os.makedirs(fastq_path)
                for rd in range(100):
                    command = '{} -ss HS25 -i {} -p -na -sam -l {} -f {} -m {} -s {} -o {}/{}\n'.format(art_path, fasta_path, read_length, coverage, mil, sil, fastq_path, str(rd).zfill(2))
                    to_print += command
                    # single_cmd = './art_454 -s -t -r {} {} {}/{}_single {}\n'.format(rn, chr_reference_path, this_sim_path, str(rd).zfill(2), coverage)
                    # paired_cmd = './art_454 -s -t -r {} {} {}/{}_paired {} {} {}\n\n'.format(rn, chr_reference_path, this_sim_path, str(rd).zfill(2), coverage, mil, sil)
                    # to_print += single_cmd
                    # to_print += paired_cmd
            with open(this_sh_path, 'w') as f:
                f.write(to_print)


def align_fastq_files():
    main_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW'
    chr21_fasta_path = '/labs/Aguiar/pHapCompass/references/EGA.GRCh37/chr21.fa'
    sh_path = os.path.join('/labs/Aguiar/pHapCompass/scripts/simulation', '02_align.sh')
    # to_print = 'module load bwa-mem2/2.1\nmodule load bwa/0.7.17\n\n'
    contig_lens = [100] # [10, 100, 300, 1000]
    ploidies =  [3] # [3, 4, 6, 8, 10]
    coverages = [10] # [10, 20, 30, 40, 50]
    for contig_len in contig_lens:
        for ploidy in ploidies:
            fasta_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'contig_{}_ploidy_{}.fa'.format(contig_len, ploidy))
            fa_index = 'bwa index {}\n\n'.format(fasta_path)
            this_sh_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '02_align_{}_{}.sh'.format(contig_len, ploidy))
            # to_print += fa_index
            to_print = 'module load bwa-mem2/2.1\nmodule load bwa/0.7.17\n' + fa_index
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                fastq_path = os.path.join(this_cov_path, 'fastq')
                bam_path = os.path.join(this_cov_path, 'bam')
                if not os.path.exists(bam_path):
                    os.makedirs(bam_path)
                for rd in range(100):
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

        # '/home/FCAM/mhosseini/HaplOrbit/extract_poly/build/extractHAIRS --bam /labs/Aguiar/pHapCompass/simulated_data_NEW/contig_100/ploidy_3/cov_10/bam/00.bam --vcf /labs/Aguiar/pHapCompass/simulated_data_NEW/contig_100/ploidy_3/merged_haplotypes.vcf --out /labs/Aguiar/pHapCompass/simulated_data_NEW/contig_100/ploidy_3/cov_10/frag/00.frag'


def extract_hairs():
    main_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW'
    extract_hairs_path = '/home/FCAM/mhosseini/HaplOrbit/extract_poly/build/extractHAIRS'
    sh_path = os.path.join('/labs/Aguiar/pHapCompass/scripts/simulation', '03_extract_hairs.sh')
    contig_lens = [100] # [10, 100, 300, 1000]
    ploidies = [3] # [3, 4, 6, 8, 10]
    coverages = [10] # [10, 20, 30, 40, 50]
    for contig_len in contig_lens:
        for ploidy in ploidies:
            this_sh_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '03_extract_hairs_{}_{}.sh'.format(contig_len, ploidy))
            to_print = 'module load bcftools/1.20\nmodule load htslib/1.20\n\n'
            this_vcf_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'maf0.01_hapref_chr21_filtered_NA12878_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
            this_vcf_sorted_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'sorted_maf0.01_hapref_chr21_filtered_NA12878_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
            unzipped_vcf = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'sorted_maf0.01_hapref_chr21_filtered_NA12878_ploidy{}_contig{}.vcf'.format(ploidy, contig_len))
            sort_cmd = 'bcftools sort {} -Oz -o {}\n'.format(this_vcf_path, this_vcf_sorted_path)
            idx_cmd = 'tabix -p vcf {}\n'.format(this_vcf_sorted_path)
            unzip_cmd = 'gunzip {}\n\n'.format(this_vcf_sorted_path)
            to_print += sort_cmd
            to_print += idx_cmd
            to_print += unzip_cmd
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                bam_path = os.path.join(this_cov_path, 'bam')
                frag_path = os.path.join(this_cov_path, 'frag')
                if not os.path.exists(frag_path):
                    os.makedirs(frag_path)
                for rd in range(100): 
                    command = '{} --bam {}/{}.bam --vcf {} --out {}/{}.frag\n'.format(extract_hairs_path, bam_path, str(rd).zfill(2), unzipped_vcf, frag_path, str(rd).zfill(2))
                    to_print += command
            with open(this_sh_path, 'w') as f:
                f.write(to_print)





if __name__ == '__main__':
    generate_short_reads()
