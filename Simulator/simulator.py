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



def simulate_from_diploid():
    vcf_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/_EGAZ00001239269_HRC.r1-1.EGA.GRCh37.chr2.haplotypes.vcf.gz'
    fasta_path = '/mnt/research/aguiarlab/data/hg19/chr2.fa'
    
    with open(fasta_path) as f:
        fasta_file = f.readlines()




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


if __name__ == '__main__':
    generate_short_reads()
