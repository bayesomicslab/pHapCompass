import os
import numpy as np
import pysam
import random
import pandas as pd 
from evaluation.evaluation import *
import pickle
import sys


def single_hpop_g_command_line(contig, ploidy, coverage, scaffold):
    cmd = f"""java -jar /mnt/research/aguiarlab/proj/HaplOrbit/comp_methods/H-PoPG/H-PoPG.jar -p {ploidy} \\
    -w 0.9 -f /mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_{contig}/ploidy_{ploidy}/cov_{coverage}/frag/{scaffold}.frag \\
    -vcf /mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_{contig}/ploidy_{ploidy}/HPOP_AWRI_ploidy{ploidy}_contig{contig}.vcf \\
    -o /mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_{contig}/ploidy_{ploidy}/cov_{coverage}/results/HPoP/{scaffold}.txt"""
    return cmd


def generate_hpop_command_runs():
    sh_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/comp_methods/hpop_runs/simulated_data_NA12878.sh'
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    hpop_path = '/mnt/research/aguiarlab/proj/HaplOrbit/comp_methods/H-PoPG/H-PoPG.jar'
    contig_lens = [100]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    to_print = ''
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # stop
            vcf_file = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'NA12878_ploidy{}_contig{}.vcf'.format(ploidy, contig_len))
            hpop_vcf_file = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'HPOP_NA12878_ploidy{}_contig{}.vcf'.format(ploidy, contig_len))
            convert_command = "\nsed 's/\\t\.\\t/\\tPASS\\t/g' {} > {}\n\n".format(vcf_file, hpop_vcf_file)
            to_print += convert_command
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                frag_path = os.path.join(this_cov_path, 'frag')
                results_path = os.path.join(this_cov_path, 'hpop_results')
                if not os.path.exists(results_path):
                    os.makedirs(results_path)
                for rd in range(n_samples):
                    frag_file = os.path.join(frag_path, str(rd).zfill(2) + '.frag')
                    output_file = os.path.join(results_path, str(rd).zfill(2) + '.txt')
                    command = 'java -jar {} -p {} -w 0.9 -f {} -vcf {} -o {}\n'.format(hpop_path, ploidy, frag_file, hpop_vcf_file, output_file)
                    to_print += command 
    with open(sh_path, 'w') as f:
        f.write(to_print)   


def parse_hpop_blocks_numpy(filename, true_haplotypes, fragment_list):
    """
    Parses a file containing blocks of haplotype information and stores data in NumPy matrices.

    Returns:
        blocks_dict (dict): 
            {
                block_index: {
                    'blocks': [pos1, pos2, ...],  # Sorted position array
                    'evaluation': {
                        'vector_error_rate': float,
                        'vector_error': int,
                        'accuracy': float,
                        'mismatch_error': float,
                        'mec': float
                    },
                    'haplotype_matrix': np.array([...])  # Haplotype matrix
                },
                ...
            }
    """
    predicted_haplotypes = pd.DataFrame(columns=true_haplotypes.columns, index=true_haplotypes.index)
    blocks_dict = {}
    block_index = 0

    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    while i < len(lines):
        line = lines[i]

        # Detect the start of a BLOCK
        if line.startswith("BLOCK: offset:"):
            parts = line.split()
            offset = int(parts[2])  # Starting position
            length = int(parts[4])  # Length of the block

            pos_list = []
            hap_list = []

            # Read the next `length` lines for this block
            start = i + 1
            end = i + 1 + length

            for j in range(start, end):
                row = lines[j].split()
                pos = int(row[0])  # Position index
                haps = [np.nan if x == '-' else int(x) for x in row[1:]]  # Replace '-' with np.nan
                
                pos_list.append(pos)
                hap_list.append(haps)
            
            # Convert lists to NumPy arrays
            block_pred_haplotype = np.array(hap_list, dtype=float).T  # Transpose to (haplotypes, positions)
            positions_ind = [p-1 for p in pos_list]
            predicted_haplotypes.loc[:, positions_ind] = block_pred_haplotype
            block_true_haplotype = true_haplotypes[positions_ind].to_numpy()

            # Select only columns that contain no NaNs
            valid_columns = list(np.where(~np.isnan(block_pred_haplotype).any(axis=0))[0])

            # Select only valid columns
            filtered_block_pred_haplotype = block_pred_haplotype[:, valid_columns] if valid_columns else np.array([])
            filtered_block_true_haplotype = block_true_haplotype[:, valid_columns] if valid_columns else np.array([])

            # Compute metrics conditionally
            if len(valid_columns) == 0:
                block_vector_error_rate = np.nan
                block_vector_error = np.nan
                block_accuracy = np.nan
                block_mismatch_error = np.nan
                block_mec_ = np.nan
            else:
                block_vector_error_rate, block_vector_error, _, _ = compute_vector_error_rate(filtered_block_pred_haplotype, filtered_block_true_haplotype)
                block_accuracy, _ = calculate_accuracy(filtered_block_pred_haplotype, filtered_block_true_haplotype)
                block_mismatch_error, _ = calculate_mismatch_error(filtered_block_pred_haplotype, filtered_block_true_haplotype) 
                block_mec_ = mec(filtered_block_pred_haplotype, fragment_list) 

            metrics = {
                'vector_error_rate': block_vector_error_rate,
                'vector_error': block_vector_error,
                'accuracy': block_accuracy,
                'mismatch_error': block_mismatch_error,
                'mec': block_mec_
            }

            # Store in dictionary
            blocks_dict[block_index] = {
                'blocks': pos_list,
                'evaluation': metrics, 
                'block_size': len(valid_columns)
            }

            block_index += 1
            i = end  # Move to the next block
        else:
            i += 1

    valid_predicted_columns = predicted_haplotypes.columns[~predicted_haplotypes.isna().any()].tolist()


    block_info = {'vector_error_rate': np.nanmean([blocks_dict[key]['evaluation']['vector_error_rate'] for key in blocks_dict.keys()]), 
                'vector_error': np.nanmean([blocks_dict[key]['evaluation']['vector_error'] for key in blocks_dict.keys()]),
                'accuracy': np.nanmean([blocks_dict[key]['evaluation']['accuracy'] for key in blocks_dict.keys()]),
                'mismatch_error': np.nanmean([blocks_dict[key]['evaluation']['mismatch_error'] for key in blocks_dict.keys()]),
                'mec': np.nanmean([blocks_dict[key]['evaluation']['mec'] for key in blocks_dict.keys()]), 
                'average_block_size': np.mean([blocks_dict[key]['block_size'] for key in blocks_dict.keys()]),
                'n_blocks': len(blocks_dict.keys()), 'length_phased': len(valid_predicted_columns), 'haplotypes': predicted_haplotypes}

    return blocks_dict, block_info


def hpop_to_haplotypes(result_file):
    # haplotypes are np.array
    haplotypes = []
    indices = []
    with open(result_file, 'r') as f:
        hpop_output = f.readlines()
    for line in hpop_output:
        line = line.strip('\n')
        positions = line.split("\t")
        
        # print(positions)
        # if positions[0] is not a number, then discard
        if positions[0].isdigit() == False:
            continue
        if "-" in positions:
            continue
        # if haplotypes is empty, then create lists for each piece within positions
        # print(f"Adding: {positions}")
        if len(haplotypes) == 0:
            for i in range(len(positions)-1):
                haplotypes.append([])
        indices.append(int(positions[0])-1)
        # add the positions to the haplotypes
        for i in range(len(positions[1:])):
            try:
                haplotypes[i].append(int(positions[i+1]))
            except:
                haplotypes[i].append("-")
    # convert haplotypes to np.array
    # print(haplotypes)
    haplotypes = np.array(haplotypes)
    # print("Haplotype array")
    # print(haplotypes)
    return haplotypes, indices


def get_fragment_list(frag_path):
    fragment_list = []
    for fragment in open(frag_path, 'r'):

        parts = fragment.split()
        positions = []
        alleles = []
        for iii in range(int(parts[0])):
            # process i+2,i+3.... i+4,i+5...
            start_idx_of_read = iii * 2 + 3
            seq_len = len(parts[start_idx_of_read])
            positions.extend(
                list(range(int(parts[start_idx_of_read - 1]), int(parts[start_idx_of_read - 1]) + seq_len)))
            [alleles.append(int(a)) for a in parts[start_idx_of_read]]
            fragment_list.append(positions)
            fragment_list.append(alleles)
    return fragment_list


def collect_results_hpop():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/hpop_results_simulated_data_NA12878_346.csv'
    contig_lens = [100]
    ploidies = [3, 4, 6]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value', 'length_phased'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*len(metrics)))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                results_path = os.path.join(this_cov_path, 'hpop_results')
                frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.txt')
                    cut_predicted_haplotypes, valid_columns = hpop_to_haplotypes(result_file)
                    cut_true_haplotypes = true_haplotypes[:, valid_columns]
                    fragment_list = get_fragment_list(os.path.join(frag_path, str(rd).zfill(2) + '.frag'))
                    hpmec = ev.mec(cut_predicted_haplotypes, fragment_list)
                    hpvector_error_rate, hpvector_error, backtracking_steps, dp_table = ev.compute_vector_error_rate(cut_predicted_haplotypes, cut_true_haplotypes)
                    hpmismatch_error, best_permutation = ev.calculate_mismatch_error(cut_predicted_haplotypes, cut_true_haplotypes)
                    hpaccuracy, _ = ev.calculate_accuracy(cut_predicted_haplotypes, cut_true_haplotypes)
                    phased_snp = cut_predicted_haplotypes.shape[1]
                    evals = {'vector_error_rate': hpvector_error_rate, 'vector_error': hpvector_error, 'accuracy': hpaccuracy, 'mismatch_error': hpmismatch_error, 'mec': hpmec}
                    for metric in metrics:
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = metric
                        result_df.loc[counter, 'Value'] = evals[metric]
                        result_df.loc[counter, 'length_phased'] = phased_snp
                        counter += 1
    result_df['Method'] = 'H-PoPG'
    result_df.to_csv(output_path, index=False)
    print("H-PoPG results collected")


def collect_results_hpop_missing():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/hpop_results_simulated_data_NA12878_346.csv'
    contig_lens = [100]
    ploidies = [3, 4, 6]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    # metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*2))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                results_path = os.path.join(this_cov_path, 'hpop_results')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.txt')
                    cut_predicted_haplotypes, valid_columns = hpop_to_haplotypes(result_file)
                    predicted_haplotypes = np.full((ploidy, contig_len), np.nan)
                    predicted_haplotypes[:, valid_columns] = cut_predicted_haplotypes
                    length_phased = cut_predicted_haplotypes.shape[1]
                    vector_error_rate, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes, predicted_haplotypes)
                    result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                    result_df.loc[counter, 'Metric'] = 'Vector Error Rate'
                    result_df.loc[counter, 'Value'] = vector_error_rate
                    result_df.loc[counter, 'Contig'] = contig_len
                    result_df.loc[counter, 'Ploidy'] = ploidy
                    result_df.loc[counter, 'Coverage'] = coverage
                    counter += 1
                    result_df.loc[counter, 'Contig'] = contig_len
                    result_df.loc[counter, 'Ploidy'] = ploidy
                    result_df.loc[counter, 'Coverage'] = coverage
                    result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                    result_df.loc[counter, 'Metric'] = '# Phased Variants'
                    result_df.loc[counter, 'Value'] = length_phased
                    counter += 1

    result_df['Method'] = 'H-PoPG'
    result_df.to_csv(output_path, index=False)
    print("H-PoPG results collected")


def collect_results_hpop_blocks():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/hpop_results_simulated_data_NA12878_346_block.csv'
    contig_lens = [100]
    ploidies = [3, 4, 6]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec', 'length_phased', 'n_blocks', 'average_block_size']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*len(metrics)))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                results_path = os.path.join(this_cov_path, 'hpop_results')
                frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.txt')
                    fragment_list = get_fragment_list(os.path.join(frag_path, str(rd).zfill(2) + '.frag'))
                    _, block_info = parse_hpop_blocks_numpy(result_file, true_haplotypes, fragment_list)

                    for metric in metrics:
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = metric
                        result_df.loc[counter, 'Value'] = block_info[metric]
                        counter += 1
    result_df['Method'] = 'H-PoPG'
    result_df.to_csv(output_path, index=False)
    print("H-PoPG results collected")


def collect_results_hpop_from_pkl():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/hpop_results_simulated_data_NA12878_8.csv'
    contig_lens = [100]
    ploidies = [8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    # metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*2))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                results_path = os.path.join(this_cov_path, 'hpop_results')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.pkl')
                    with open(result_file, "rb") as f:
                        evals = pickle.load(f)
                    result_df.loc[counter, 'Sample'] = evals['Sample']
                    result_df.loc[counter, 'Metric'] = 'Vector Error Rate'
                    result_df.loc[counter, 'Value'] = evals['vector_error_rate']
                    result_df.loc[counter, 'Contig'] = contig_len
                    result_df.loc[counter, 'Ploidy'] = ploidy
                    result_df.loc[counter, 'Coverage'] = coverage
                    counter += 1
                    result_df.loc[counter, 'Contig'] = contig_len
                    result_df.loc[counter, 'Ploidy'] = ploidy
                    result_df.loc[counter, 'Coverage'] = coverage
                    result_df.loc[counter, 'Sample'] = evals['Sample']
                    result_df.loc[counter, 'Metric'] = '# Phased Variants'
                    result_df.loc[counter, 'Value'] = evals['length_phased']
                    counter += 1

    result_df['Method'] = 'H-PoPG'
    result_df.to_csv(output_path, index=False)
    print("H-PoPG results collected")


def collect_results_hpop_from_pkl_block():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/hpop_results_simulated_data_NA12878_8_block.csv'
    contig_lens = [100]
    ploidies = [8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec', 'length_phased', 'n_blocks', 'average_block_size']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*len(metrics)))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                results_path = os.path.join(this_cov_path, 'hpop_results')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.pkl')
                    with open(result_file, "rb") as f:
                        evals = pickle.load(f)
                    for metric in metrics:
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = metric
                        result_df.loc[counter, 'Value'] = evals[metric]
                        counter += 1
    result_df['Method'] = 'H-PoPG'
    result_df.to_csv(output_path, index=False)
    print("H-PoPG results collected")


def make_inputs_for_evals_hpop():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_dir = '/mnt/research/aguiarlab/proj/HaplOrbit/hpopg_inputs'
    # output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/hpop_results_simulated_data_NA12878.csv'
    contig_lens = [100]
    ploidies = [8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    inputs = []
    # metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    # result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value', 'length_phased'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*len(metrics)))
    # counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            true_haplotypes_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                results_path = os.path.join(this_cov_path, 'hpop_results')
                frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.txt')
                    inputs.append([result_file, true_haplotypes_path, frag_path, contig_len, ploidy, coverage, rd])

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")


def eval_one_input_hpop_missing(inp):
    result_file, true_haplotypes_path, frag_path, contig_len, ploidy, coverage, rd = inp
    true_haplotypes = pd.read_csv(true_haplotypes_path).T.to_numpy()
    cut_predicted_haplotypes, valid_columns = hpop_to_haplotypes(result_file)

    predicted_haplotypes = np.full((ploidy, contig_len), np.nan)
    predicted_haplotypes[:, valid_columns] = cut_predicted_haplotypes
    length_phased = cut_predicted_haplotypes.shape[1]
    vector_error_rate, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes, predicted_haplotypes)

    pkl_name = result_file.replace('.txt', '.pkl')

    evals = {'Contig': contig_len, 'Ploidy': ploidy, 'Coverage': coverage , 'Sample': str(rd).zfill(2), 'vector_error_rate': vector_error_rate, 'length_phased': length_phased}    
    with open(pkl_name, "wb") as f:
        pickle.dump(evals, f)
    print(f"Done with {pkl_name}")


def eval_one_input_hpop_block(inp):
    result_file, true_haplotypes_path, frag_path, contig_len, ploidy, coverage, rd = inp
    true_haplotypes = pd.read_csv(true_haplotypes_path).T
    fragment_list = get_fragment_list(os.path.join(frag_path, str(rd).zfill(2) + '.frag'))

    pkl_name = result_file.replace('.txt', '.pkl')
    _, block_info = parse_hpop_blocks_numpy(result_file, true_haplotypes, fragment_list)
    
    evals = {'Contig': contig_len, 'Ploidy': ploidy, 'Coverage': coverage , 'Sample': str(rd).zfill(2), 'vector_error_rate': block_info['vector_error_rate'], 
             'vector_error': block_info['vector_error'], 'accuracy': block_info['accuracy'], 'mismatch_error': block_info['mismatch_error'], 'mec': block_info['mec'], 
             'length_phased': block_info['length_phased'], 'n_blocks': block_info['n_blocks'], 'average_block_size': block_info['average_block_size']}
    
    with open(pkl_name, "wb") as f:
        pickle.dump(evals, f)
    print(f"Done with {pkl_name}")


def eval_one_input_from_input_hpop(input_file):
    with open(input_file, "rb") as f:
        inp = pickle.load(f)

    # eval_one_input_hpop_block(inp)
    eval_one_input_hpop_missing(inp)


if __name__ == '__main__':
    # generate_hpop_runs()
    # collect_results_hpop()


    if len(sys.argv) != 2:
        print("Usage: python3 compare_hpop.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    eval_one_input_from_input_hpop(input_file)