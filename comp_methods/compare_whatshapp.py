import os
import numpy as np
import pysam
import random
import pandas as pd 
from evaluation.evaluation import *
import pickle
import sys


def generate_whatshapp_runs():
    sh_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/comp_methods/whatshapp_runs/simulated_data_NA12878.sh'
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    contig_lens = [100]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    to_print = ''
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # stop
            vcf_file = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'NA12878_ploidy{}_contig{}.vcf'.format(ploidy, contig_len))
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                bam_path = os.path.join(this_cov_path, 'bam')
                results_path = os.path.join(this_cov_path, 'whatshapp_results')
                if not os.path.exists(results_path):
                    os.makedirs(results_path)
                for rd in range(n_samples):
                    bam_file = os.path.join(bam_path, str(rd).zfill(2) + '.bam')
                    output_file = os.path.join(results_path, str(rd).zfill(2) + '.vcf.gz')
                    command = 'whatshap polyphase --ignore-read-groups {} {} --ploidy {} --output {}\n'.format(vcf_file, bam_file, ploidy, output_file)
                    to_print += command 
    with open(sh_path, 'w') as f:
        f.write(to_print)   


def vcf_to_haplotypes(vcf_path):
    haplotypes = []
    indices = []
    vcf = pysam.VariantFile(vcf_path)
    i = 0
    for rec in vcf.fetch():
        
        for sample in rec.samples:
            
            # print(sample)
            # Access the genotype field for the sample (GT)
            is_phased = rec.samples[sample].phased
            
            gt = rec.samples[sample]["GT"]

            # if there's a non-digit in gt tuple, skip
            if all(isinstance(x, int) for x in gt):
                if not is_phased:
                    gt = [np.nan for x in gt]
                    
                haplotypes.append(list(gt))
                indices.append(i)
            # print(gt)
            i += 1
    haplotypes = np.array(haplotypes)
    haplotypes = haplotypes.transpose()
    # print(haplotypes)
    return haplotypes, indices


def parse_vcf_consecutive_blocks(vcf_path, true_haplotypes, fragment_list):
    """
    Reads a single-sample VCF, identifies blocks of consecutive, phased records
    that share the same PS tag, and computes metrics vs. a 'true_haplotypes' DataFrame.

    Each "block" is formed by lines i, i+1, i+2, ... in the VCF that:
      (1) Are phased (phased==True),
      (2) Have the same PS tag,
      (3) Appear consecutively until a line breaks the condition 
          (different PS or unphased).
    Unphased lines are ignored (not part of any block).

    Args:
        vcf_path (str):
            Path to the input VCF file, single-sample.
        true_haplotypes (pd.DataFrame):
            shape = (ploidy, n_variants),
            columns = [0, 1, 2, ..., n_variants-1],
            index = haplotype IDs. This is the "true" genotype/haplotype data.
        fragment_list (object):
            Additional data needed by some metric function (ev.mec).

    Returns:
        blocks_dict (dict):
            {
                block_index: {
                    'blocks': [i, i+1, ...],  # global consecutive indices in the VCF
                    'evaluation': {
                        'vector_error_rate': float,
                        'vector_error': float,
                        'accuracy': float,
                        'mismatch_error': float,
                        'mec': float
                    },
                    'block_size': int
                },
                ...
            }

        block_info (dict):
            {
                'vector_error_rate': float,
                'vector_error': float,
                'accuracy': float,
                'mismatch_error': float,
                'mec': float,
                'average_block_size': float,
                'n_blocks': int,
                'length_phased': int,
                'haplotypes': predicted_haplotypes (DataFrame)
            }
    """

    # Make sure predicted_haplotypes has the same shape as true_haplotypes
    # (ploidy x number_of_variants), initially filled with NaN
    n_variants = true_haplotypes.shape[1]  # total columns
    predicted_haplotypes = pd.DataFrame(
        data=np.nan, 
        index=true_haplotypes.index, 
        columns=true_haplotypes.columns
    )

    vcf = pysam.VariantFile(vcf_path)
    sample = None

    # We'll store the relevant info for each phased variant in a list.
    # We skip unphased variants entirely.
    variant_records = []  # each item: dict(pos=..., global_idx=..., ps=..., gt=[...])

    # Global counter i for "consecutive index" across *all* lines (phased or not).
    # We only increment for each line in the VCF that we *examine*, whether phased or unphased.
    i = 0  

    for rec in vcf.fetch():
        # We assume single-sample VCF, so get that sample name (only once)
        if sample is None:
            sample = next(iter(rec.samples))

        # Check if genotype is phased
        is_phased = rec.samples[sample].phased
        gt = rec.samples[sample].get("GT", None)  # e.g. (1, 0, 1)
        ps = rec.samples[sample].get("PS", None)  # e.g. 9880113

        if gt is not None and is_phased and ps is not None:
            # Additional check: must all be integers
            if all(isinstance(allele, int) for allele in gt):
                # Keep this variant in our record
                variant_records.append({
                    "global_idx": i,     # i => the consecutive index
                    "ps": ps,            # phase set tag
                    "gt": list(gt)       # genotype array
                })
        # If unphased or missing PS, we do NOT store it in variant_records

        i += 1  # increment for each line processed (phased or not)

    # Now we have a list of "phased" records only, each with a "global_idx" and "ps".
    # Next, we identify consecutive blocks with the same PS.

    blocks = []
    current_ps = None
    current_block = []

    # We'll keep track of the last global_idx we handled, so we know if consecutive in VCF lines
    last_idx = None

    for record in variant_records:
        ps_val = record["ps"]
        idx_val = record["global_idx"]

        # If block is empty => start a new one
        if not current_block:
            current_ps = ps_val
            current_block = [record]
            last_idx = idx_val
            continue

        # If the same PS and "consecutive in terms of VCF lines" => belongs to current block
        # i.e. idx_val == last_idx + 1
        if ps_val == current_ps and idx_val == (last_idx + 1):
            current_block.append(record)
            last_idx = idx_val
        else:
            # finalize old block
            blocks.append(current_block)
            # start a new block
            current_block = [record]
            current_ps = ps_val
            last_idx = idx_val

    # finalize last block if not empty
    if current_block:
        blocks.append(current_block)

    # Each element of blocks is a list of consecutive variants with the same PS.
    # Now let's fill predicted_haplotypes columns, compute metrics, etc.

    blocks_dict = {}
    block_index = 0

    for block_records in blocks:
        # block_records is a list of dicts with keys: "ps", "global_idx", "gt"
        global_indices = [r["global_idx"] for r in block_records]
        # store these consecutive indices as "blocks": [...]
        # We also need to fill predicted_haplotypes for these columns.

        # Build the genotype matrix for this block => shape: (ploidy, number_of_positions_in_block)
        # We'll gather them in order, i.e. sorted by global_idx
        block_records_sorted = sorted(block_records, key=lambda x: x["global_idx"])
        block_gts = [r["gt"] for r in block_records_sorted]  # shape (num_positions, ploidy)
        block_gts = np.array(block_gts, dtype=float).T  # => shape (ploidy, num_positions)

        # Fill predicted haplotypes
        # Because we sorted by global_idx, let's extract them again in sorted order
        sorted_indices = [r["global_idx"] for r in block_records_sorted]
        predicted_haplotypes.loc[:, sorted_indices] = block_gts

        # Now get the corresponding columns from true_haplotypes
        # (assuming these columns exist: 0..(n_variants-1)).
        block_true = true_haplotypes.loc[:, sorted_indices].to_numpy()

        # Filter out columns from block_gts that contain any NaN
        valid_cols = np.where(~np.isnan(block_gts).any(axis=0))[0]

        if len(valid_cols) == 0:
            # no valid columns => all metrics are NaN
            block_vector_error_rate = np.nan
            block_vector_error = np.nan
            block_accuracy = np.nan
            block_mismatch_error = np.nan
            block_mec = np.nan
        else:
            block_pred_filt = block_gts[:, valid_cols]
            block_true_filt = block_true[:, valid_cols]

            # Now compute your metrics:
            block_vector_error_rate, block_vector_error, _, _ = ev.compute_vector_error_rate(block_pred_filt, block_true_filt)
            block_accuracy, _ = ev.calculate_accuracy(block_pred_filt, block_true_filt)
            block_mismatch_error, _ = ev.calculate_mismatch_error(block_pred_filt, block_true_filt)
            block_mec = ev.mec(block_pred_filt, fragment_list)

        metrics = {
            'vector_error_rate': block_vector_error_rate,
            'vector_error': block_vector_error,
            'accuracy': block_accuracy,
            'mismatch_error': block_mismatch_error,
            'mec': block_mec
        }

        blocks_dict[block_index] = {
            # "blocks" => the global consecutive indices, not the real genomic positions
            'blocks': global_indices,
            'evaluation': metrics,
            'block_size': len(valid_cols)
        }

        block_index += 1

    # After processing all blocks, compute aggregated block_info
    # We also can figure out how many columns have no NaNs in predicted_haplotypes
    valid_pred_columns = predicted_haplotypes.columns[~predicted_haplotypes.isna().any()].tolist()

    block_info = {
        'vector_error_rate': np.nanmean([blocks_dict[k]['evaluation']['vector_error_rate'] for k in blocks_dict]),
        'vector_error': np.nanmean([blocks_dict[k]['evaluation']['vector_error'] for k in blocks_dict]),
        'accuracy': np.nanmean([blocks_dict[k]['evaluation']['accuracy'] for k in blocks_dict]),
        'mismatch_error': np.nanmean([blocks_dict[k]['evaluation']['mismatch_error'] for k in blocks_dict]),
        'mec': np.nanmean([blocks_dict[k]['evaluation']['mec'] for k in blocks_dict]),
        'average_block_size': np.mean([blocks_dict[k]['block_size'] for k in blocks_dict]),
        'n_blocks': len(blocks_dict),
        'length_phased': len(valid_pred_columns),
        'haplotypes': predicted_haplotypes
    }

    return blocks_dict, block_info


def randomize_orders():
    vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test/contig_100/ploidy_6/NA12878_ploidy6_contig100.vcf'
    output_vcf_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test/contig_100/ploidy_6/randomized_NA12878_ploidy6_contig100.vcf"

    # Open the input VCF file
    vcf = pysam.VariantFile(vcf_path)

    # Create an output VCF file with the same header
    output_vcf = pysam.VariantFile(output_vcf_path, "w", header=vcf.header)

    for rec in vcf.fetch():
        for sample in rec.samples:
            # Access the genotype field for the sample (GT)
            gt = list(rec.samples[sample]["GT"])
            
            # Ensure GT does not contain None values
            if None in gt:
                continue  # Skip records with missing genotypes
            
            # Shuffle the genotype values while preserving the count of 0s and 1s
            random.shuffle(gt)
            
            # Assign the new randomized genotype
            rec.samples[sample]["GT"] = tuple(gt)
        
        # Write the modified record to the output VCF
        output_vcf.write(rec)

    # Close the files
    vcf.close()
    output_vcf.close()


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


def collect_results_whatshap_missing():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/whatshap_results_simulated_data_NA12878_346.csv'
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
                results_path = os.path.join(this_cov_path, 'whatshapp_results')
                frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.vcf.gz')
                    predicted_haplotypes, _ = vcf_to_haplotypes(result_file)
                    
                    valid_columns = [c for c in range(predicted_haplotypes.shape[1]) if not np.any(np.isnan(predicted_haplotypes[:, c]))]

                    length_phased = len(valid_columns)
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

    result_df['Method'] = 'WhatsHap'
    result_df.to_csv(output_path, index=False)
    print("WhatsHap results collected")



def collect_results_whatshap_block():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/whatshap_results_simulated_data_NA12878_346_block.csv'
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
                results_path = os.path.join(this_cov_path, 'whatshapp_results')
                frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.vcf.gz')
                    fragment_list = get_fragment_list(os.path.join(frag_path, str(rd).zfill(2) + '.frag'))
                    _, block_info = parse_vcf_consecutive_blocks(result_file, true_haplotypes, fragment_list)
                    for metric in metrics:
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = metric
                        result_df.loc[counter, 'Value'] = block_info[metric]
                        counter += 1

    result_df['Method'] = 'WhatsHap'
    result_df.to_csv(output_path, index=False)


def collect_results_whatshap_from_pkl():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/whatshap_results_simulated_data_NA12878_ploidy8.csv'
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
                results_path = os.path.join(this_cov_path, 'whatshapp_results')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.pkl')

                    with open(result_file, 'rb') as f:
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

    result_df['Method'] = 'WhatsHap'
    result_df.to_csv(output_path, index=False)


def collect_results_whatshap_from_pkl_block():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/whatshap_results_simulated_data_NA12878_ploidy8_block.csv'
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
                results_path = os.path.join(this_cov_path, 'whatshapp_results')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.pkl')
                    with open(result_file, 'rb') as f:
                        evals = pickle.load(f)
                    for metric in metrics:
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = metric
                        result_df.loc[counter, 'Value'] = evals[metric]
                        counter += 1
    result_df['Method'] = 'WhatsHap'
    result_df.to_csv(output_path, index=False)



def make_inputs_for_evals():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_dir = '/mnt/research/aguiarlab/proj/HaplOrbit/whatshap_inputs'
    # output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/whatshap_results_simulated_data_NA12878.csv'
    inputs = []
    contig_lens = [100]
    ploidies = [8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            true_haplotypes_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                results_path = os.path.join(this_cov_path, 'whatshapp_results')
                frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.vcf.gz')
                    inputs.append([result_file, true_haplotypes_path, frag_path, rd, contig_len, ploidy, coverage])
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")
    # return inputs


def eval_one_input_missing(inp):
    result_file, true_haplotypes_path, frag_path, rd, contig_len, ploidy, coverage = inp
    true_haplotypes = pd.read_csv(true_haplotypes_path).T.to_numpy()
    predicted_haplotypes, _ = vcf_to_haplotypes(result_file)
    valid_columns = [c for c in range(predicted_haplotypes.shape[1]) if not np.any(np.isnan(predicted_haplotypes[:, c]))]
    
    length_phased = len(valid_columns)
    vector_error_rate, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes, predicted_haplotypes)
    pkl_name = result_file.replace('.vcf.gz', '.pkl')
    
    evals = {'Contig': contig_len, 'Ploidy': ploidy, 'Coverage': coverage , 'Sample': str(rd).zfill(2),
             'vector_error_rate': vector_error_rate, 'length_phased': length_phased}
    
    with open(pkl_name, 'wb') as f:
        pickle.dump(evals, f)
    print(f"Done with {pkl_name}")


def eval_one_input_block(inp):
    result_file, true_haplotypes_path, frag_path, rd, contig_len, ploidy, coverage = inp
    true_haplotypes = pd.read_csv(true_haplotypes_path).T
    fragment_list = get_fragment_list(os.path.join(frag_path, str(rd).zfill(2) + '.frag'))

    pkl_name = result_file.replace('.vcf.gz', '.pkl')

    _, block_info = parse_vcf_consecutive_blocks(result_file, true_haplotypes, fragment_list)
    
    evals = {'Contig': contig_len, 'Ploidy': ploidy, 'Coverage': coverage , 'Sample': str(rd).zfill(2), 'vector_error_rate': block_info['vector_error_rate'], 
             'vector_error': block_info['vector_error'], 'accuracy': block_info['accuracy'], 'mismatch_error': block_info['mismatch_error'], 'mec': block_info['mec'], 
             'length_phased': block_info['length_phased'], 'n_blocks': block_info['n_blocks'], 'average_block_size': block_info['average_block_size']}
    
    with open(pkl_name, "wb") as f:
        pickle.dump(evals, f)
    print(f"Done with {pkl_name}")


def eval_one_input_from_input(input_file):
    with open(input_file, "rb") as f:
        inp = pickle.load(f)
    eval_one_input_missing(inp)


if __name__ == '__main__':
    # collect_results_whatshap_missing()


    if len(sys.argv) != 2:
        print("Usage: python3 compare_whatshapp.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    eval_one_input_from_input(input_file)