import os
import numpy as np
import pysam
import random
import pandas as pd 
from evaluation import evaluation as ev

def generate_whatshapp_runs():
    sh_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/comp_methods/whatshapp_runs/simulated_data_test.sh'
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test'
    contig_lens = [100]
    ploidies = [6]
    coverages = [10, 50, 100]
    n_samples = 10
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

# vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test/contig_100/ploidy_6/cov_10/whatshapp_results/00.vcf.gz'
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


def collect_results_whatshap():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/whatshap_results_simulated_data_test.csv'
    contig_lens = [100]
    ploidies = [6]
    coverages = [10, 50, 100]
    n_samples = 10
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value', 'length_phased'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*len(metrics)))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                results_path = os.path.join(this_cov_path, 'whatshapp_results')
                frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    result_file = os.path.join(results_path, str(rd).zfill(2) + '.vcf.gz')
                    predicted_haplotypes, indices = vcf_to_haplotypes(result_file)
                    valid_columns = [c for c in range(predicted_haplotypes.shape[1]) if not np.any(np.isnan(predicted_haplotypes[:, c]))]
                    cut_predicted_haplotypes = predicted_haplotypes[:, valid_columns]
                    cut_true_haplotypes = true_haplotypes[:, valid_columns]
                    fragment_list = get_fragment_list(os.path.join(frag_path, str(rd).zfill(2) + '.frag'))
                    wpmec = ev.mec(cut_predicted_haplotypes, fragment_list)
                    wpvector_error_rate, wpvector_error, backtracking_steps, dp_table = ev.compute_vector_error_rate(cut_predicted_haplotypes, cut_true_haplotypes)
                    wpmismatch_error, best_permutation = ev.calculate_mismatch_error(cut_predicted_haplotypes, cut_true_haplotypes)
                    wpaccuracy, _ = ev.calculate_accuracy(cut_predicted_haplotypes, cut_true_haplotypes)
                    phased_snp = cut_predicted_haplotypes.shape[1]
                    evals = {'vector_error_rate': wpvector_error_rate, 'vector_error': wpvector_error, 'accuracy': wpaccuracy, 'mismatch_error': wpmismatch_error, 'mec': wpmec}
                    for metric in metrics:
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = metric
                        result_df.loc[counter, 'Value'] = evals[metric]
                        result_df.loc[counter, 'length_phased'] = phased_snp
                        counter += 1
    result_df['Method'] = 'WhatsHap'
    result_df.to_csv(output_path, index=False)


