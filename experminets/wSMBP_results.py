import pickle
import os
import pandas as pd
import itertools

def produce_results(results_path, args):
    results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    eval_df = pd.DataFrame(columns=['Contig', 'Ploidy', 'Coverage', 'Sample', 'Weighting', 'Vector Error Rate', 'Accuracy', 'MEC', 'Mismatch Error', 'Pairwise Phasing Accuracy'])
    counter = 0
    for contig, ploidy, coverage, sample, weighting in args:
        try:
            print(f"Processing contig {contig}, ploidy {ploidy}, coverage {coverage}, sample {sample}, weighting {weighting}")
            sample_result = os.path.join(main_path, 'contig_{}/ploidy_{}/cov_{}/results_wSMBP/weighting_{}/wSMBP_{}.pkl'.format(contig, ploidy, coverage, str(weighting), sample))
            with open(sample_result, 'rb') as f:
                this_results = pickle.load(f)
                f.close()
            eval_df.loc[counter, 'Contig'] = contig
            eval_df.loc[counter, 'Ploidy'] = ploidy
            eval_df.loc[counter, 'Coverage'] = coverage
            eval_df.loc[counter, 'Sample'] = sample
            eval_df.loc[counter, 'Weighting'] = str(weighting)
            eval_df.loc[counter, 'Vector Error Rate'] = this_results['evaluation']['vector_error_rate']
            eval_df.loc[counter, 'Accuracy'] = this_results['evaluation']['accuracy']
            eval_df.loc[counter, 'MEC'] = this_results['evaluation']['mec']
            eval_df.loc[counter, 'Mismatch Error'] = this_results['evaluation']['mismatch_error']
            eval_df.loc[counter, 'Pairwise Phasing Accuracy'] = this_results['evaluation']['ffbs_acc']
            counter += 1
        except FileNotFoundError:
            print(f"File not found for contig {contig}, ploidy {ploidy}, coverage {coverage}, weighting {weighting}, sample {sample}. Skipping.")
            continue
                    
    # eval_groups = eval_df.groupby(['Contig', 'Ploidy', 'Coverage'])['FFBS_Acc'].mean().reset_index()
    eval_df.to_csv(os.path.join(results_path, 'wSMBP_full_results_v2_NA12878.csv'), index=False)

if __name__ == '__main__':
    results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    contigs = [100]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70]
    weightings = [None, 'standard', 'mean', 'max']
    samples = [str(i).zfill(2) for i in range(100)]

    args = list(itertools.product(contigs, ploidies, coverages, samples, weightings))
    produce_results(results_path, args)