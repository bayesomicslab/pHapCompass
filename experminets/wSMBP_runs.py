import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from algorithm.wSMBP import pipeline
import pickle
import pandas as pd
from multiprocessing import Pool

def evaluate_ffbs_acc(results_path, contigs=[100, 1000], ploidies=[3, 4, 6, 8], coverages=[10, 30, 50, 70, 100], samples=None, weightings=[None, 'standard', 'mean', 'max']):
    results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    eval_df = pd.DataFrame(columns=['Contig', 'Ploidy', 'Coverage', 'Sample', 'Weighting', 'FFBS_Acc'])
    counter = 0 
    for contig in contigs:
        for ploidy in ploidies:
            for cov in coverages:
                for weighting in weightings:
                    for sample_name in samples:
                        try:
                            # print('Contig:', contig, 'Ploidy:', ploidy, 'Coverage:', cov)
                            # ploidy = 4
                            sample_result = os.path.join(main_path, 'contig_{}/ploidy_{}/cov_{}/results_wSMBP/weighting_{}/wSMBP_{}.pkl'.format(contig, ploidy, cov, str(weighting), sample_name)) 
                            with open(sample_result, 'rb') as f:
                                this_results = pickle.load(f)
                                f.close()
                            ffbs_acc = this_results['evaluation']['ffbs_acc']
                            eval_df.loc[counter, 'Contig'] = contig
                            eval_df.loc[counter, 'Ploidy'] = ploidy
                            eval_df.loc[counter, 'Coverage'] = cov
                            eval_df.loc[counter, 'Sample'] = sample_name
                            eval_df.loc[counter, 'Weighting'] = weighting
                            eval_df.loc[counter, 'FFBS_Acc'] = ffbs_acc
                            counter += 1
                        except FileNotFoundError:
                            print(f"File not found for contig {contig}, ploidy {ploidy}, coverage {cov}, divide {divide}, sample {sample_name}. Skipping.")
                            continue
    # eval_groups = eval_df.groupby(['Contig', 'Ploidy', 'Coverage'])['FFBS_Acc'].mean().reset_index()
    eval_df.to_csv(os.path.join(results_path, 'wSMBP_sample_acc_NA12878.csv'), index=False)

def run_pipeline(args):
    contig, ploidy, coverage, sample, data, weighting = args
    try:
        pipeline(ploidy=ploidy, coverage=coverage, contig=contig, sample=sample, data=data, weighting=weighting)
    except Exception as e:
        print(f"Error running pipeline for contig {contig}, ploidy {ploidy}, coverage {coverage}, sample {sample}, weighting {weighting}: {e}")

def run_pipelines(data, contigs, ploidies, coverages, samples, weightings, num_processes=15):
    args_list = [(contig, ploidy, coverage, sample, data, weighting) for contig in contigs for ploidy in ploidies for coverage in coverages for sample in samples for weighting in weightings]
    processes_pool = Pool(num_processes)
    processes_pool.map(run_pipeline, args_list)
    processes_pool.close()
    processes_pool.join()
    print("All pipelines completed.")

if __name__ == '__main__':
    data = 'NA12878'
    contigs = [100]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70]
    weightings = [None, 'standard', 'mean', 'max']
    samples = [str(i).zfill(2) for i in range(100)]

    run_pipelines(data, contigs, ploidies, coverages, samples, weightings, num_processes=15)
