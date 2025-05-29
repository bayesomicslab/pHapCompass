import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from algorithm.SMBP import pipeline
import pickle
import pandas as pd

def evaluate_ffbs_acc(results_path, ploidy=3, cov=10, contig=10):
    results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri'
    eval_df = pd.DataFrame(columns=['Contig', 'Ploidy', 'Coverage', 'Sample', 'FFBS_Acc'])
    counter = 0 
    for contig in [contig]:
        for ploidy in [ploidy]:
            for cov in [cov]:
                if True:
                    sample_name = '00'
                    print('Contig:', contig, 'Ploidy:', ploidy, 'Coverage:', cov)
                    # ploidy = 4
                    sample_result = os.path.join(main_path, 'contig_{}/ploidy_{}/cov_{}/results_SMBP/SMBP_{}.pkl'.format(contig, ploidy, cov, sample_name)) 
                    with open(sample_result, 'rb') as f:
                        this_results = pickle.load(f)
                        f.close()
                    ffbs_acc = this_results['evaluation']['ffbs_acc']
                    eval_df.loc[counter, 'Contig'] = contig
                    eval_df.loc[counter, 'Ploidy'] = ploidy
                    eval_df.loc[counter, 'Coverage'] = cov
                    eval_df.loc[counter, 'Sample'] = sample_name
                    eval_df.loc[counter, 'FFBS_Acc'] = ffbs_acc
                    print("Predicted Haplotypes")
                    print(this_results['predicted_haplotypes'])
                    print("True Haplotypes")
                    print(this_results['true_haplotypes'])
                    counter += 1
    # eval_groups = eval_df.groupby(['Contig', 'Ploidy', 'Coverage'])['FFBS_Acc'].mean().reset_index()
    print(eval_df)

if __name__ == '__main__':
    data = 'awri'
    contigs = [100, 1000]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70, 100]
    # samples = 00 - 99
    samples = [str(i).zfill(2) for i in range(100)]

    contigs = [10]
    ploidies = [3]
    coverages = [10]
    samples = ['00']

    for contig in contigs:
        for ploidy in ploidies:
            for coverage in coverages:
                for sample in samples:
                    pipeline(ploidy=ploidy, coverage=coverage, contig=contig, sample=sample, data=data)
    
    evaluate_ffbs_acc('/mnt/research/aguiarlab/proj/HaplOrbit/results', ploidy=3, cov=10, contig=10)