import os
import pickle
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri'
metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
contigs = ['10', '100'] # 100, 1000
ploidies = ['3', '4', '6', '8']
coverages = ['10', '20', '30', '40', '50']
results_dfs = []
for contig in contigs:
    for ploidy in ploidies:
        for coverage in coverages:
            # results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results')
            results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results_likelihood')
            if os.path.exists(results_path):
                samples = [f for f in os.listdir(results_path) if 'FFBS' in f]
                this_result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(samples)*len(metrics)))
                this_result_df['Contig'] = contig
                this_result_df['Ploidy'] = ploidy
                this_result_df['Coverage'] = coverage
                this_result_df['Method'] = 'pHapCompass'                
                counter = 0
                for sample in samples:
                    sample_result = os.path.join(results_path, sample)
                    with open(sample_result, 'rb') as f:
                        this_results = pickle.load(f)
                    sample_name = sample.split('.pkl')[0].split('_')[-1]
                    evals = this_results['evaluation']
                    for metric in metrics:
                        this_result_df.loc[counter, 'Sample'] = sample_name
                        this_result_df.loc[counter, 'Metric'] = metric
                        this_result_df.loc[counter, 'Value'] = evals[metric]
                        counter += 1
                results_dfs.append(this_result_df)

results_df = pd.concat(results_dfs, ignore_index=True)
results_df.to_csv(os.path.join(agg_results_path, 'sim_awri_results_likelihood100.csv'), index=False)


results100 = pd.read_csv(os.path.join(agg_results_path, 'sim_awri_results_likelihood100.csv'))
results10 = pd.read_csv(os.path.join(agg_results_path, 'sim_awri_results_likelihood10.csv'))

all_results = pd.concat([results100, results10], ignore_index=True)

# results_df_100 = results_df[results_df['contig'] == '100'].reset_index(drop=True)
results_df_10 = all_results[all_results['Contig'] == 10].reset_index(drop=True)
results_df_100 = all_results[all_results['Contig'] == 100].reset_index(drop=True)

for metric in metrics:
    print('Metric:', metric)
    metric_df = results_df_10[results_df_10['Metric'] == metric].reset_index(drop=True)
    sns.catplot(x="Coverage", y="Value",hue="Ploidy", data=metric_df, kind="violin", height=6, aspect=1.5)
    plt.title(f"{metric.capitalize()}")
    plt.xlabel("Coverage")
    plt.ylabel("Value")
    plt.tight_layout()
    plt.savefig(os.path.join(agg_results_path, f"{metric}_10.png"))
    plt.close()
