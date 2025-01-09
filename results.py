import os
import pickle
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri'
metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
contigs = ['10'] # 100, 1000
ploidies = ['3', '4', '6', '8', '10']
coverages = ['10', '20', '30', '40', '50']
results_dfs = []
for contig in contigs:
    for ploidy in ploidies:
        for coverage in coverages:
            results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results')
            if os.path.exists(results_path):
                samples = [f for f in os.listdir(results_path) if 'FFBS' in f]
                this_result_df = pd.DataFrame(columns=['contig', 'ploidy', 'coverage', 'sample', 'metric', 'value'], index=range(len(samples)*len(metrics)))
                this_result_df['contig'] = contig
                this_result_df['ploidy'] = ploidy
                this_result_df['coverage'] = coverage
                counter = 0
                for sample in samples:
                    sample_result = os.path.join(results_path, sample)
                    with open(sample_result, 'rb') as f:
                        this_results = pickle.load(f)
                    sample_name = sample.split('.pkl')[0].split('_')[-1]
                    evals = this_results['evaluation']
                    for metric in metrics:
                        this_result_df.loc[counter, 'sample'] = sample_name
                        this_result_df.loc[counter, 'metric'] = metric
                        this_result_df.loc[counter, 'value'] = evals[metric]
                        counter += 1
                results_dfs.append(this_result_df)

results_df = pd.concat(results_dfs, ignore_index=True)
results_df.to_csv(os.path.join(agg_results_path, 'sim_awri_results.csv'), index=False)


for metric in metrics:
    print('Metric:', metric)
    metric_df = results_df[results_df['metric'] == metric].reset_index(drop=True)
    sns.catplot(x="coverage", y="value",hue="ploidy", data=metric_df, kind="box", height=6, aspect=1.5)
    plt.title(f"{metric.capitalize()}")
    plt.xlabel("Coverage")
    plt.ylabel("Value")
    plt.tight_layout()
    plt.savefig(os.path.join(agg_results_path, f"{metric}.png"))
    plt.close()
