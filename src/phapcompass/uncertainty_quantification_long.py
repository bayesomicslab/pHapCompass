from .evaluations import pair_uncertainty_quantification
import gzip
import pickle
from .long.crf_gibbs_log_space_vectorized import *
import os
import pandas as pd
import matplotlib.pyplot as plt


def make_inputs_for_create_pair_phase_df():

    results_path = '/archive/labs/aguiar/pHapCompass/results/results_long_auto'
    save_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/long_auto'
    output_dir = '/labs/Aguiar/pHapCompass/scripts/long_quant/input'

    if not os.path.exists(save_path):
        os.makedirs(save_path)

    coverages = [3, 5, 10, 20]
    ploidies = [2, 3, 4, 6]
    mut_rates = [0.001, 0.005]
    samples = range(10)
    inputs = []
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                for sample in samples:
                    this_results = os.path.join(results_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2), 'delta5_ep1e-05_lr0.02.pkl.gz')
                    this_save_file = os.path.join(save_path, f"uncertainty_quant_ploidy{ploidy}_mut_{mr}_cov_{coverage}_sample{str(sample).zfill(2)}.csv.gz")
                    if os.path.exists(this_results): # and os.path.exists(this_save_file):
                        this_inp = [this_results, this_save_file, ploidy, coverage, mr, sample]
                        inputs.append(this_inp)
                        print(this_save_file)

    print(len(inputs))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")


def collect_quant_results():
    # results_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto'
    results_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/long_auto'
    coverages = [3, 5, 10, 20]
    ploidies = [2, 3, 4, 6]
    mut_rates = [0.001]
    samples = range(2)
    dfs = []
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                for sample in samples:
                    # this_results = os.path.join(results_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2))
                    result_file = os.path.join(results_path, f'uncertainty_quant_ploidy{ploidy}_mut_{mr}_cov_{coverage}_sample{str(sample).zfill(2)}.csv.gz')
                    if os.path.exists(result_file):
                        this_df = pd.read_csv(result_file, compression="gzip")
                        dfs.append(this_df)
    all_df = pd.concat(dfs)
    print(len(all_df))
    all_df.to_csv(os.path.join(results_path, 'quant_uncert.csv'), compression="gzip", index=False)
    print('saved:', os.path.join(results_path, 'quant_uncert.csv'))


def plot_mean_hamming_vs_distance():
    """
    Plot mean Hamming distance vs genomic distance with proper aggregation.
    Creates three separate plots:
    1. Ploidy + Coverage (averaged across all mutation rates)
    2. Ploidy + Mutation Rate (averaged across all coverages)  
    3. Ploidy only (averaged across all coverages and mutation rates)
    """
    
    # Paths
    df_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/long_auto/quant_uncert.csv'
    plot_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/long_auto/plots'
    if not os.path.exists(plot_path):
        os.makedirs(plot_path, exist_ok=True)
    
    # Read data
    df = pd.read_csv(df_path, compression='gzip')
    df["snp difference"] = df["ell_2"] - df["ell_1"]

    # Get unique values for reporting
    unique_coverages = sorted(df['Coverage'].unique())
    unique_mut_rates = sorted(df['mutation_rate'].unique())
    unique_ploidies = sorted(df['Ploidy'].unique())
    
    print(f"Unique coverages: {unique_coverages}")
    print(f"Unique mutation rates: {unique_mut_rates}")
    print(f"Unique ploidies: {unique_ploidies}")
    
    # ========== PLOT 1: Ploidy + Coverage (averaged across mutation rates) ==========
    # Group by Coverage, Ploidy, and genomic_distance, averaging across Mutation Rate
    grouped_cov_ploidy = df.groupby(['Coverage', 'Ploidy', 'snp difference'])['pair_phase_hamming_distance'].mean().reset_index()
    
    plt.figure(figsize=(14, 6))
    
    for cov in unique_coverages:
        for ploidy in unique_ploidies:
            subset = grouped_cov_ploidy[
                (grouped_cov_ploidy['Coverage'] == cov) & 
                (grouped_cov_ploidy['Ploidy'] == ploidy)
            ].sort_values('snp difference')
            
            if len(subset) > 0:
                label = f"cov={cov}, ploidy={ploidy}"
                plt.plot(subset['snp difference'], 
                        subset['pair_phase_hamming_distance'], 
                        label=label, 
                        linewidth=2,
                        alpha=0.8)
    
    plt.xlabel('Distance between SNPs', fontsize=12)
    plt.ylabel('Mean Hamming Distance', fontsize=12)
    plt.title('Mean Hamming Distance vs Distance between SNPs\n(Averaged across all mutation rates)', fontsize=14)
    plt.legend(fontsize=9, loc='best', ncol=2)
    plt.grid(True, alpha=0.3)
    plt.ylim(bottom=0)
    plt.tight_layout()
    
    save_path_1 = os.path.join(plot_path, 'hamming_vs_distance_coverage_ploidy.png')
    plt.savefig(save_path_1, dpi=300, bbox_inches='tight')
    print(f"✓ Plot 1 (Coverage + Ploidy) saved to {save_path_1}")
    plt.close()
    
    # ========== PLOT 2: Ploidy + Mutation Rate (averaged across coverages) ==========
    # Group by Mutation Rate, Ploidy, and genomic_distance, averaging across Coverage
    grouped_mut_ploidy = df.groupby(['mutation_rate', 'Ploidy', 'snp difference'])['pair_phase_hamming_distance'].mean().reset_index()
    
    plt.figure(figsize=(14, 6))
    
    for mut in unique_mut_rates:
        for ploidy in unique_ploidies:
            subset = grouped_mut_ploidy[
                (grouped_mut_ploidy['mutation_rate'] == mut) & 
                (grouped_mut_ploidy['Ploidy'] == ploidy)
            ].sort_values('snp difference')
            
            if len(subset) > 0:
                label = f"mut_rate={mut}, ploidy={ploidy}"
                plt.plot(subset['snp difference'], 
                        subset['pair_phase_hamming_distance'], 
                        label=label, 
                        linewidth=2,
                        alpha=0.8)
    
    plt.xlabel('Distance between SNPs', fontsize=12)
    plt.ylabel('Mean Hamming Distance', fontsize=12)
    plt.title('Mean Hamming Distance vs Distance between SNPs\n(Averaged across all coverages)', fontsize=14)
    plt.legend(fontsize=9, loc='best', ncol=2)
    plt.grid(True, alpha=0.3)
    plt.ylim(bottom=0)
    plt.tight_layout()
    
    save_path_2 = os.path.join(plot_path, 'hamming_vs_distance_mutation_ploidy.png')
    plt.savefig(save_path_2, dpi=300, bbox_inches='tight')
    print(f"✓ Plot 2 (Mutation Rate + Ploidy) saved to {save_path_2}")
    plt.close()
    
    # ========== PLOT 3: Only Ploidy (averaged across ALL coverages and mutation rates) ==========
    # Group by Ploidy and genomic_distance only, averaging across Coverage and Mutation Rate
    grouped_ploidy_only = df.groupby(['Ploidy', 'snp difference'])['pair_phase_hamming_distance'].mean().reset_index()
    
    # Custom colors for different ploidies
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    plt.figure(figsize=(6, 4))
    
    for idx, ploidy in enumerate(unique_ploidies):
        subset = grouped_ploidy_only[
            grouped_ploidy_only['Ploidy'] == ploidy
        ].sort_values('snp difference')
        
        if len(subset) > 0:
            label = f"ploidy={ploidy}"
            plt.plot(subset['snp difference'], 
                    subset['pair_phase_hamming_distance'], 
                    label=label, 
                    linewidth=3,
                    color=colors[idx % len(colors)],
                    alpha=0.9)
    
    plt.xlabel('Distance between SNPs', fontsize=12)
    plt.ylabel('Mean Hamming Distance', fontsize=12)
    # plt.title('Mean Hamming Distance vs Distance between SNPs\n(Averaged across all coverages and mutation rates)', fontsize=14)
    plt.legend(fontsize=14, loc='best')
    plt.grid(True, alpha=0.3)
    plt.ylim(bottom=0)
    plt.tight_layout()
    # plt.show()

    save_path_3 = os.path.join(plot_path, 'hamming_vs_distance_ploidy_only.png')
    plt.savefig(save_path_3, dpi=300, bbox_inches='tight')
    print(f"✓ Plot 3 (Ploidy only) saved to {save_path_3}")
    plt.close()
    
    print("\n✓ All 3 plots generated successfully!")


def create_pair_phase_df(inp):
    result_file, save_path, ploidy, coverage, mr, sample = inp
    
    columns = ['Ploidy', 'Coverage', 'mutation_rate', 'sample', 'ell_1', 'ell_2', "pair_phase_hamming_distance"]
    df = pd.DataFrame(columns=columns, index=range(2000000))
    counter = 0
    with gzip.open(result_file, "rb") as f:
        result_pkl = pickle.load(f)
    true_haplotypes = result_pkl["true_haplotypes"]
    sampler = result_pkl["sampler"]
    ffbs_samples = sampler._ffbs_decode_joint_phase_multiple(10)
    results, cache = pair_uncertainty_quantification(true_haplotypes, ffbs_samples)
    for k, v in results.items():
        ell_1, ell_2 = k
        # v is now a single aggregated value over all ffbs samples
        df.loc[counter, :] = ploidy, coverage, mr, sample, ell_1, ell_2, v
        counter += 1
    print(len(df))
    df = df.dropna()
    df.to_csv(save_path, compression="gzip", index=False)
    print(f"[Done] Completed ploidy {ploidy}, mut_rate {mr}, coverage {coverage}, sample {sample}", len(df))


def plot_long_read_model_pair_phase_uncertainty_quantification_plot():
    ploidy = 2
    df_gzipped_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/long_auto/quant_uncert.csv'
    plot_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/long_auto/plots'
    df = pd.read_csv(df_gzipped_path, compression="gzip")
    print(df.columns.values)
    df = df[df["Ploidy"] == ploidy]
    
    # Find max ell_2 for each unique combination
    max_ell2_per_combo = df.groupby(['sample', 'Coverage', 'Ploidy', 'mutation_rate'])['ell_2'].max()

    # Find the smallest value among those maximums
    smallest_max_ell2 = max_ell2_per_combo.min()
    
    df_trunc = df[df["ell_2"] <= smallest_max_ell2]
    
    df_trunc["snp difference"] = df_trunc["ell_2"] - df_trunc["ell_1"]
    
    summary_df = df_trunc.groupby(['Ploidy', 'Coverage', 'mutation_rate', 'snp difference']).agg(
        mean_hamming_distance=('pair_phase_hamming_distance', 'mean'),
        std_hamming_distance=('pair_phase_hamming_distance', 'std'),
        count=('pair_phase_hamming_distance', 'count')
    ).reset_index()
    
    # Get unique combinations of coverage and mutation_rate
    unique_combos = summary_df[['Coverage', 'mutation_rate']].drop_duplicates()

    # Create the plot
    plt.figure(figsize=(10, 6))

    for _, row in unique_combos.iterrows():
        cov = row['Coverage']
        mut = row['mutation_rate']
        
        # Filter data for this combination
        subset = summary_df[(summary_df['Coverage'] == cov) & 
                            (summary_df['mutation_rate'] == mut)]
        
        # Sort by distance for proper line plotting
        subset = subset.sort_values('snp difference')
        
        # Plot the line
        label = f'cov={cov}, mut_rate={mut}'
        plt.plot(subset['snp difference'], 
                subset['mean_hamming_distance'], 
                label=label)

    plt.xlabel('Distance between SNPs')
    plt.ylabel('Mean Hamming Distance')
    plt.title('Mean Hamming Distance vs Distance between SNPs')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    # plt.savefig()

    save_path = os.path.join(plot_path, 'hamming_vs_distance_ploidy_plot.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

  
def create_pair_phase_df_from_input(input_file):
    # input_file = '/mnt/research/aguiarlab/proj/HaplOrbit/inputs_paper/input_0.pkl'
    with open(input_file, "rb") as f:
        inp = pickle.load(f)

    create_pair_phase_df(inp)


if __name__ == '__main__':

    # if len(sys.argv) != 2:
    #     print("Usage: python3 simulator_paper.py <input_file>")
    #     sys.exit(1)
    
    # input_file = sys.argv[1]
    # create_pair_phase_df_from_input(input_file)

    # collect_quant_results()
    plot_long_read_model_pair_phase_uncertainty_quantification_plot()
    plot_mean_hamming_vs_distance()
