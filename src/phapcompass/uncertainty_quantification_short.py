import argparse, os, sys, time, pickle
from pHapCompass.combinatorials import *
from pHapCompass.FFBS import *
from pHapCompass.viterbi import *
from evaluations import *
from utils import *
from read_input import *
from matching import *
import pandas as pd
import numpy as np
from mec_based_predict_hap import *
from itertools import combinations
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


def make_inputs_for_uncertainty_allo():
    # data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_auto_short'
    data_path = '/labs/Aguiar/pHapCompass/dataset/simulated_data_allo_short_new'
    # output_dir = '/mnt/research/aguiarlab/proj/HaplOrbit/inputs_paper'
    output_dir = '/labs/Aguiar/pHapCompass/scripts/short_allo/input'
    # results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results_short'
    results_path = '/labs/Aguiar/pHapCompass/results/short_allo_new'
    # reference_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    reference_path = '/labs/Aguiar/pHapCompass/reference/simulated_haplotypes/allo'
    inputs = []
    ffb = 1
    mw = 10
    li = 1
    subgenome_configs = ['0.0005', '0.0001']
    within_mutation_rates = ['0.00005', '0.0001']
    ploidies = [3, 4, 6]
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C'} # , 8:'8ploid_4A+4B'}
    coverages = [5, 10, 20, 40]
    samples = range(10)
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for coverage in coverages:
                    for sample in samples:
                        this_results = os.path.join(results_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(1))
                        if not os.path.exists(this_results):
                            os.makedirs(this_results, exist_ok=True)
                            
                        frag_path = os.path.join(data_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(1) + '.frag')
                        bam_path = os.path.join(data_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(1) + '.bam')                        
                        vcf_path = os.path.join(reference_path, 'subgenome_config_mut' + str(sgc), ploidy_folder_name[ploidy], 'within_mut_' + str(mr), str(sample).zfill(1), 'Chr1.vcf')
                        genotype_path = os.path.join(reference_path, 'subgenome_config_mut' + str(sgc), ploidy_folder_name[ploidy], 'within_mut_' + str(mr), str(sample).zfill(1), 'genotype.csv')
                        result_file = os.path.join(this_results, f'viterbi_mec_based_single_weighted_mw{mw}_li{li}_ffbs{ffb}.pkl.gz')
                        this_inp = [frag_path, bam_path, vcf_path, genotype_path, 0.001, 0.0001, ploidy, this_results]
                        if not os.path.exists(result_file):
                            inputs.append(this_inp)
    print(len(inputs))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")


def make_inputs_for_uncertainty_auto():
    # data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_auto_short'
    data_path = '/labs/Aguiar/pHapCompass/dataset/simulated_data_auto_short'
    # output_dir = '/mnt/research/aguiarlab/proj/HaplOrbit/inputs_paper'
    output_dir = '/labs/Aguiar/pHapCompass/scripts/quant/input'
    # results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results_short'
    results_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto'
    # reference_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    reference_path = '/labs/Aguiar/pHapCompass/reference/simulated_haplotypes/auto'

    coverages = [3, 5, 10, 20]
    ploidies = [2, 3, 4, 6]
    mut_rates = [0.001, 0.005]
    samples = range(10)
    inputs = []
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                for sample in samples:
                    this_results = os.path.join(results_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2))
                    if not os.path.exists(this_results):
                        os.makedirs(this_results, exist_ok=True)
                    frag_path = os.path.join(data_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2) + '.frag')
                    bam_path = os.path.join(data_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2) + '.bam')
                    vcf_path = os.path.join(reference_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2), 'Chr1.vcf')
                    genotype_path = os.path.join(reference_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2), 'genotype.csv')
                    result_file = os.path.join(this_results, f'uncertainty_quant_ploidy{ploidy}_mut{mr}_coverage{coverage}_sample{str(sample).zfill(2)}.csv.gz')
                    # this_inp = [frag_path, bam_path, vcf_path, genotype_path, 0.001, 0.0001, ploidy, this_results]
                    this_inp = [frag_path, bam_path, vcf_path, genotype_path, 0.001, 0.0001, ploidy, result_file, coverage, mr, sample]

                    if not os.path.exists(result_file):
                       inputs.append(this_inp)
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
    results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results_quant/short_auto'
    coverages = [3, 5, 10, 20]
    ploidies = [2, 3, 4, 6]
    mut_rates = [0.001]
    samples = range(2)
    dfs = []
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                for sample in samples:
                    this_results = os.path.join(results_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2))
                    result_file = os.path.join(this_results, f'uncertainty_quant_ploidy{ploidy}_mut{mr}_coverage{coverage}_sample{str(sample).zfill(2)}.csv.gz')
                    if os.path.exists(result_file):
                        this_df = pd.read_csv(result_file)
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
    df_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/quant_uncert.csv'
    plot_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/plots'
    if not os.path.exists(plot_path):
        os.makedirs(plot_path, exist_ok=True)
    
    # Read data
    df = pd.read_csv(df_path, compression='gzip')
    
    # Get unique values for reporting
    unique_coverages = sorted(df['Coverage'].unique())
    unique_mut_rates = sorted(df['Mutation Rate'].unique())
    unique_ploidies = sorted(df['Ploidy'].unique())
    
    print(f"Unique coverages: {unique_coverages}")
    print(f"Unique mutation rates: {unique_mut_rates}")
    print(f"Unique ploidies: {unique_ploidies}")
    
    # ========== PLOT 1: Ploidy + Coverage (averaged across mutation rates) ==========
    # Group by Coverage, Ploidy, and genomic_distance, averaging across Mutation Rate
    grouped_cov_ploidy = df.groupby(['Coverage', 'Ploidy', 'genomic_distance'])['pair_phase_hamming_distance'].mean().reset_index()
    
    plt.figure(figsize=(14, 6))
    
    for cov in unique_coverages:
        for ploidy in unique_ploidies:
            subset = grouped_cov_ploidy[
                (grouped_cov_ploidy['Coverage'] == cov) & 
                (grouped_cov_ploidy['Ploidy'] == ploidy)
            ].sort_values('genomic_distance')
            
            if len(subset) > 0:
                label = f"cov={cov}, ploidy={ploidy}"
                plt.plot(subset['genomic_distance'], 
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
    grouped_mut_ploidy = df.groupby(['Mutation Rate', 'Ploidy', 'genomic_distance'])['pair_phase_hamming_distance'].mean().reset_index()
    
    plt.figure(figsize=(14, 6))
    
    for mut in unique_mut_rates:
        for ploidy in unique_ploidies:
            subset = grouped_mut_ploidy[
                (grouped_mut_ploidy['Mutation Rate'] == mut) & 
                (grouped_mut_ploidy['Ploidy'] == ploidy)
            ].sort_values('genomic_distance')
            
            if len(subset) > 0:
                label = f"mut_rate={mut}, ploidy={ploidy}"
                plt.plot(subset['genomic_distance'], 
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
    grouped_ploidy_only = df.groupby(['Ploidy', 'genomic_distance'])['pair_phase_hamming_distance'].mean().reset_index()
    
    # Custom colors for different ploidies
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

    plt.figure(figsize=(6, 4))

    for idx, ploidy in enumerate(unique_ploidies):
        subset = grouped_ploidy_only[
            grouped_ploidy_only['Ploidy'] == ploidy
        ].sort_values('genomic_distance')
        
        if len(subset) > 0:
            label = f"ploidy={ploidy}"
            plt.plot(
                subset['genomic_distance'],
                subset['pair_phase_hamming_distance'],
                label=label,
                linewidth=3,
                color=colors[idx % len(colors)],
                alpha=0.9
            )

    # ----- Bigger labels -----
    plt.xlabel('Distance between SNPs', fontsize=12)
    plt.ylabel('Mean Hamming Distance', fontsize=12)

    # ----- Bigger ticks -----
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    # ----- Bigger legend -----
    plt.legend(fontsize=14, loc='best')

    plt.grid(True, alpha=0.3)
    plt.ylim(bottom=0)
    plt.xlim(0, 150)
    plt.tight_layout()
    # plt.show()


    save_path_3 = os.path.join(plot_path, 'hamming_vs_distance_ploidy_only_auto_short.png')
    plt.savefig(save_path_3, dpi=300, bbox_inches='tight')
    print(f"✓ Plot 3 (Ploidy only) saved to {save_path_3}")
    plt.close()
    
    print("\n✓ All 3 plots generated successfully!")


def plot_uncertainty_bands():
    """
    Plot mean with confidence intervals showing uncertainty.
    """
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    
    df_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/quant_uncert.csv'
    plot_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/plots'
    os.makedirs(plot_path, exist_ok=True)
    
    df = pd.read_csv(df_path, compression='gzip')
    
    # Group and compute statistics
    stats = df.groupby(['Ploidy', 'genomic_distance'])['pair_phase_hamming_distance'].agg([
        'mean', 'std', 'count',
        ('q25', lambda x: x.quantile(0.25)),
        ('q75', lambda x: x.quantile(0.75))
    ]).reset_index()
    
    # Compute confidence intervals
    stats['ci_lower'] = stats['mean'] - 1.96 * stats['std'] / np.sqrt(stats['count'])
    stats['ci_upper'] = stats['mean'] + 1.96 * stats['std'] / np.sqrt(stats['count'])
    
    unique_ploidies = sorted(stats['Ploidy'].unique())
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Top plot: Mean with confidence intervals
    for idx, ploidy in enumerate(unique_ploidies):
        subset = stats[stats['Ploidy'] == ploidy].sort_values('genomic_distance')
        
        ax1.plot(subset['genomic_distance'], subset['mean'], 
                label=f"ploidy={ploidy}", color=colors[idx], linewidth=2)
        ax1.fill_between(subset['genomic_distance'], 
                         subset['ci_lower'], subset['ci_upper'],
                         alpha=0.2, color=colors[idx])
    
    ax1.set_xlabel('Distance between SNPs', fontsize=12)
    ax1.set_ylabel('Mean Hamming Distance', fontsize=12)
    ax1.set_title('Mean Hamming Distance with 95% Confidence Intervals', fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(bottom=0)
    
    # Bottom plot: Standard deviation (uncertainty itself)
    for idx, ploidy in enumerate(unique_ploidies):
        subset = stats[stats['Ploidy'] == ploidy].sort_values('genomic_distance')
        ax2.plot(subset['genomic_distance'], subset['std'], 
                label=f"ploidy={ploidy}", color=colors[idx], linewidth=2)
    
    ax2.set_xlabel('Distance between SNPs', fontsize=12)
    ax2.set_ylabel('Standard Deviation of Hamming Distance', fontsize=12)
    ax2.set_title('Uncertainty (Variability) in Phasing Quality', fontsize=14)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, 'uncertainty_with_bands.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Uncertainty bands plot saved")
    plt.close()


def plot_uncertainty_heatmap():
    """
    Heatmap showing uncertainty across ploidy and distance.
    """

    
    df_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/quant_uncert.csv'
    plot_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/plots'
    os.makedirs(plot_path, exist_ok=True)
    
    df = pd.read_csv(df_path, compression='gzip')
    
    # Bin distances for cleaner visualization
    df['distance_bin'] = pd.cut(df['genomic_distance'], bins=20)
    
    # Compute statistics
    heatmap_data = df.groupby(['Ploidy', 'distance_bin'])['pair_phase_hamming_distance'].agg([
        'mean', 'std',
        ('perfect_rate', lambda x: (x == 0).mean())
    ]).reset_index()
    
    # Create 3 heatmaps
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Pivot for heatmap format
    for idx, (metric, title) in enumerate([
        ('mean', 'Mean Hamming Distance'),
        ('std', 'Standard Deviation (Uncertainty)'),
        ('perfect_rate', 'Perfect Match Rate')
    ]):
        pivot = heatmap_data.pivot(index='Ploidy', columns='distance_bin', values=metric)
        
        sns.heatmap(pivot, ax=axes[idx], cmap='YlOrRd' if metric != 'perfect_rate' else 'RdYlGn',
                   cbar_kws={'label': metric}, annot=False)
        axes[idx].set_title(title, fontsize=12)
        axes[idx].set_xlabel('Genomic Distance', fontsize=11)
        axes[idx].set_ylabel('Ploidy', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, 'uncertainty_heatmap.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Uncertainty heatmap saved")
    plt.close()


def plot_error_distributions():
    """
    Show distribution of Hamming distances as violin plots.
    """

    
    df_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/quant_uncert.csv'
    plot_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/plots'
    os.makedirs(plot_path, exist_ok=True)
    
    df = pd.read_csv(df_path, compression='gzip')
    
    # Bin distances
    df['distance_category'] = pd.cut(df['genomic_distance'], 
                                     bins=[0, 50, 100, 200, 500],
                                     labels=['0-50', '50-100', '100-200', '200+'])
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    for idx, ploidy in enumerate(sorted(df['Ploidy'].unique())):
        subset = df[df['Ploidy'] == ploidy]
        
        sns.violinplot(data=subset, x='distance_category', y='pair_phase_hamming_distance',
                      ax=axes[idx], inner='box')
        axes[idx].set_title(f'Ploidy = {ploidy}', fontsize=12)
        axes[idx].set_xlabel('Distance Category', fontsize=11)
        axes[idx].set_ylabel('Hamming Distance', fontsize=11)
        axes[idx].grid(True, alpha=0.3, axis='y')
    
    plt.suptitle('Distribution of Phasing Errors by Distance', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, 'error_distributions.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Error distribution plot saved")
    plt.close()


def plot_sample_diversity():
    """
    Show diversity among FFBS samples - are they actually different?
    """

    
    df_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/quant_uncert.csv'
    plot_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/plots'
    os.makedirs(plot_path, exist_ok=True)
    
    df = pd.read_csv(df_path, compression='gzip')
    
    # Compute per-position statistics across samples
    position_stats = df.groupby(['Ploidy', 'ell_1', 'ell_2'])['pair_phase_hamming_distance'].agg([
        'mean', 'std', 'min', 'max',
        ('n_perfect', lambda x: (x == 0).sum()),
        ('n_samples', 'count')
    ]).reset_index()
    
    position_stats['diversity_score'] = position_stats['std']
    position_stats['confidence'] = position_stats['n_perfect'] / position_stats['n_samples']
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    unique_ploidies = sorted(position_stats['Ploidy'].unique())
    
    for idx, ploidy in enumerate(unique_ploidies):
        subset = position_stats[position_stats['Ploidy'] == ploidy]
        
        ax = axes[idx // 2, idx % 2]
        scatter = ax.scatter(subset['mean'], subset['diversity_score'], 
                           c=subset['confidence'], cmap='RdYlGn',
                           s=20, alpha=0.6, vmin=0, vmax=1)
        ax.set_xlabel('Mean Hamming Distance', fontsize=11)
        ax.set_ylabel('Diversity (Std Dev)', fontsize=11)
        ax.set_title(f'Ploidy = {ploidy}', fontsize=12)
        ax.grid(True, alpha=0.3)
        
        plt.colorbar(scatter, ax=ax, label='Confidence (% perfect)')
    
    plt.suptitle('Sample Diversity: Mean Error vs Uncertainty', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, 'sample_diversity.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Sample diversity plot saved")
    plt.close()


def plot_smoothed_trends():
    """
    Smoothed version with moving average to show clearer trends.
    """

    df_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/quant_uncert.csv'
    plot_path = '/archive/labs/aguiar/pHapCompass/results/uncertainty_quant/short_auto/plots'
    os.makedirs(plot_path, exist_ok=True)
    
    df = pd.read_csv(df_path, compression='gzip')
    
    grouped = df.groupby(['Ploidy', 'genomic_distance'])['pair_phase_hamming_distance'].mean().reset_index()
    
    unique_ploidies = sorted(grouped['Ploidy'].unique())
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    window = 20  # Moving average window
    
    for idx, ploidy in enumerate(unique_ploidies):
        subset = grouped[grouped['Ploidy'] == ploidy].sort_values('genomic_distance')
        
        # Raw data
        ax1.plot(subset['genomic_distance'], subset['pair_phase_hamming_distance'],
                label=f"ploidy={ploidy}", color=colors[idx], alpha=0.3, linewidth=1)
        
        # Smoothed data
        smoothed = subset['pair_phase_hamming_distance'].rolling(window=window, center=True).mean()
        ax2.plot(subset['genomic_distance'], smoothed,
                label=f"ploidy={ploidy}", color=colors[idx], linewidth=2.5)
    
    ax1.set_title('Raw Data (Noisy)', fontsize=12)
    ax1.set_xlabel('Distance between SNPs', fontsize=11)
    ax1.set_ylabel('Mean Hamming Distance', fontsize=11)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(bottom=0)
    
    ax2.set_title(f'Smoothed (Moving Average, window={window})', fontsize=12)
    ax2.set_xlabel('Distance between SNPs', fontsize=11)
    ax2.set_ylabel('Mean Hamming Distance', fontsize=11)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(bottom=0)
    
    plt.suptitle('Smoothed Trends Show Clearer Patterns', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, 'smoothed_trends.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Smoothed trends plot saved")
    plt.close()


def pair_uncertainty_quantification_from_ffbs(
    samples: List[Dict[int, Dict[str, str]]],
    true_haplotypes: pd.DataFrame,
    ploidy: int,
    coverage: float = None,
    mutation_rate: float = None,
    sample_id: int = None) -> Tuple[pd.DataFrame, Dict]:
    """
    Compute uncertainty quantification for pairwise positions in your graph.
    
    Args:
        samples: List of FFBS samples, each is {slice_id: {node: phasing_str}}
                 e.g., samples[0][5]['306-307'] = '0011'
        true_haplotypes: DataFrame with true haplotypes (K rows, L columns)
                         Columns are 0-indexed but represent 1-indexed SNP positions
        ploidy: Haplotype ploidy (K)
        coverage: Optional coverage value for output DataFrame
        mutation_rate: Optional mutation rate for output DataFrame
        sample_id: Optional sample ID for output DataFrame
        
    Returns:
        df: DataFrame with columns [Ploidy, Coverage, mutation_rate, sample, 
                                    ell_1, ell_2, ffbs_sample_iter, pair_phase_hamming_distance]
        cache: Dict mapping standardized phasing pairs -> distance (for efficiency)
    """

    # Extract unique nodes across all samples and slices
    all_nodes = set()
    for sample in samples:
        for slice_nodes in sample.values():
            all_nodes.update(slice_nodes.keys())
    
    # Parse node positions (convert '306-307' to tuple (306, 307))
    node_positions = {}
    for node in all_nodes:
        pos1, pos2 = node.split('-')
        node_positions[node] = (int(pos1), int(pos2))
    
    cache = {}
    df_data = []
    
    # For each node (pair of positions)
    for node, (pos1, pos2) in node_positions.items():
        # Extract true phasing for this pair (convert 1-indexed to 0-indexed)
        u = true_haplotypes.iloc[:, [pos1-1, pos2-1]].values  # (K, 2)
        u_sorted = sort_matrix_rows(u)
        u_key = tuple(map(tuple, u_sorted))
        
        # For each FFBS sample
        for ffbs_iter, sample in enumerate(samples):
            # Find this node in any slice of this sample (it's the same across slices)
            phasing_str = None
            for slice_nodes in sample.values():
                if node in slice_nodes:
                    phasing_str = slice_nodes[node]
                    break
            
            if phasing_str is None:
                # This node doesn't exist in this sample (shouldn't happen but handle it)
                continue
            
            # Convert string to matrix
            v = str_2_phas_1(phasing_str, ploidy)  # (K, 2)
            v_sorted = sort_matrix_rows(v)
            v_key = tuple(map(tuple, v_sorted))
            
            # Check cache (symmetric)
            cache_key = (u_key, v_key)
            reverse_key = (v_key, u_key)
            
            if cache_key in cache:
                distance = cache[cache_key]
            elif reverse_key in cache:
                distance = cache[reverse_key]
            else:
                # Compute and cache
                distance = optimal_hamming_distance(u, v)
                cache[cache_key] = distance
                cache[reverse_key] = distance
            
            # Append to results in colleague's format
            df_data.append({
                'Ploidy': ploidy,
                'Coverage': coverage,
                'mutation_rate': mutation_rate,
                'sample': sample_id,
                'n_snps': true_haplotypes.shape[1],  # NEW
                'n_nodes': None,  # NEW
                'n_edges': None,  # NEW
                'ell_1': pos1,
                'ell_2': pos2,
                'genomic_distance': pos2 - pos1,
                'ffbs_sample_iter': ffbs_iter,
                'pair_phase_hamming_distance': float(distance)
            })
    
    df = pd.DataFrame(df_data)
    return df, cache


def pair_uncertainty_quantification_from_full_haplotypes(
    predicted_haplotypes: List[np.ndarray],  # List of (K, L) arrays, may contain NaN
    true_haplotypes: pd.DataFrame,  # Ground truth (K, L)
    ploidy: int,
    coverage: float = None,
    mutation_rate: float = None,
    sample_id: int = None) -> Tuple[pd.DataFrame, Dict]:
    """
    Compute uncertainty quantification for all pairs of positions in full haplotypes.
    Similar to colleague's approach but handles NaN values from incomplete phasing.
    
    Args:
        predicted_haplotypes: List of m predicted haplotype matrices, each (K, L)
                              May contain NaN for unphased positions
        true_haplotypes: DataFrame with true haplotypes (K rows, L columns)
        ploidy: Haplotype ploidy (K)
        coverage: Optional coverage value
        mutation_rate: Optional mutation rate
        sample_id: Optional sample ID
        
    Returns:
        df: DataFrame with columns [Ploidy, Coverage, mutation_rate, sample, 
                                    ell_1, ell_2, genomic_distance, 
                                    ffbs_sample_iter, pair_phase_hamming_distance]
        cache: Dict mapping standardized phasing pairs -> distance
    """
    
    K, L = true_haplotypes.shape
    m = len(predicted_haplotypes)
    
    # Validate input shapes
    for i, pred_H in enumerate(predicted_haplotypes):
        assert pred_H.shape == (K, L), f"Sample {i} shape mismatch: {pred_H.shape} vs ({K}, {L})"
    
    # Convert true haplotypes to numpy
    H_true = true_haplotypes.values  # (K, L)
    
    cache = {}
    df_data = []
    
    # Process all pairs of positions
    valid_pairs = 0
    nan_pairs = 0
    
    for ell_1, ell_2 in combinations(range(L), 2):
        # True phasing for this pair
        u = H_true[:, [ell_1, ell_2]]  # (K, 2)
        u_sorted = sort_matrix_rows(u)
        u_key = tuple(map(tuple, u_sorted))
        
        # For each FFBS sample
        for ffbs_iter, pred_H in enumerate(predicted_haplotypes):
            # Extract predicted phasing for this pair
            v = pred_H[:, [ell_1, ell_2]]  # (K, 2)
            
            # Skip if either position has NaN (unphased)
            if np.any(np.isnan(v)):
                nan_pairs += 1
                continue
            
            valid_pairs += 1
            
            # Convert to int (in case they're floats)
            v = v.astype(int)
            
            # Sort for caching
            v_sorted = sort_matrix_rows(v)
            v_key = tuple(map(tuple, v_sorted))
            
            # Check cache
            cache_key = (u_key, v_key)
            reverse_key = (v_key, u_key)
            
            if cache_key in cache:
                distance = cache[cache_key]
            elif reverse_key in cache:
                distance = cache[reverse_key]
            else:
                distance = optimal_hamming_distance(u, v)
                cache[cache_key] = distance
                cache[reverse_key] = distance
            
            # Append result
            df_data.append({
                'Ploidy': ploidy,
                'Coverage': coverage,
                'Mutation Rate': mutation_rate,
                'Sample': sample_id,
                'ell_1': ell_1,  # 0-indexed position
                'ell_2': ell_2,  # 0-indexed position
                'genomic_distance': ell_2 - ell_1,
                'ffbs_sample_iter': ffbs_iter,
                'pair_phase_hamming_distance': float(distance)
            })
    
    df = pd.DataFrame(df_data)
    # # Summary statistics
    # total_possible_pairs = L * (L - 1) // 2
    # print(f"\n=== Full Haplotype Uncertainty Quantification ===")
    # print(f"Total SNP positions (L): {L}")
    # print(f"Total possible pairs: {total_possible_pairs}")
    # print(f"Valid comparisons: {valid_pairs} ({100*valid_pairs/(total_possible_pairs*m):.1f}% of all)")
    # print(f"Skipped (NaN): {nan_pairs} ({100*nan_pairs/(total_possible_pairs*m):.1f}% of all)")
    # print(f"Perfect match rate: {100 * (df['pair_phase_hamming_distance'] == 0).mean():.1f}%")
    # print(f"Mean Hamming distance: {df['pair_phase_hamming_distance'].mean():.2f}")
    return df, cache


def sort_matrix_rows(matrix):
    """Sort rows lexicographically. Columns stay in original order."""
    return matrix[np.lexsort(matrix.T[::-1])]


def optimal_hamming_distance(u, v):
    """
    Compute minimum Hamming distance between matrices u and v
    considering all possible row permutations.
    """
    
    K = u.shape[0]
    
    # Cost matrix: cost[i,j] = hamming distance between u[i] and v[j]
    cost = np.zeros((K, K))
    for i in range(K):
        for j in range(K):
            cost[i, j] = np.sum(u[i] != v[j])
    
    # Find optimal assignment
    row_ind, col_ind = linear_sum_assignment(cost)
    min_distance = cost[row_ind, col_ind].sum()
    
    return int(min_distance)


def compute_ffbs_sample_diversity(
    samples: List[Dict[int, Dict[str, str]]],
    ploidy: int) -> pd.DataFrame:
    """
    Compute pairwise Hamming distances between FFBS samples.
    This measures how different the samples are from EACH OTHER (not from ground truth).
    
    Returns:
        DataFrame with columns: [node, sample_i, sample_j, hamming_distance]
    """

    
    # Extract all nodes
    all_nodes = set()
    for sample in samples:
        for slice_nodes in sample.values():
            all_nodes.update(slice_nodes.keys())
    
    diversity_data = []
    cache = {}
    
    for node in all_nodes:
        # Extract phasings from all samples for this node
        phasings = []
        for sample_idx, sample in enumerate(samples):
            phasing_str = None
            for slice_nodes in sample.values():
                if node in slice_nodes:
                    phasing_str = slice_nodes[node]
                    break
            if phasing_str is not None:
                phasings.append((sample_idx, phasing_str))
        
        # Compute pairwise distances between samples
        for i in range(len(phasings)):
            for j in range(i + 1, len(phasings)):
                idx_i, str_i = phasings[i]
                idx_j, str_j = phasings[j]
                
                # Convert to matrices
                v_i = str_2_phas_1(str_i, ploidy)
                v_j = str_2_phas_1(str_j, ploidy)
                
                # Sort for caching
                v_i_sorted = sort_matrix_rows(v_i)
                v_j_sorted = sort_matrix_rows(v_j)
                key_i = tuple(map(tuple, v_i_sorted))
                key_j = tuple(map(tuple, v_j_sorted))
                
                cache_key = (key_i, key_j)
                reverse_key = (key_j, key_i)
                
                if cache_key in cache:
                    distance = cache[cache_key]
                elif reverse_key in cache:
                    distance = cache[reverse_key]
                else:
                    distance = optimal_hamming_distance(v_i, v_j)
                    cache[cache_key] = distance
                    cache[reverse_key] = distance
                
                diversity_data.append({
                    'node': node,
                    'sample_i': idx_i,
                    'sample_j': idx_j,
                    'hamming_distance': float(distance)
                })
    
    return pd.DataFrame(diversity_data)


def str_2_phas_1(phasing, ploidy):
    """Convert phasing string to numpy array."""
    return np.array([int(p) for p in [*phasing]]).reshape(ploidy, -1)


def uncertainty_quantification_short_model(inp):
    """
    End-to-end short-reads pipeline with canonical naming:
      - Nodes: '{i}-{j}' with 1-based SNP positions, i<j
      - Edges: '{i}-{j}--{j}-{k}' (topological)
    Uses vectorized emissions/transitions, forward pass, and FFBS sampler.
    """
    ffb = 1
    mw = 10
    li = 1

    this_frag_path, bam_file_path, vcf_file_path, genotype_path, error_rate, epsilon, ploidy, result_path, coverage, mr, sample = inp

    # ------------------- 1) Args + genotype (derive from VCF if needed) -------------------
    class Args:
        def __init__(self):
            self.vcf_path      = vcf_file_path
            self.data_path     = this_frag_path
            self.bam_path      = bam_file_path
            self.genotype_path = genotype_path
            self.ploidy        = ploidy
            self.error_rate    = error_rate
            self.epsilon       = epsilon
            self.result_path   = result_path
    
    args = Args()
    print(f'Working on {this_frag_path} ...')

    gen_df = pd.read_csv(genotype_path)
    true_haplotypes = gen_df.T
    n_snps = len(gen_df)

    # ------------------- 2) Read fragments (.frag) into sparse CSR ------------------------
    cfg  = InputConfigSparse(data_path=this_frag_path, genotype_path=genotype_path, ploidy=ploidy)
    frag = SparseFragment(cfg, positions_from_genotype=list(range(n_snps)))         # NOTE: _ingest_block must use START_IS_ONE_BASED=True
    data_matrix = frag.csr                    # (reads × SNPs) with {0,1(REF),2(ALT)}

    # SNP-index mappings:
    #   frag.col_to_snp[col] -> 0-based SNP position index
    #   frag.snp_to_col[pos] -> 0-based column index for 0-based position
    col_to_snp = np.array(frag.col_to_snp, dtype=np.int64)

    g = gen_df.sum(axis=1).to_numpy(dtype=np.int16)   # length = #SNP positions
    g_per_snp = g.astype(np.int16)

    # ------------------- 4) Build pair layer (nodes, counts4, edges, triple index) -------
    pair_layer = build_pair_layer(data_matrix, min_cocov=1)     # uses O^T O etc.

    # genotype per node, aligned by POSITIONS (not raw columns!)
    nodes_cols = pair_layer.nodes.astype(int)         # (P,2) column indices (0-based)
    gi = g_per_snp[col_to_snp[nodes_cols[:, 0]]]      # g at SNP position for col_i
    gj = g_per_snp[col_to_snp[nodes_cols[:, 1]]]
    genotype_pairs = list(zip(gi.tolist(), gj.tolist()))

    # ------------------- 5) Bank + node emissions (vectorized) ---------------------------
    bank = build_node_phasing_bank(ploidy, error_rate, genotype_pairs)
    state_names = build_state_names_from_bank(pair_layer, gi, gj, bank, ploidy, True, frag.col_to_snp)

    transitions_dict, transitions_extra = compute_transitions_optimized(data_matrix, pair_layer, state_names, ploidy, 
    error_rate,col_to_snp, frag.snp_to_col, return_extra=True)
    emissions = build_pair_emissions_from_state_names(state_names=state_names, ploidy=ploidy, error_rate=error_rate)

    # Build nodes/edges lists from canonical dictionaries
    nodes  = list(emissions.keys())
    edges  = [(u, v) for (u, v) in (k.split("--") for k in transitions_dict.keys())]

    # ------------------- 8) Slices (DAG layering) ----------------------------------------
    slices, _ = assign_slices_and_interfaces(nodes, edges)

    # (Optional) evidence indices per node/edge from M using position→column map
    assignment_dict = assign_evidence_to_states_and_transitions_from_M(nodes, edges, data_matrix, frag.snp_to_col)
    forward_messages = compute_forward_messages_from_M(slices,edges, assignment_dict, emissions, transitions_dict, data_matrix, frag.snp_to_col)
    predicted_haplotypes = []
    # samples = []
    for _ in range(10):
        this_sample = sample_states_book(slices, edges, forward_messages, transitions_dict)
        predicted_combined, block_ids = predict_haplotypes_mec_based(nodes, edges, this_sample, ploidy, genotype_path, data_matrix, 
        frag.snp_to_col, transitions_extra, args, priority="combined", likelihood_weight= li, mec_weight=mw, ffbs_weight=ffb, 
        allow_ffbs_override=True, verbose=False)
        # samples.append(this_sample)

        # Convert to numpy immediately if it's a DataFrame
        if isinstance(predicted_combined, pd.DataFrame):
            predicted_combined = predicted_combined.values
        
        predicted_haplotypes.append(predicted_combined)

    df, cache = pair_uncertainty_quantification_from_full_haplotypes(predicted_haplotypes, true_haplotypes, ploidy)

    df['Ploidy'] = ploidy
    df['Sample'] = sample
    df['Mutation Rate'] = mr
    df['Coverage'] = coverage
    df['n_snps'] = n_snps
    df['n_nodes'] = len(nodes)
    df['n_edges'] = len(edges)
    
    # pairwise_df, cache = pair_uncertainty_quantification_from_ffbs(
    #     samples=samples,  # Your 10 FFBS samples
    #     true_haplotypes=true_haplotypes,  # Ground truth DataFrame
    #     ploidy=ploidy)

    df.to_csv(result_path, compression="gzip", index=False)
    print('[Done]', result_path)
        
        
def run_pHapcompass_from_input(input_file):
    # input_file = '/mnt/research/aguiarlab/proj/HaplOrbit/inputs_paper/input_0.pkl'
    with open(input_file, "rb") as f:
        inp = pickle.load(f)

    uncertainty_quantification_short_model(inp)


if __name__ == '__main__':

    # collect_quant_results()
    # plot_mean_hamming_vs_distance()

    plot_uncertainty_bands()
    plot_uncertainty_heatmap()
    plot_error_distributions()
    plot_sample_diversity()
    plot_smoothed_trends()


    # if len(sys.argv) != 2:
    #     print("Usage: python3 simulator_paper.py <input_file>")
    #     sys.exit(1)
    
    # input_file = sys.argv[1]
    # run_pHapcompass_from_input(input_file)

    