import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import pysam
import os
import seaborn as sns
import graph_tool.all as gt
from scipy.spatial.distance import pdist, squareform


def snp_distance_distribution():
    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/plots'
    # input_vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/SRR942191/vcf_files/SRR942191_AHIQ01000001.1.vcf'
    input_vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/sassafras/vcf/filtered_CP142451.1_RagTag.vcf'
    
    vcf_file = pysam.VariantFile(input_vcf_path)
    original_snps = [record.pos - 1 for record in vcf_file.fetch()]
    original_snps = sorted(original_snps)
    snp_differences = np.diff(original_snps)

    # Plot the histogram
    plt.figure(figsize=(8, 4))
    # plt.hist(snp_differences, bins=20, edgecolor='black', alpha=0.75)
    plt.hist(snp_differences, edgecolor='black', alpha=0.75)
    # plt.title("Distribution of Distances Between SNP Sites", fontsize=16)
    plt.xlabel("Distance Between SNP Sites", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.yscale('log')
    plt.tight_layout()
    plt.show()
    plt.savefig(os.path.join(plot_path, 'snp_distance_histogram_AWRI.png'))

    vcf_file_human_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/hapref_chr21_filtered.vcf.bgz'
    vcf_file = pysam.VariantFile(vcf_file_human_path)
    original_snps = [record.pos - 1 for record in vcf_file.fetch()]

    original_snps = sorted(original_snps)
    snp_differences = np.diff(original_snps)
    snp_differences = [x for x in snp_differences if x > 0]

    # Plot the histogram
    plt.figure(figsize=(8, 6))
    plt.hist(snp_differences, bins=20, edgecolor='black', alpha=0.75)
    # plt.title("Distribution of Distances Between SNP Sites", fontsize=16)
    plt.xlabel("Distance Between SNP Sites", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.xscale('log')  # Set the x-axis to logarithmic scale
    # plt.show()
    plt.savefig(os.path.join(plot_path, 'snp_distance_histogram_hg19_chr21.png'))




    snp_differences = [x for x in snp_differences if x > 0]
    log_bins = np.logspace(np.log10(min(snp_differences)), np.log10(max(snp_differences)), 20)
    plt.figure(figsize=(8, 6))
    sns.histplot(snp_differences, bins=log_bins, kde=False, edgecolor='black', alpha=0.75)
    # Set the x-axis to logarithmic scale
    plt.xscale('log')

    # Customize the plot
    # plt.title("Distribution of Distances Between SNP Sites (Logarithmic Scale)", fontsize=16)
    plt.xlabel("Distance Between SNPs (log space)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Show the plot
    plt.tight_layout()
    plt.show()
    plt.savefig(os.path.join(plot_path, 'snp_distance_histogram_hg19_chr21.png'))


def plot_example_graph():
    # Initialize the graph
    g = gt.Graph(directed=True)

    # Add vertices (nodes)
    nodes = {
        '1-2': g.add_vertex(),
        '1-3': g.add_vertex(),
        '2-3': g.add_vertex(),
        '3-4': g.add_vertex(),
        '4-5': g.add_vertex(),
        '5-6': g.add_vertex(),
        '4-7': g.add_vertex(),
        '6-8': g.add_vertex(),
        '7-8': g.add_vertex()
    }

    # Define edges (directed)
    edges = [
        ('1-2', '1-3'), 
        ('1-2', '2-3'),
        ('2-3', '3-4'),
        ('1-3', '3-4'),
        ('1-3', '2-3'),
        ('4-5', '4-7'),
        ('3-4', '4-5'),
        ('3-4', '4-7'),
        ('4-5', '5-6'),
        ('5-6', '6-8'),
        ('4-7', '7-8'),
        ('6-8', '7-8')
    ]

    # Add edges to the graph
    for source, target in edges:
        g.add_edge(nodes[source], nodes[target])

    # Create a vertex property to label the nodes
    vertex_labels = g.new_vertex_property("string")
    for node, vertex in nodes.items():
        vertex_labels[vertex] = node

    # Draw the graph and save it as a PNG file
    # gt.graph_draw(
    #     g, 
    #     vertex_text=vertex_labels,
    #     vertex_font_size=12, 
    #     output_size=(800, 800)
    # )


    gt.graph_draw(
        g, 
        vertex_text=vertex_labels,
        vertex_font_size=18,  # Increase font size for node labels
        output_size=(800, 800),  # Increase output size to make edges appear longer
        layout="arf",  # Use the ARF layout for better spacing
        edge_pen_width=5 # Slightly thicker edges
    )


def plot_haplotypes_similarity_heatmap():
    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/plots'
    for ploidy in [3,4,6,8]:
        true_haplotypes_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878/contig_100/ploidy_{}/haplotypes.csv'.format(ploidy)

        true_haplotype = pd.read_csv(true_haplotypes_path).T.to_numpy()

        # Compute pairwise Hamming distances
        hamming_distances = pdist(true_haplotype, metric='hamming')  # Compute condensed distance matrix
        hamming_matrix = squareform(hamming_distances)  # Convert to square form
        hamming_similarity = 1 - squareform(hamming_distances)  # Convert to similarity (1 - distance)
        labels = [f"$h_{{{i+1}}}$" for i in range(true_haplotype.shape[0])]

        # Plot heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(hamming_similarity, annot=True, cmap="Blues", xticklabels=labels, yticklabels=labels, annot_kws={"size": 16})
        # plt.title("Pairwise Hamming Distance Heatmap")
        plt.xlabel("Haplotypes", fontsize=26, labelpad=10)
        plt.ylabel("Haplotypes", fontsize=26, labelpad=10)
        plt.xticks(rotation=0, fontsize=22)  # Rotate x-tick labels for better readability
        plt.yticks(rotation=0, fontsize=22)   # Keep y-tick labels horizontal
        # plt.show()
        plt.savefig(os.path.join(plot_path, f'hamming_similarity_ploidy_{ploidy}.png'))


def plot_sampling_accuracy():

    pairwise_samples_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pairwise_samples_accuracy.csv'
    pairwise_samples = pd.read_csv(pairwise_samples_path)
    
    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/plots/'

    method_palette = {
    "FFBS": "tab:blue",  # Blue
    "Likelihood": "tab:orange",  # Orange
    "Beliefs": "tab:green", 
    "LBP Samples": "tab:red",  # Red
    "LBP + FFBS": "tab:purple"  # Purple}
    }
    # compare metrics for different methods:
    for contig in [100]:
        for ploidy in [3, 4, 6, 8]:
            print('Contig:', contig, 'Ploidy:', ploidy)
            metric_df = pairwise_samples[(pairwise_samples['Ploidy'] == ploidy)].reset_index(drop=True)
            
            # Create figure and axis
            plt.figure(figsize=(6, 6))
            ax = sns.lineplot(x="Coverage", y="Value", hue="Sampling", data=metric_df, marker="o", linewidth=4, palette=method_palette)

            # Add the title and labels
            plt.title(f"Ploidy: {str(ploidy).capitalize()}", y=1.05, fontdict={"size": 20})
            plt.xlabel("Coverage", fontdict={"size": 20})
            # plt.ylabel("Value", fontdict={"size": 20})
            plt.ylabel('Sampling Accuracy', fontdict={"size": 20})

            # Ensure x-ticks match only the unique Coverage values
            unique_coverage = sorted(metric_df["Coverage"].unique())
            plt.xticks(unique_coverage)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            # Move the legend to the top-right, inside the plot area
            plt.legend(title="Pairwise Sampling", loc="upper right", bbox_to_anchor=(0.95, 0.8), frameon=True, fontsize=16, title_fontsize=20)

            # Adjust layout
            plt.tight_layout()

            # Save the figure
            plt.savefig(os.path.join(plot_path, f"Sampling_acc_{contig}_{ploidy}.png"), bbox_inches="tight", dpi=300)
            plt.close()


def plot_metric_methods_line(all_results, agg_results_path):
    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/plots/NA12878/'
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)

    all_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapcompass_simulated_NA12878_count_evidence_mv.csv'
    all_results = pd.read_csv(all_results_path)
    method_palette = {
    "pHapCompass": "tab:blue",  # Blue
    'pHapCompass + LBP': "tab:green",
    'pHapCompass + M.V.': "tab:red", 
    "WhatsHap": "tab:orange",  # Orange
    "H-PoPG": "tab:purple", 
    'Edge count': "tab:blue", 
    "Evidence count": "tab:orange",  # Orange
    "M.V.": "tab:green", }
    
    # compare metrics for different methods:
    for contig in [100]:
        for metric in ['# Phased Variants', 'Vector Error Rate']:
            for ploidy in [2, 3, 4, 6, 8]:
                print('Contig:', contig, 'Metric:', metric, 'Ploidy:', ploidy)
                metric_df = all_results[(all_results['Metric'] == metric) & (all_results['Ploidy'] == ploidy)].reset_index(drop=True)
                
                # Create figure and axis
                plt.figure(figsize=(6, 6))
                ax = sns.lineplot(x="Coverage", y="Value", hue="Method", data=metric_df, marker="o", linewidth=4, palette=method_palette)

                # Add the title and labels
                plt.title(f"Ploidy: {str(ploidy).capitalize()}", y=1.05, fontdict={"size": 20})
                plt.xlabel("Coverage", fontdict={"size": 20})
                # plt.ylabel("Value", fontdict={"size": 20})
                plt.ylabel(metric, fontdict={"size": 20})

                # Ensure x-ticks match only the unique Coverage values
                unique_coverage = sorted(metric_df["Coverage"].unique())
                plt.xticks(unique_coverage)
                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)
                # Move the legend to the top-right, inside the plot area
                plt.legend(title="Method", loc="upper right", bbox_to_anchor=(0.95, 0.9), frameon=True, fontsize=16, title_fontsize=20)

                # Adjust layout
                plt.tight_layout()

                # Save the figure
                plt.savefig(os.path.join(plot_path, f"line_compare_{contig}_{metric}_{ploidy}.png"), bbox_inches="tight", dpi=300)
                plt.close()


def plot_metric_methods_line_broken_axis(all_results, agg_results_path):

    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/plots/multi_variants/'
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)

    all_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/all_methods_results_LBP_multi_included.csv'
    all_results = pd.read_csv(all_results_path)
    method_palette = {
    "pHapCompass": "tab:blue",  # Blue
    'pHapCompass + LBP': "tab:green",
    'pHapCompass + M.V.': "tab:red", 
    "WhatsHap": "tab:orange",  # Orange
    "H-PoPG": "tab:purple"}
    
    # Function to determine y-axis breaks
    def determine_breaks(values):
        """
        Identify where to break the y-axis.
        The break happens between a small positive value (fixed at 0.05) and 
        just below the minimum actual data value.
        """
        values = np.array(values)
        min_val = values.min()
        max_val = values.max()
        
        lower_break = 0.01  # A small positive value just above zero
        upper_break = min_val * 0.9  # Just below min_val

        return lower_break, upper_break, max_val


    # Loop through ploidy and metrics
    for contig in [100]:
        for metric in ['# Phased Variants', 'Vector Error Rate']:
            for ploidy in [3, 4, 6, 8]:
                print('Contig:', contig, 'Metric:', metric, 'Ploidy:', ploidy)
                metric_df = all_results[(all_results['Metric'] == metric) & (all_results['Ploidy'] == ploidy)].reset_index(drop=True)
                
                if metric_df.empty:
                    print(f"Skipping ploidy {ploidy} for {metric} due to empty data.")
                    continue

                # Determine y-axis breakpoints
                lower_break, upper_break, max_value = determine_breaks(metric_df["Value"])

                # Create subplots with two axes for broken y-axis
                fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6), 
                                            gridspec_kw={'height_ratios': [1, 3]})
                fig.subplots_adjust(hspace=0.1)  # Slightly increase space between axes

                # Plot on both axes
                sns.lineplot(x="Coverage", y="Value", hue="Method", data=metric_df, 
                            marker="o", linewidth=4, palette=method_palette, ax=ax1)
                sns.lineplot(x="Coverage", y="Value", hue="Method", data=metric_df, 
                            marker="o", linewidth=4, palette=method_palette, ax=ax2)

                # Adjust y-limits
                ax1.set_ylim(upper_break, max_value)  # Upper part for main data
                ax2.set_ylim(0, lower_break)  # Lower part for near-zero values

                # Hide spines between axes
                ax1.spines.bottom.set_visible(False)
                ax2.spines.top.set_visible(False)
                ax1.xaxis.tick_top()
                ax1.tick_params(labeltop=False)  # Hide upper labels
                ax2.xaxis.tick_bottom()

                # Add break marks
                d = 0.5  # Proportion of vertical to horizontal extent of the slanted line
                kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none", 
                            color='k', mec='k', mew=1, clip_on=False)
                ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
                ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

                # Add title above both axes
                fig.suptitle(f"Ploidy: {str(ploidy).capitalize()}", fontsize=20, y=0.98)

                # Set axis labels
                ax2.set_xlabel("Coverage", fontsize=20)
                ax2.set_ylabel(metric, fontsize=20)

                # Ensure x-ticks match unique Coverage values
                unique_coverage = sorted(metric_df["Coverage"].unique())
                ax2.set_xticks(unique_coverage)
                ax2.tick_params(axis='x', labelsize=16)  # X-tick font size
                ax1.tick_params(axis='y', labelsize=16)  # Match y-tick font sizes
                ax2.tick_params(axis='y', labelsize=16)

                # Move the legend to the bottom axis (ax2 only)
                ax2.legend(title="Method", loc="upper right", bbox_to_anchor=(0.95, 0.9), 
                        frameon=True, fontsize=16, title_fontsize=20)

                # Save the figure
                plt.savefig(os.path.join(plot_path, f"line_compare_{contig}_{metric}_{ploidy}.png"), 
                            bbox_inches="tight", dpi=300)
                plt.close()


def plot_FFBS_entropy(all_results, agg_results_path):
    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/plots/ffbs_entropies/'
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)

    all_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/ffbs_simulated_NA12878.csv'
    all_results = pd.read_csv(all_results_path)
    # all_results.fillna(0, inplace=True)

    method_palette = {
    "entropy_correct_mean": "tab:green",  # Blue
    'entropy_wrong_mean': "tab:red",
    'confidence_correct_mean': "tab:green", 
    "confidence_wrong_mean": "tab:red"  # Orange
    }
    metrics = ['entropy', 'confidence']

    # compare metrics for different methods:
    for contig in [100]:
        for metric in metrics:
            for ploidy in [3, 4, 6, 8]:
                print('Contig:', contig, 'Metric:', metric, 'Ploidy:', ploidy)
                included_metrics = [metric + '_correct_mean', metric + '_wrong_mean']
                metric_df = all_results[(all_results['Metric'].isin(included_metrics)) & (all_results['Ploidy'] == ploidy)].reset_index(drop=True)
                
                # Create figure and axis
                plt.figure(figsize=(12, 6))
                ax = sns.boxplot(x="Coverage", y="Value", hue="Metric", data=metric_df, palette=method_palette)
                
                # Add the title and labels
                plt.title(f"Ploidy: {str(ploidy).capitalize()}", y=1.05, fontdict={"size": 20})
                plt.xlabel("Coverage", fontdict={"size": 20})
                # plt.ylabel("Value", fontdict={"size": 20})
                plt.ylabel(metric, fontdict={"size": 20})

                # Ensure x-ticks match only the unique Coverage values
                # unique_coverage = sorted(metric_df["Coverage"].unique())
                # plt.xticks(unique_coverage)
                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)
                # Move the legend to the top-right, inside the plot area
                plt.legend(title="Method", loc="upper right", bbox_to_anchor=(0.95, 0.9), frameon=True, fontsize=16, title_fontsize=20)

                # Adjust layout
                plt.tight_layout()
                
                # plt.show()

                # Save the figure
                plt.savefig(os.path.join(plot_path, f"ffbs_{contig}_{metric}_{ploidy}.png"), bbox_inches="tight", dpi=300)
                plt.close()

