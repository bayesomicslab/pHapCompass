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


def plot_samples_confidence():
    results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapcompass_simulated_NA12878_solutions_sampled.csv'

    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/plots/solution_samples'
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
    
    results_df = pd.read_csv(results_path)
    results_df = results_df.dropna()

    
    # Define color palette
    palette = sns.color_palette("tab10")

    # Loop through specific contigs and ploidies
    for contig in [100]:
        for ploidy in [3, 4, 6, 8]:
            print('Contig:', contig, 'Ploidy:', ploidy)
            
            # Subset data for this ploidy
            ploidy_df = results_df[(results_df['Ploidy'] == ploidy) & (results_df['Contig'] == contig)].reset_index(drop=True)
            df_pivot = ploidy_df.pivot(index=["Method", "Contig", "Ploidy", "Coverage", "Sample"], 
                     columns="Metric", values="Value").reset_index()
            # Create scatter plot with regression lines
            plt.figure(figsize=(8, 6))

            # Scatter plot (swapped axes)
            scatter = sns.scatterplot(data=df_pivot, x="Vector Error Rate", y="Confidence",
                                    hue="Coverage", palette=palette, s=10, edgecolor="black", alpha=0.5)

            # Regression lines (swapped axes) - NO LEGEND
            coverages = sorted(df_pivot['Coverage'].unique())  # Ensure same order for colors
            for i, coverage in enumerate(coverages):
                subset = df_pivot[df_pivot["Coverage"] == coverage]
                sns.regplot(data=subset, x="Vector Error Rate", y="Confidence",
                            scatter=False, color=palette[i], ci=95)  # Removed label to avoid extra legend

            # Title and labels (swapped)
            plt.title(f"Ploidy {ploidy}, Contig {contig}: Vector Error Rate vs. Confidence")
            plt.xlabel("Vector Error Rate")  
            plt.ylabel("Confidence")  

            # Keep only scatter plot legend
            handles, labels = scatter.get_legend_handles_labels()
            plt.legend(handles, labels, title="Coverage")

            # Save plot (optional)
            plt.savefig(os.path.join(plot_path, f"confidence_vs_vector_error_rate_ploidy_{ploidy}_contig_{contig}.png"), dpi=300, bbox_inches="tight")

            # Show plot
            # plt.show()


def plot_samples_entropy():
    results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapcompass_simulated_NA12878_solutions_sampled.csv'

    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/plots/solution_samples'
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
    
    results_df = pd.read_csv(results_path)
    results_df = results_df.dropna()

    # Define color palette
    palette = sns.color_palette("tab10")

    # Loop through specific contigs and ploidies
    for contig in [100]:
        for ploidy in [3, 4, 6, 8]:
            print('Contig:', contig, 'Ploidy:', ploidy)
            
            # Subset data for this ploidy
            ploidy_df = results_df[(results_df['Ploidy'] == ploidy) & (results_df['Contig'] == contig)].reset_index(drop=True)
            df_pivot = ploidy_df.pivot(index=["Method", "Contig", "Ploidy", "Coverage", "Sample"], 
                     columns="Metric", values="Value").reset_index()
            # Create scatter plot with regression lines
            plt.figure(figsize=(8, 6))

            # Scatter plot (swapped axes)
            scatter = sns.scatterplot(data=df_pivot, x="Vector Error Rate", y="Entropy",
                                    hue="Coverage", palette=palette, s=10, edgecolor="black", alpha=0.5)

            # Regression lines (swapped axes) - NO LEGEND
            coverages = sorted(df_pivot['Coverage'].unique())  # Ensure same order for colors
            for i, coverage in enumerate(coverages):
                subset = df_pivot[df_pivot["Coverage"] == coverage]
                sns.regplot(data=subset, x="Vector Error Rate", y="Entropy",
                            scatter=False, color=palette[i], ci=95)  # Removed label to avoid extra legend

            # Title and labels (swapped)
            plt.title(f"Ploidy {ploidy}, Contig {contig}: Vector Error Rate vs. Entropy")
            plt.xlabel("Vector Error Rate")  
            plt.ylabel("Entropy")  

            # Keep only scatter plot legend
            handles, labels = scatter.get_legend_handles_labels()
            plt.legend(handles, labels, title="Coverage")

            # Save plot (optional)
            plt.savefig(os.path.join(plot_path, f"entropy_vs_vector_error_rate_ploidy_{ploidy}_contig_{contig}.png"), dpi=300, bbox_inches="tight")

            # Show plot
            # plt.show()


def plot_viterbi_vs_samples():
    viterbi_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapcompass_simulated_NA12878_single_viterbi_true22.csv'
    samples_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapcompass_simulated_NA12878_solutions_sampled.csv'
    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/plots/viterbi22'
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
    viterbi_df = pd.read_csv(viterbi_path)
    samples_df = pd.read_csv(samples_path)
    viterbi_df = viterbi_df.dropna()
    samples_df = samples_df.dropna()
    samples_df = samples_df[samples_df['Metric'].isin(['Vector Error Rate', '# Phased Variants'])].reset_index(drop=True)

    samples_df.loc[samples_df['Method'] == 'Multi Variant', 'Method'] = 'FFBS + Multi'
    samples_df.loc[samples_df['Method'] == 'Single Variant', 'Method'] = 'FFBS + Single'

    all_results = pd.concat([viterbi_df, samples_df], ignore_index=True)

    method_palette = {
    "Viterbi + Single": "tab:blue",  # Blue
    "Viterbi + Multi": "tab:orange",  # Orange
    "FFBS + Single": "tab:green",  # Orange
    "FFBS + Multi": "tab:red"  # Purple}
    }


    # compare metrics for different methods:
    for contig in [100]:
        for metric in ['# Phased Variants', 'Vector Error Rate']:
            for ploidy in [3, 4, 6, 8]:
                
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
                
                # plt.show() 
                # Save the figure
                plt.savefig(os.path.join(plot_path, f"line_compare_{contig}_{metric}_{ploidy}.png"), bbox_inches="tight", dpi=300)
                plt.close()



def plot_autocorrelation():

    import os
    import numpy as np
    import matplotlib.pyplot as plt

    # Paths
    CLEANED_DIR = "/home/mah19006/projects/t2t/test_plot"
    PLOTS_DIR = os.path.join(CLEANED_DIR, "plots")
    os.makedirs(PLOTS_DIR, exist_ok=True)

    # Configuration
    CHROMOSOMES = ["chr3"]

    NONB_MOTIFS = ['A_Phased_Repeat', 'G_Quadruplex_Motif', 'Inverted_Repeat', 
                'Mirror_Repeat', 'Direct_Repeat', 'Short_Tandem_Repeat', 'Z_DNA_Motif']
    NONB_NAMES = ['A Phased Repeat', 'G Quadruplex', 'Inverted Repeat', 
                'Mirror Repeat', 'Direct Repeat', 'Short Tandem Repeat', 'Z DNA']

    # QUANTILES = [0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
    QUANTILES = [0.05, 0.25, 0.5, 0.75, 0.95, 0.99]

    COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'black']
    X_POS = np.arange(100)
    DIRECTIONS = ['same', 'opposite']
    DIRECTION_NAMES = {'same': 'Forward', 'opposite': 'Reverse'}


def quantile_plots_for_chrom(chrom):
    for i, motif in enumerate(NONB_MOTIFS):
        motif_pretty = NONB_NAMES[i]
        for dire in DIRECTIONS:
            print(motif, dire)
            nonb_file = os.path.join(CLEANED_DIR, f"noNaN_{chrom}_{motif}_{dire}_matrix.npz")
            control_file = os.path.join(CLEANED_DIR, f"noNaN_{chrom}_Control_{dire}_matrix.npz")
            
            if not (os.path.isfile(nonb_file) and os.path.isfile(control_file)):
                continue

            with np.load(nonb_file) as data:
                nonb_matrix = data[list(data.keys())[0]]
            print(motif, nonb_matrix.shape[0])
            with np.load(control_file) as data:
                control_matrix = data[list(data.keys())[0]]

            fig, ax = plt.subplots(figsize=(10, 7))

            for quant, color in zip(QUANTILES, COLORS):
                nonb_quantile = np.nanquantile(nonb_matrix, quant, axis=0)

                control_quantile = np.nanquantile(control_matrix, quant, axis=0)
                ax.plot(X_POS, nonb_quantile, color=color, linewidth=3, label=f"{quant}")
                ax.plot(X_POS, control_quantile, color=color, linestyle='--', linewidth=3)

            ax.set_xlabel('Position in Window', fontsize=18)
            ax.set_ylabel('Normalized Current', fontsize=18)
            ax.set_title(f"{chrom} {motif_pretty} ({DIRECTION_NAMES[dire]}) Quantile Plot", fontsize=20)

            quant_leg = ax.legend(loc=(1.03, 0.50), title="Quantile", prop={'size': 14}, title_fontsize=15)
            ax.add_artist(quant_leg)
            line_handles = [plt.Line2D([], [], color="gray", linestyle=ls, linewidth=3) for ls in ['-', '--']]
            ax.legend(handles=line_handles, labels=['Non-B DNA', 'Control'], loc=(1.03, 0.85),
                      title="Structure", prop={'size': 14}, title_fontsize=15)

            plt.setp(ax.get_xticklabels(), fontsize=13)
            plt.setp(ax.get_yticklabels(), fontsize=13)
            plt.tight_layout()

            plot_filename = os.path.join(PLOTS_DIR, f"{chrom}_{motif}_{DIRECTION_NAMES[dire]}_quantile.png")
            plt.savefig(plot_filename)
            plt.close()
            print(f"Saved plot: {plot_filename}")

def main():
    print("Generating quantile plots for chr3...")
    quantile_plots_for_chrom("chr3")
    print("Quantile plot generation complete.")





if __name__ == "__main__":
    main()
