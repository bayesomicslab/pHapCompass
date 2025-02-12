import numpy as np
from matplotlib import pyplot as plt
import pysam
import os
import seaborn as sns
import graph_tool.all as gt


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


# Call the function to plot the graph
