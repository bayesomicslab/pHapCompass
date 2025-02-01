import os
import pickle
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from utils.utils import sort_nodes, str_2_phas_1
import numpy as np
import itertools
from collections import defaultdict
import pysam


def prepare_results_ismb():
    agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test'
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    contigs = ['100'] 
    ploidies = ['6']
    coverages = ['10', '50', '100']
    results_dfs = []
    for contig in contigs:
        for ploidy in ploidies:
            for coverage in coverages:
                # results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results')
                results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results_likelihood')
                if os.path.exists(results_path):
                    samples = [f for f in os.listdir(results_path) if 'FFBS' in f]
                    this_result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value', 'length_phased'], index=range(len(samples)*len(metrics)))
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
                        phased_snp = this_results['predicted_haplotypes'].shape[1]
                        for metric in metrics:
                            this_result_df.loc[counter, 'Sample'] = sample_name
                            this_result_df.loc[counter, 'Metric'] = metric
                            this_result_df.loc[counter, 'Value'] = evals[metric]
                            this_result_df.loc[counter, 'length_phased'] = phased_snp
                            counter += 1
                    results_dfs.append(this_result_df)

    results_df = pd.concat(results_dfs, ignore_index=True)
    results_df.to_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_test.csv'), index=False)


    # server3_results = pd.read_csv(os.path.join(agg_results_path, 'sim_awri_results_likelihood10_100_server3.csv'))
    beagle_results = pd.read_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_test.csv'))
    # ser3_beagle = pd.concat([server3_results, beagle_results], ignore_index=True)
    # ser3_beagle_non_redundant = ser3_beagle.groupby(["Contig", "Ploidy", "Coverage"], group_keys=False).apply(lambda group: group.drop_duplicates(subset="Sample")).reset_index(drop=True)
    # ser3_beagle_non_redundant = ser3_beagle_non_redundant[ser3_beagle_non_redundant['Ploidy'] != 8].reset_index(drop=True)
    # non_redundant_sorted_df = (
    #     ser3_beagle_non_redundant.groupby(["Contig", "Ploidy", "Coverage"], group_keys=False)
    #     .apply(lambda group: group.drop_duplicates(subset="Sample").sort_values(by="Sample"))
    #     .reset_index(drop=True)
    # )
    group_sizes = beagle_results.groupby(["Contig", "Ploidy", "Coverage"]).size()
    small_groups = group_sizes[group_sizes < 100]


    capital = {'method': 'Method', 'contig': 'Contig', 'ploidy': 'Ploidy', 'coverage': 'Coverage', 'sample': 'Sample', 'metric': 'Metric', 'value': 'Value'}
    hpopg_sim10 = pd.read_csv(os.path.join(agg_results_path, 'hpopg_sim_awri_results10_346.csv'))
    hpopg_sim100 = pd.read_csv(os.path.join(agg_results_path, 'hpopg_sim_awri_results100_3_4_6.csv'))
    hpopg_sim10 = hpopg_sim10.rename(columns=capital)
    hpopg_sim100 = hpopg_sim100.rename(columns=capital)
    whatshapp_10 = pd.read_csv(os.path.join(agg_results_path, 'whatshap_sim_awri_results_10_3_4_6.csv'))
    whatshapp_100 = pd.read_csv(os.path.join(agg_results_path, 'whatshap_sim_awri_results_100_3_4_6.csv'))
    whatshapp_10 = whatshapp_10.rename(columns=capital)
    whatshapp_100 = whatshapp_100.rename(columns=capital)

    all_results = pd.concat([beagle_results, hpopg_sim10, hpopg_sim100, whatshapp_10, whatshapp_100], ignore_index=True)
    all_results['Contig'] = all_results['Contig'].astype(int)
    all_results['Ploidy'] = all_results['Ploidy'].astype(int)
    all_results['Coverage'] = all_results['Coverage'].astype(int)
    all_results['length_phased'] = all_results['length_phased'].astype(int)

    all_results = all_results[all_results['Ploidy'] != 8].reset_index(drop=True)
    all_results = all_results[all_results['Metric'].isin(['vector_error_rate', 'mec', 'mismatch_error'])].reset_index(drop=True)
    all_results = all_results.sort_values(by=['Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric']).reset_index(drop=True)
    methods_order = ['pHapCompass', 'HPoP-G', 'WhatsHap']
    all_results = all_results.sort_values(by=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric'],
        key=lambda col: col.map({method: i for i, method in enumerate(methods_order)})).reset_index(drop=True)

    all_results.to_csv(os.path.join(agg_results_path, 'all_methods_results.csv'), index=False)

    # results_df_100 = results_df[results_df['contig'] == '100'].reset_index(drop=True)
    results_df_10 = all_results[all_results['Contig'] == 10].reset_index(drop=True)
    results_df_10.to_csv(os.path.join(agg_results_path, 'all_methods_results_10.csv'), index=False)

    results_df_100 = all_results[all_results['Contig'] == 100].reset_index(drop=True)
    results_df_100.to_csv(os.path.join(agg_results_path, 'all_methods_results_100.csv'), index=False)


def prepare_results():
    agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test'
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    contigs = ['100'] 
    ploidies = ['6']
    coverages = ['10', '50', '100']
    results_dfs = []
    for contig in contigs:
        for ploidy in ploidies:
            for coverage in coverages:
                # results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results')
                results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results_likelihood')
                if os.path.exists(results_path):
                    samples = [f for f in os.listdir(results_path) if 'FFBS' in f]
                    this_result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value', 'length_phased'], index=range(len(samples)*len(metrics)))
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
                        phased_snp = this_results['predicted_haplotypes'].shape[1]
                        for metric in metrics:
                            this_result_df.loc[counter, 'Sample'] = sample_name
                            this_result_df.loc[counter, 'Metric'] = metric
                            this_result_df.loc[counter, 'Value'] = evals[metric]
                            this_result_df.loc[counter, 'length_phased'] = phased_snp
                            counter += 1
                    results_dfs.append(this_result_df)

    results_df = pd.concat(results_dfs, ignore_index=True)
    results_df.to_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_test.csv'), index=False)


    # server3_results = pd.read_csv(os.path.join(agg_results_path, 'sim_awri_results_likelihood10_100_server3.csv'))
    phapcompass_results = pd.read_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_test.csv'))
    whatshapp_results = pd.read_csv(os.path.join(agg_results_path, 'whatshap_results_simulated_data_test.csv'))


    all_results = pd.concat([phapcompass_results, whatshapp_results], ignore_index=True)
    all_results['Contig'] = all_results['Contig'].astype(int)
    all_results['Ploidy'] = all_results['Ploidy'].astype(int)
    all_results['Coverage'] = all_results['Coverage'].astype(int)
    all_results['length_phased'] = all_results['length_phased'].astype(int)

    all_results = all_results.sort_values(by=['Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric']).reset_index(drop=True)
    methods_order = ['pHapCompass', 'WhatsHap' , 'HPoP-G']
    all_results = all_results.sort_values(by=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric'],
        key=lambda col: col.map({method: i for i, method in enumerate(methods_order)})).reset_index(drop=True)

    all_results.to_csv(os.path.join(agg_results_path, 'all_methods_results.csv'), index=False)


def compare_metric_methods(all_results, agg_results_path):
    # compare metrics for different methods:
    for contig in [100]:
        for metric in ['vector_error_rate', 'mismatch_error', 'mec', 'accuracy', 'vector_error']:
            for ploidy in [6]:
                print('Contig:', contig, 'Metric:', metric, 'Ploidy:', ploidy)
                metric_df = all_results[(all_results['Metric'] == metric) & (all_results['Ploidy'] == ploidy)].reset_index(drop=True)
                g = sns.catplot(x="Coverage", y="Value", hue="Method", data=metric_df, kind="box", height=6, aspect=1.5)

                # Add the title and labels
                g.fig.suptitle(f"Contig: {str(contig)}, Metric: {metric.capitalize()}, Ploidy: {str(ploidy).capitalize()}", y=1.05)
                g.set_axis_labels("Coverage", "Value")

                # Move the legend to the top-right, inside the plot area
                g._legend.set_bbox_to_anchor((0.95, 0.9))  # Adjust the position
                g._legend.set_frame_on(True)  # Optional: Add a frame around the legend
                g._legend.set_title("Method")  # Optional: Customize legend title

                # Adjust the layout to ensure everything fits within the figure
                g.fig.subplots_adjust(top=0.85, right=0.9)

                # Save the figure
                # plt.show()
                plt.savefig(os.path.join(agg_results_path, f"compare_{contig}_{metric}_{ploidy}.png"), bbox_inches="tight", dpi=300)
                plt.close()


def plot_contig_length(all_results, agg_results_path):
    for contig in [100]:
        for metric in ['vector_error_rate']:
            for ploidy in [6]:
                # stop
                print('Contig:', contig, 'Ploidy:', ploidy)
                metric_df = all_results[(all_results['Metric'] == metric) & (all_results['Ploidy'] == ploidy) & (all_results['Contig'] == contig)].reset_index(drop=True)

                metric_df['length_phased'] = metric_df['length_phased'].astype(int)
                metric_df['length_phased'] = metric_df['length_phased']/contig * 100

                g = sns.catplot(
                    x="Coverage", y="length_phased", hue="Method",
                    data=metric_df, kind="box", height=6, aspect=1.5
                )

                # Add the title and labels
                g.fig.suptitle("Percentage of phased variants (%)," + f"Contig: {str(contig)}" + f"Ploidy: {str(ploidy)}", y=1.05)
                g.set_axis_labels("Coverage", "phased variants (%)")

                # Move the legend to the top-right, inside the plot area
                g._legend.set_bbox_to_anchor((0.95, 0.9))  # Adjust the position
                g._legend.set_frame_on(True)  # Optional: Add a frame around the legend
                g._legend.set_title("Method")  # Optional: Customize legend title

                # Adjust the layout to ensure everything fits within the figure
                # g.fig.subplots_adjust(top=0.85, right=0.9)
                g.fig.subplots_adjust(top=1, right=0.9)

                # Save the figure
                # plt.show()
                plt.savefig(os.path.join(agg_results_path, f"phased_length_{ploidy}_{contig}.png"), bbox_inches="tight", dpi=300)
                plt.close()


def generate_fragmentfile_test():
    vcf_file_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/AWRI_ploidy6_contig10.vcf'
    modified_vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/modified_AWRI_ploidy6_contig10.vcf'
    bam_file_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/cov_50/bam/00.bam'
    frag_file_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/cov_50/frag/00.frag'

    vcf_file = pysam.VariantFile(vcf_file_path)
    original_snps = [record.pos for record in vcf_file.fetch()] # 1-based positions
    orig_base = [record.ref for record in vcf_file.fetch()]
    alt_base = [record.alts[0] for record in vcf_file.fetch()]

    # Open the BAM file
    read_cnt = 0
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file.fetch():
            # Get the read name
            read_name = read.query_name

            # Get the aligned positions (list of genomic coordinates)
            aligned_positions = read.get_reference_positions(full_length=True)

            # Extract the base at each position of interest
            for pos in original_snps:
                # Check if the position is covered by the read
                if pos in aligned_positions:
                    base_index = aligned_positions.index(pos)
                    base = read.query_sequence[base_index]
                    # print(f"Read: {read_name}, Position: {pos}, Base: {base}")
                    if base == orig_base[original_snps.index(pos)]:
                        # print(0)
                        print(f"Read: {read_name}, Position: {pos}, Base: {base}: 0")
                    elif base == alt_base[original_snps.index(pos)]:
                        # print(1)
                        print(f"Read: {read_name}, Position: {pos}, Base: {base}: 1")
                    else:
                        print(f"Read: {read_name}, Position: {pos}, Base: {base}: -")
            print('------------------------------------------------')
            read_cnt += 1    # else:
                    # print(f"Read: {read_name}, Position: {pos} not covered.")

    fasta_file_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/contig_10_ploidy_6.fa'
    with open(fasta_file_path, 'r') as f:
        haplotypes = f.readlines()
        f.close()

    haplotypes = [h.strip() for h in haplotypes if h.startswith('>') == False]
    for ii, haplotype in enumerate(haplotypes):
        print([haplotype[pos-1] for pos in original_snps])


    readlist = []
    fragment_list = []

    for fragment in open(frag_file_path, 'r'):
        # print(fragment)
        # read index when there are pairs, readname, start pos (in variants), allele sequence, alleles
        parts = fragment.split()
        readlist.append(parts[1])
        positions = []
        alleles = []
        for iii in range(int(parts[0])):
            # process i+2,i+3.... i+4,i+5...
            start_idx_of_read = iii * 2 + 3
            seq_len = len(parts[start_idx_of_read])
            positions.extend(list(range(int(parts[start_idx_of_read - 1]), int(parts[start_idx_of_read - 1]) + seq_len)))
            [alleles.append(int(a)) for a in parts[start_idx_of_read]]
            fragment_list.append(positions)
            fragment_list.append(alleles)


    # List of missing chromosomes to add
    missing_chromosomes = ["haplotype_1", "haplotype_2", "haplotype_3", "haplotype_4", "haplotype_5", "haplotype_6"]

    # Read the original VCF file
    with open(vcf_file_path, "r") as infile:
        lines = infile.readlines()

    # Open the output VCF file for writing
    with open(modified_vcf_path, "w") as outfile:
        header_written = False

        # Track if contigs are already in the VCF
        existing_contigs = {line.split("=")[1].split(",")[0] for line in lines if line.startswith("##contig")}

        for line in lines:
            # Write the original lines
            if not header_written and not line.startswith("#CHROM"):
                outfile.write(line)
                # If the line is a contig line, mark it
                if line.startswith("##contig"):
                    last_contig_line = line
            elif not header_written and line.startswith("#CHROM"):
                # Add missing chromosomes just before the data starts
                for chrom in missing_chromosomes:
                    if chrom not in existing_contigs:
                        outfile.write(f"##contig=<ID={chrom}>\n")
                header_written = True
                outfile.write(line)  # Write the column headers
            else:
                # Write the rest of the data lines
                outfile.write(line)

    print(f"Updated VCF file saved to: {modified_vcf_path}")




    input_vcf = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/modified_AWRI_ploidy6_contig10.vcf"
    output_vcf = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/1_modified_AWRI_ploidy6_contig10.vcf"

    # List of haplotypes to add
    haplotypes = ["haplotype_1", "haplotype_2", "haplotype_3", "haplotype_4", "haplotype_5", "haplotype_6"]

    # Open the input and output files
    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            # Write header lines as-is
            if line.startswith("#"):
                outfile.write(line)
            else:
                # Split the original line into columns
                columns = line.strip().split("\t")
                original_chrom, pos = columns[0], columns[1]
                # Write the original line for AHIQ01000001.1
                outfile.write(line)
                # Add corresponding lines for haplotypes
                for hap in haplotypes:
                    columns[0] = hap  # Replace chromosome name
                    outfile.write("\t".join(columns) + "\n")

    print(f"Expanded VCF file saved to: {output_vcf}")





    # File paths
    bam_file_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/test_out.bam"
    vcf_file_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/1_modified_AWRI_ploidy6_contig10.vcf"

    # Extract contig names from BAM
    bam_header = pysam.AlignmentFile(bam_file_path, "rb").header
    bam_contigs = [ref['SN'] for ref in bam_header['SQ']]

    # Extract contig names from VCF
    vcf_contigs = []
    with open(vcf_file_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith("##contig"):
                contig_name = line.split("ID=")[1].split(",")[0].strip(">")
                vcf_contigs.append(contig_name)
            elif not line.startswith("#"):
                break  # Stop after header

    # Check for contig mismatches
    mismatched_contigs = set(bam_contigs) - set(vcf_contigs)

    print("BAM Contigs:", bam_contigs)
    print("VCF Contigs:", vcf_contigs)
    print("Mismatched Contigs:", mismatched_contigs)





    # Updated BAM file path
    bam_file_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/test_out.bam"

    # Gather summary statistics

    # Open the BAM file
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")

    # Gather summary statistics
    bam_stats = {
        "total_reads": bam_file.mapped + bam_file.unmapped,
        "mapped_reads": bam_file.mapped,
        "unmapped_reads": bam_file.unmapped,
        "references": bam_file.references,  # Contig names
    }

    print("BAM File Statistics:", bam_stats)


    input_vcf = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/1_modified_AWRI_ploidy6_contig10.vcf"
    output_vcf = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/2_modified_AWRI_ploidy6_contig10.vcf" # Output VCF file

    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                # Write header lines as is
                outfile.write(line)
            else:
                # Split the line into columns
                columns = line.strip().split("\t")
                format_field = columns[8]  # FORMAT column
                sample_data = columns[9]  # Sample data
                # stop
                # Modify the GT field if it exists in the FORMAT column
                if "GT" in format_field:
                    gt_index = format_field.split(":").index("GT")
                    sample_data_fields = sample_data.split(":")
                    sample_data_fields[gt_index] = sample_data_fields[gt_index].replace("|", "/")
                    columns[9] = ":".join(sample_data_fields)

                # Write the modified line to the output
                outfile.write("\t".join(columns) + "\n")

    print(f"Fixed VCF file saved as {output_vcf}")


def extract_fragments_0(vcf_file, bam_file, output_file):
    # Parse VCF file and collect variant information
    vcf_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/1_modified_AWRI_ploidy6_contig10.vcf'
    bam_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/test_out.bam'
    output_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/test_out7.frag'
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        bam_contigs = set(bam.references)

    # Parse VCF file and collect variant information
    variants = []
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip header lines
            columns = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, _, format_field, sample_data = columns[:10]
            
            # Skip variants on the excluded chromosome
            if chrom == "AHIQ01000001.1":
                continue

            # Skip variants on contigs not in the BAM file
            if chrom not in bam_contigs:
                continue

            gt_index = format_field.split(":").index("GT")
            genotype = sample_data.split(":")[gt_index]
            if "|" in genotype or "/" in genotype:  # Keep heterozygous positions
                variants.append({
                    "chrom": chrom,
                    "pos": int(pos),  # Keep 1-based position as is
                    "ref": ref,
                    "alt": alt.split(",")
                })

    # Open BAM file and process reads for each variant
    fragments = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for variant in variants:
            chrom, pos, ref, alt_alleles = variant["chrom"], variant["pos"], variant["ref"], variant["alt"]
            
            # Fetch reads overlapping the variant position
            for read in bam.fetch(chrom, pos - 1, pos):  # BAM fetch uses 0-based start, 1-based end
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue  # Skip unmapped or non-primary alignments

                aligned_positions = read.get_reference_positions(full_length=True)
                if pos in aligned_positions:  # Directly use 1-based position for checking
                    # Get the base at the variant position
                    base_index = aligned_positions.index(pos)
                    base = read.query_sequence[base_index]

                    # Determine haplotype informative base (0 for REF, 1 for ALT)
                    if base == ref:
                        haplotype_base = "0"
                    elif base in alt_alleles:
                        haplotype_base = "1"
                    else:
                        haplotype_base = "-"  # Non-informative or mismatched base

                    # Create fragment entry
                    fragment_entry = f"{read.query_name}\t{chrom}-{pos}\t{haplotype_base}\t{base}\n"
                    fragments.append(fragment_entry)

    # Write fragments to output file
    with open(output_file, "w") as output:
        output.writelines(fragments)

    print(f"Fragments saved to {output_file}")


def extract_fragments_old(vcf_file, bam_file, output_file):
    # Parse VCF file and collect variant information into a dictionary
    variants = {}
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip header lines
            columns = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, _, _, _ = columns[:10]
            variants.setdefault(chrom, {})[int(pos)] = {
                "ref": ref,
                "alt": alt.split(","),
            }
    variants = {chrom: variants[chrom] for chrom in variants if chrom != "AHIQ01000001.1"}
    # Open BAM file and process reads
    fragments = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue  # Skip unmapped or non-primary alignments

            read_name = read.query_name
            chrom = bam.get_reference_name(read.reference_id)  # Chromosome for this read
            if chrom not in variants:
                continue  # Skip reads from chromosomes without variants

            aligned_positions = read.get_reference_positions(full_length=True)
            read_sequence = read.query_sequence

            # Process only positions in the current chromosome's variants
            fragment_positions = []
            fragment_bases = []
            for pos in variants[chrom]:
                if pos - 1 in aligned_positions:  # BAM positions are 0-based
                    base_index = aligned_positions.index(pos - 1)
                    base = read_sequence[base_index]
                    if base == variants[chrom][pos]["ref"]:
                        fragment_bases.append("0")
                    elif base in variants[chrom][pos]["alt"]:
                        fragment_bases.append("1")
                    else:
                        fragment_bases.append("-")
                    fragment_positions.append(pos)  # Keep as 1-based for output

            # Write a fragment if the read spans multiple variants
            if len(fragment_positions) > 1:
                fragments.append(
                    f"{read_name}\t{chrom}\t{len(fragment_positions)}\t{''.join(fragment_bases)}\n"
                )


    # Write fragments to output file
    with open(output_file, "w") as output:
        output.writelines(fragments)

    print(f"Fragments saved to {output_file}")



    import pysam


def extract_custom_fragments(vcf_file, bam_file, output_file, distance_threshold=10):
    # Parse VCF file and collect variant information into a dictionary
    # Parse VCF file and collect variant information into a dictionary
    variants = {}
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip header lines
            columns = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, _, _, _ = columns[:10]
            variants.setdefault(chrom, {})[int(pos)] = {
                "ref": ref,
                "alt": alt.split(","),
            }

    # Open BAM file and process reads
    fragments = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue  # Skip unmapped or non-primary alignments

            read_name = read.query_name
            chrom = bam.get_reference_name(read.reference_id)  # Chromosome for this read
            if chrom not in variants:
                continue  # Skip reads from chromosomes without variants

            aligned_positions = read.get_reference_positions(full_length=True)
            read_sequence = read.query_sequence

            # Process only positions in the current chromosome's variants
            variant_positions = []
            haplotype_bases = []
            for pos in variants[chrom]:
                if pos - 1 in aligned_positions:  # BAM positions are 0-based
                    base_index = aligned_positions.index(pos - 1)
                    base = read_sequence[base_index]
                    if base == variants[chrom][pos]["ref"]:
                        haplotype_bases.append("0")
                    elif base in variants[chrom][pos]["alt"]:
                        haplotype_bases.append("1")
                    else:
                        haplotype_bases.append("-")
                    variant_positions.append(pos)  # Keep as 1-based for output

            # Split into parts based on gap_threshold = 1
            parts = []
            current_part = []
            current_bases = []
            for i in range(len(variant_positions)):
                if current_part and variant_positions[i] != current_part[-1] + 1:
                    parts.append((current_part, current_bases))
                    current_part = []
                    current_bases = []
                current_part.append(variant_positions[i])
                current_bases.append(haplotype_bases[i])
            if current_part:
                parts.append((current_part, current_bases))

            # Write the fragment parts
            if parts:
                parts_count = len(parts)
                fragment_line = f"{parts_count}\t{chrom}-{variant_positions[0]}_MP\t"
                for part_positions, part_bases in parts:
                    fragment_line += f"{part_positions[0]}\t{len(part_bases)}\t{''.join(part_bases)}\t"
                fragment_line += "GGCGGGG\n"  # Arbitrary sequence
                fragments.append(fragment_line)

    # Write fragments to output file
    with open(output_file, "w") as output:
        output.writelines(fragments)

    print(f"Fragments saved to {output_file}")


def change_haplotype_orientation():
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test/contig_100/ploidy_6/haplotypes.csv'
    hap_df = pd.read_csv(hap_path).T
    columns = ['haplotype_' + str(i + 1) for i in range(len(hap_df.columns))]
    hap_df.columns = columns
    hap_df.to_csv(hap_path, index=False)



def extract_fragments(vcf_file, bam_file, output_file):
    # Parse VCF file and collect variant information into a dictionary
    vcf_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/AWRI_ploidy6_contig10.vcf'
    bam_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/cov_50/bam/00.bam'
    variants = {}
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip header lines
            columns = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, _, _, _ = columns[:10]
            variants.setdefault(chrom, {})[int(pos)] = {
                "ref": ref,
                "alt": alt.split(","),
            }

    # Open BAM file and process reads
    fragments = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue  # Skip unmapped or non-primary alignments
            
            read_name = read.query_name
            chrom = bam.get_reference_name(read.reference_id)  # Chromosome for this read
            
            if chrom not in variants:
                continue  # Skip reads from chromosomes without variants

            aligned_positions = read.get_reference_positions(full_length=True)
            read_sequence = read.query_sequence

            # Process only positions in the current chromosome's variants
            fragment_positions = []
            fragment_bases = []
            for pos in variants[chrom]:
                if pos - 1 in aligned_positions:  # BAM positions are 0-based
                    base_index = aligned_positions.index(pos - 1)
                    base = read_sequence[base_index]
                    if base == variants[chrom][pos]["ref"]:
                        fragment_bases.append("0")
                    elif base in variants[chrom][pos]["alt"]:
                        fragment_bases.append("1")
                    else:
                        fragment_bases.append("-")
                    fragment_positions.append(pos)  # Keep as 1-based for output
            
            # Write a fragment if the read spans multiple variants
            if len(fragment_positions) > 1:
                # stop
                chrom_variants = sorted(variants[chrom].keys())
                fragment_indices = [chrom_variants.index(pos) + 1 for pos in fragment_positions]

                # Segment into consecutive parts
                parts = []
                current_part_indices = []
                current_part_bases = []

                for idx, base in zip(fragment_indices, fragment_bases):
                    if current_part_indices and idx != current_part_indices[-1] + 1:
                        # Append the current part when a gap is found
                        parts.append((current_part_indices, current_part_bases))
                        current_part_indices = []
                        current_part_bases = []
                    # Add the current index and base to the current part
                    current_part_indices.append(idx)
                    current_part_bases.append(base)
                            # Append the last part
                if current_part_indices:
                    parts.append((current_part_indices, current_part_bases))

            # Construct the fragment line
                fragment_line = f"{len(parts)}\t{chrom}-{read_name}\t"
                for part_indices, part_bases in parts:
                    print(part_indices, part_bases)
                    fragment_line += f"{part_indices[0]}\t{''.join(part_bases)}\t"
                fragment_line += "GG\n"


                # fragments.append(
                #     f"{read_name}\t{chrom}-{'-'.join(map(str, fragment_positions))}\t{len(fragment_positions)}\t{''.join(fragment_bases)}\n"
                # )
                fragments.append(fragment_line)

    # Write fragments to output file
    with open(output_file, "w") as output:
        output.writelines(fragments)

    print(f"Fragments saved to {output_file}")