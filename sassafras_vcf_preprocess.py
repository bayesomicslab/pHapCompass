
##### IMPORT LIBRARIES
import io
import os
import pandas as pd
import re
import gzip
import pysam

######################################################################################################
####################################### VCF READER ###################################################
######################################################################################################

## I. VCF TO PANDAS READER FUNCTION

# Adapted from source code publicly available from Kevin Blighe (https://www.biostars.org/p/302940/).
# Defines the read_vcf function.
def read_vcf(path):
    if path.endswith('.gz'):
        with gzip.open(path, 'rt') as file:
            lines = [line for line in file if not line.startswith('##')]
    else:
        with open(path, 'r') as file:
            lines = [line for line in file if not line.startswith('##')]
    '''
    # Opens the file in read mode.
    with open(path, 'r') as file:
        # Omits header lines (VCF headers start with ##).
        lines = [line for line in file if not line.startswith('##')]
        #INFO_col = [column for column in file if column.startswith('INFO')].split(';')
    '''
    # Returns the VCF file as a pandas dataframe.
    vcf = pd.read_csv(
        #imports the text and joins the standard VCF category columns together separated by tab.
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    # Renames the #CHROM column to CHROM.
    ).rename(columns={'#CHROM': 'CHROM'})
    # Returns the VCF variable to the global space.
    return vcf


# Function to extract the required fields from a VCF file and save to CSV
    # Open the VCF file using pysam
    vcf = pysam.VariantFile(vcf_file)

    # List to hold the extracted data
    extracted_data = []

    # Iterate through each record in the VCF file
    for record in vcf:
        # Extract CHROM, POS, and GT (from FORMAT column)
        chrom = record.chrom
        pos = record.pos

        # Extract GT for the sample (assuming a single sample in the VCF file)
        sample = list(record.samples.values())[0]  # Get the sample data
        gt = sample.get("GT")  # Extract the GT field
        
        # Convert GT tuple to a string representation
        gt_str = "/".join(map(str, gt)) if gt else "./."

        # Append the extracted information to the list
        extracted_data.append([chrom, pos, gt_str])

    # Create a DataFrame from the extracted data
    df = pd.DataFrame(extracted_data, columns=["CHROM", "POS", "GT"])

    # Save the DataFrame to a CSV file
    df.to_csv(output_csv, index=False)


# Function to extract the required fields from a VCF file and save to CSV
def extract_vcf_to_csv(vcf_file, output_csv):
    # Open the VCF file using pysam
    vcf = pysam.VariantFile(vcf_file)

    # List to hold the extracted data
    extracted_data = []

    # Iterate through each record in the VCF file
    for record in vcf:
        # Extract CHROM, POS, REF, ALT, and GT (from FORMAT column)
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alt = ",".join(record.alts) if record.alts else ""

        # Extract GT for the sample (assuming a single sample in the VCF file)
        sample = list(record.samples.values())[0]  # Get the sample data
        gt = sample.get("GT")  # Extract the GT field
        
        # Convert GT tuple to a string representation
        gt_str = "/".join(map(str, gt)) if gt else "./."

        # Append the extracted information to the list
        extracted_data.append([chrom, pos, ref, alt, gt_str])

    # Create a DataFrame from the extracted data
    df = pd.DataFrame(extracted_data, columns=["CHROM", "POS", "REF", "ALT", "GT"])

    # Save the DataFrame to a CSV file
    df.to_csv(output_csv, index=False)


def heterozygous_filter():
    folder_path = '/labs/Aguiar/pHapCompass/datasets/sassafras/vcf'
    vcfs = [f for f in os.listdir(folder_path) if f.endswith('.vcf') and 'filter' not in f]
    for f in vcfs:
        input_vcf  = os.path.join(folder_path, f)
        output_vcf  = os.path.join(folder_path, 'filtered_' + f)

                # Open the input VCF file
        with pysam.VariantFile(input_vcf, "r") as vcf_in:
            # Create an output VCF file with the same header as the input
            with pysam.VariantFile(output_vcf, "w", header=vcf_in.header) as vcf_out:
                for record in vcf_in:
                    # Access the genotype information
                    genotype = record.samples[0]["GT"]
                    
                    # Check if genotype is not (0, 0, 0, 0) or (1, 1, 1, 1)
                    if genotype not in [(0, 0, 0, 0), (1, 1, 1, 1)]:
                        vcf_out.write(record)

def get_genotype():
    folder_path = '/labs/Aguiar/pHapCompass/datasets/sassafras/vcf'
    vcfs = [f for f in os.listdir(folder_path) if f.endswith('.vcf') and 'filtered' in f]
    for f in vcfs:
        vcf_file = os.path.join(folder_path, f)
        output_csv = os.path.join(folder_path, f.replace('.vcf', '.csv'))
        extract_vcf_to_csv(vcf_file, output_csv)
        print(f"Data extracted and saved to {output_csv}")
        df = pd.read_csv(output_csv)
        df = df.sort_values(by=['CHROM', 'POS'])
        df["genotype"] = df["GT"].apply(lambda gt: sum(map(int, gt.split("/"))))

        gt_split = df["GT"].str.split("/", expand=True)
        gt_split.columns = ["haplotype_1", "haplotype_2", "haplotype_3", "haplotype_4"]

        # Convert the columns to integers (optional, based on your use case)
        gt_split = gt_split.astype(int)
        haplotype_name = os.path.join(folder_path, 'haplotype_' + f.replace('.vcf', '.csv'))
        gt_split.to_csv(haplotype_name, index=False)
        


if __name__ == "__main__":


    # vcf_file = '/core/globus/cgi/sassafras_UT-Knoxville/VCF/sassafras_albidum_hap1_combined_biallelic_snps-only.vcf'
    # output_csv = '/labs/Aguiar/pHapCompass/datasets/sassafras/sassafras_albidum_hap1_combined_biallelic_snps_only.csv'
    # # vcf = read_vcf(vcf_file)
    # # vcf.to_csv(output_csv, index=False)

    # extract_vcf_to_csv(vcf_file, output_csv)

    # print(f"Data extracted and saved to {output_csv}")
    # df = pd.read_csv(output_csv)
    # chromosomes = ' '.join(sorted(list(df['CHROM'].unique())))

    # orig_vcf_path = '/core/globus/cgi/sassafras_UT-Knoxville/VCF/sassafras_albidum_hap1_combined_biallelic_snps-only.vcf'
    # orig_csv_path = '/labs/Aguiar/pHapCompass/datasets/sassafras/vcf/sassafras_albidum_hap1_combined_biallelic_snps_only.csv'
    # extract_vcf_to_csv(orig_vcf_path, orig_csv_path)



    # folder_path = '/labs/Aguiar/pHapCompass/datasets/sassafras/vcf'
    # vcfs = [f for f in os.listdir(folder_path)]
    # for f in vcfs:
    #     vcf_file = os.path.join(folder_path, f)
    #     output_csv = os.path.join(folder_path, f.replace('.vcf', '.csv'))
    #     extract_vcf_to_csv(vcf_file, output_csv)
    #     print(f"Data extracted and saved to {output_csv}")
    # csvs = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
    # for f in csvs:
    #     csv_file = os.path.join(folder_path, f)
    #     df = pd.read_csv(csv_file)
    #     df = df.sort_values(by=['CHROM', 'POS'])
    #     df["genotype"] = df["GT"].apply(lambda gt: sum(map(int, gt.split("/"))))
    #     # df["genotype"].unique()
    #     df = df[(df['genotype'] != 0)&(df['genotype'] != 4)].reset_index(drop=True)
    #     # Add a new column 'SNP_id' starting from 1 up to the number of rows
    #     df["SNP_id"] = range(1, len(df) + 1)
    #     df['genotype'] = df['genotype'].astype(str)
    #     genotype = ''.join(list(df['genotype']))
    #     # df = df.drop(columns=['CHROM', 'REF', 'ALT'])
    #     df.to_csv(csv_file, index=False)
    #     output_file = os.path.join(folder_path, f.replace('.csv', '.txt'))
    #     with open(output_file, "w") as file:
    #         file.write(genotype)

    # heterozygous_filter()



# # Function to extract the required fields from a VCF file and save to CSV
# # Example usage
# vcf_file = "example.vcf"  # Replace with your VCF file path
# output_csv = "output.csv"  # Replace with desired output CSV file path
# extract_vcf_to_csv(vcf_file, output_csv)

# print(f"Data extracted and saved to {output_csv}")
