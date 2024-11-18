import os
import pandas as pd
import random
import pysam



def extract_column_NA12878(file_path):
    # file_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.1/windows/50000/sample_NA12878.txt'
    df = pd.read_csv(file_path, skiprows=48, sep=' ')
    hetero_df = df[df['NA12878'].isin(['1|0', '0|1'])].reset_index(drop=True)
    hetero_df = hetero_df.sort_values(by='POS').reset_index(drop=True)
    hetero_df['sim_pos'] = hetero_df.index
    return hetero_df


def extract_column_NA12878_v2(file_path):


    # file_path = "/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/windows/10/hapref_chr21_filtered.vcf.bgz_sample0_len10.bcf"
    
    file_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/hapref_chr21_filtered.vcf.bgz'
    bcf_file = pysam.VariantFile(file_path)


    output_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.vcf.gz'
    output_vcf = pysam.VariantFile(output_file, 'w', header=bcf_file.header)

    # Specify the sample you want to extract
    target_sample = "NA12878"

    # Create an empty list to store the records
    records = []

    # Iterate through the records in the BCF file
    for record in bcf_file.fetch():
        if target_sample in record.samples:
            # Extract genotype information
            genotype = record.samples[target_sample]["GT"]
            # Add record to DataFrame if it's relevant
            records.append({
                "CHROM": record.chrom,
                "POS": record.pos,
                "REF": record.ref,
                "ALT": ",".join(record.alts) if record.alts else None,
                "GENOTYPE": genotype
            })
        # Write the record to the new VCF file
        output_vcf.write(record)

    # Close the output VCF file
    output_vcf.close()
    

    # Convert the list of records into a DataFrame
    df = pd.DataFrame(records)

    # Filter for heterozygous genotypes (1|0 or 0|1)
    df["GENOTYPE"] = df["GENOTYPE"].apply(lambda x: "|".join(map(str, x)))  # Convert tuples to string
    # hetero_df = df[df["GENOTYPE"].isin(["1|0", "0|1"])].reset_index(drop=True)

    # # Sort by position and add simulated positions
    # hetero_df = hetero_df.sort_values(by="POS").reset_index(drop=True)
    # hetero_df["sim_pos"] = hetero_df.index
    # hetero_df.to_csv('/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv', index=False)
    df = df.sort_values(by="POS").reset_index(drop=True)
    # df['difference'] = df['POS'].diff()
    df.to_csv('/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv', index=False)


def check_compatibility_vcf_fasta():
    # vcf_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/_EGAZ00001239288_HRC.r1-1.EGA.GRCh37.chr21.haplotypes.vcf.gz'
    vcf_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/hapref_chr21_filtered.vcf.gz'

    vcf_df = pd.read_csv(vcf_path, skiprows=48, nrows=10, sep='\t')

    hetero_df = extract_column_NA12878()
    filtered_positions = hetero_df['POS'].values

    # Path to the FASTA file
    # chr21_fasta_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/chr21.fna'
    chr21_fasta_path = '/mnt/research/aguiarlab/data/hg19/chr21.fa'
    chr21_legend_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/_EGAZ00001239288_HRC.r1-1.EGA.GRCh37.chr21.legend.gz'
    # Read the FASTA file, skipping the first line (chromosome name line)
    with open(chr21_fasta_path, 'r') as f:
        fasta_file = f.readlines()

    # Combine the sequence lines into a single string, skipping the first line
    fasta_sequence = "".join(line.strip() for line in fasta_file if not line.startswith(">"))

    legend_df = pd.read_csv(chr21_legend_path, sep=' ', compression='gzip')

    for i in range(10):
        pos = legend_df.loc[i, 'position']
        ref_pos = legend_df.loc[i, 'a0']
        # alt_pos = legend_df.loc[i, 'a1']    
        print(pos, ref_pos, fasta_sequence[pos - 1], legend_df.loc[i, 'a1'])
        # print('----------------------------------')

    legend_df = legend_df[legend_df['position'].isin(list(filtered_positions))].reset_index(drop=True)
    for i in range(10):
        pos = legend_df.loc[i, 'position']
        ref_pos = legend_df.loc[i, 'a0']
        # alt_pos = legend_df.loc[i, 'a1']
        print(pos, ref_pos, fasta_sequence[pos - 1], legend_df.loc[i, 'a1'])


def generate_genomes_fasta():
    # vcf_in = pysam.VariantFile(input_vcf_path)
    # output_vcf_prefix = 'vcf_out'
    # ploidy = 3 
    # input_vcf_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/_EGAZ00001239288_HRC.r1-1.EGA.GRCh37.chr21.haplotypes.vcf.gz'
    # snp_df_NA12878_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv'
    snp_df_NA12878_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.csv'
    snp_df_NA12878 = pd.read_csv(snp_df_NA12878_path)
    snp_positions = list(snp_df_NA12878['POS'].values)
    # input_vcf_path = '/mnt/research/aguiarlab/data/haprefconsort/hap_ref_consort/corephase_data/maf0.01/hapref_chr21_filtered.vcf.bgz'
    input_vcf_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW/maf0.01_hapref_chr21_filtered.vcf.bgz'
    # contig_fasta = '/mnt/research/aguiarlab/data/hg19/chr21.fa'
    contig_fasta = '/labs/Aguiar/pHapCompass/references/EGA.GRCh37/chr21.fa'
    # output_fasta_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW'
    output_fasta_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW'
    is_phased = True
    # output_vcf = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW/maf0.01_hapref_chr21_filtered_NA12878.vcf.gz'
    with open(contig_fasta) as f:
        contig_name = f.readline().strip()[1:]  # Get contig name without '>'
        contig_seq = f.read().replace("\n", "")  # Get sequence without newlines

    contig_lens = [100] # [10, 100, 300, 1000] 
    ploidies = [3] #, 4, 6, 8, 10]
    for contig_len in contig_lens:
        contig_path = os.path.join(output_fasta_path, 'contig_{}'.format(contig_len))
        if not os.path.exists(contig_path):
            os.makedirs(contig_path)
        for ploidy in ploidies:
            contig_ploidy_path = os.path.join(contig_path, 'ploidy_{}'.format(ploidy))
            if not os.path.exists(contig_ploidy_path):
                os.makedirs(contig_ploidy_path)
            snp_counter = 0
            last_snp = snp_positions[0:contig_len][-1]
            output_vcf = os.path.join(contig_ploidy_path, 'maf0.01_hapref_chr21_filtered_NA12878_ploidy{}_contig{}.vcf.gz'.format(ploidy, contig_len))
            haplotypes = [[] for _ in range(ploidy)]  # Initialize ploidy haplotype sequences

            # Initialize ploidy genomes with the reference sequence
            genomes = [list(contig_seq[:last_snp]) for _ in range(ploidy)]  # List of lists for mutable strings

            vcf_in = pysam.VariantFile(input_vcf_path)
            # Create a new header for the output VCF file with only NA12878
            new_header = pysam.VariantHeader()
            for line in vcf_in.header.records:
                new_header.add_line(str(line))  # Copy metadata lines
            new_header.add_sample("NA12878")  # Add only NA12878

            # Open the output VCF file with the new header
            vcf_out = pysam.VariantFile(output_vcf, "wz", header=new_header)

            for record in vcf_in.fetch():
                # print(record)
                
                pos = record.pos - 1  # Convert 1-based VCF position to 0-based for indexing
                ref = record.ref
                alts = record.alts[:1]  # Only use the first alternative allele

                # Generate a random heterozygous genotype for ploidy = 3
                if len(alts) > 0 and pos + 1 in snp_positions and snp_counter < contig_len:
                    # print(pos)
                    snp_counter += 1
                    print(snp_counter)
                    # stop
                    g = random.randint(1, ploidy - 1)
                    selected_genomes = random.sample(range(ploidy), g)
                    genotype = tuple(1 if i in selected_genomes else 0 for i in range(ploidy))

                    # Assign alleles in the genomes based on the selected genotype
                    for i in range(ploidy):
                        genomes[i][pos] = alts[0] if genotype[i] == 1 else ref
                        haplotypes[i].append(genotype[i])

                    # Format the genotype as VCF-style (e.g., 0|1|0)
                    # genotype_str = "|".join(map(str, genotype))
                    qual = 75 # round(random.uniform(3, 75), 4)
                    
                    for k in range(ploidy):
                        # Create a new VCF record with only NA12878
                        new_record = vcf_out.new_record(
                            # contig=record.chrom,
                            contig=f'haplotype_{k + 1}',
                            start=record.start,
                            stop=record.stop,
                            alleles=(ref, alts[0]),
                            id=record.id,
                            qual=qual,
                            filter=record.filter.keys(),
                            info=record.info
                        )
                        new_record.samples["NA12878"]["GT"] = genotype
                        new_record.samples["NA12878"].phased = is_phased
                        vcf_out.write(new_record)

                    # Write the record line directly to the output VCF file
                    # vcf_out.write(f"{record.chrom}\t{record.pos}\t.\t{ref}\t{alts[0]}\t{qual}\tPASS\t*\tGT:GQ\t{genotype_str}:100\n")

            vcf_out.close()

            output_fasta = os.path.join(contig_ploidy_path, f'contig_{contig_len}_ploidy_{ploidy}.fa')
            with open(output_fasta, "w") as f_out:
                for i in range(ploidy):
                    f_out.write(f">haplotype_{i + 1}\n")
                    f_out.write("".join(genomes[i]) + "\n")

            # write haplotypes to a DataFrame
            haplotype_df = pd.DataFrame(haplotypes).T  # Transpose so each column represents a haplotype

            # Rename columns for clarity
            haplotype_df.columns = [f"haplotype_{i + 1}" for i in range(len(haplotypes))]

            # Save the DataFrame as a CSV if needed
            haplotype_df.to_csv(os.path.join(contig_ploidy_path, 'haplotypes.csv'), index=False)


def simulate_fastq_art():
    # main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NEW'
    main_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW'
    # art_path = '/home/mah19006/downloads/art_bin_MountRainier/art_illumina'
    art_path = '/labs/Aguiar/pHapCompass/ART/art_bin_MountRainier/art_illumina'
    # sh_path = os.path.join('/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation', '01_simulate_illumina.sh')
    sh_path = os.path.join('/labs/Aguiar/pHapCompass/scripts/simulation', '01_simulate_illumina.sh')
    
    # to_print = ''
    contig_lens = [10, 100, 300, 1000]
    ploidies = [3, 4, 6, 8, 10]
    coverages = [10, 20, 30, 40, 50]
    mil = 800
    sil = 150
    read_length = 150
    for contig_len in contig_lens:
        for ploidy in ploidies:
            fasta_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'contig_{}_ploidy_{}.fa'.format(contig_len, ploidy))
            this_sh_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '01_simulate_{}_{}.sh'.format(contig_len, ploidy))
            to_print = ''
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                if not os.path.exists(this_cov_path):
                    os.makedirs(this_cov_path)
                fastq_path = os.path.join(this_cov_path, 'fastq')
                if not os.path.exists(fastq_path):
                    os.makedirs(fastq_path)
                for rd in range(100):
                    command = '{} -ss HS25 -i {} -p -l {} -f {} -m {} -s {} -o {}/{}\n'.format(art_path, fasta_path, read_length, coverage, mil, sil, fastq_path, str(rd).zfill(2))
                    to_print += command
                    # single_cmd = './art_454 -s -t -r {} {} {}/{}_single {}\n'.format(rn, chr_reference_path, this_sim_path, str(rd).zfill(2), coverage)
                    # paired_cmd = './art_454 -s -t -r {} {} {}/{}_paired {} {} {}\n\n'.format(rn, chr_reference_path, this_sim_path, str(rd).zfill(2), coverage, mil, sil)
                    # to_print += single_cmd
                    # to_print += paired_cmd
            with open(this_sh_path, 'w') as f:
                f.write(to_print)


def align_fastq_files():
    main_path = '/labs/Aguiar/pHapCompass/simulated_data_NEW'
    chr21_fasta_path = '/labs/Aguiar/pHapCompass/references/EGA.GRCh37/chr21.fa'
    sh_path = os.path.join('/labs/Aguiar/pHapCompass/scripts/simulation', '02_align.sh')
    # to_print = 'module load bwa-mem2/2.1\nmodule load bwa/0.7.17\n\n'
    contig_lens = [10, 100, 300, 1000]
    ploidies =  [3, 4, 6, 8, 10]
    coverages = [10, 20, 30, 40, 50]
    for contig_len in contig_lens:
        for ploidy in ploidies:
            fasta_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'contig_{}_ploidy_{}.fa'.format(contig_len, ploidy))
            fa_index = 'bwa index {}\n\n'.format(fasta_path)
            this_sh_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), '02_align_{}_{}.sh'.format(contig_len, ploidy))
            # to_print += fa_index
            to_print = 'module load bwa-mem2/2.1\nmodule load bwa/0.7.17\n' + fa_index
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                fastq_path = os.path.join(this_cov_path, 'fastq')
                bam_path = os.path.join(this_cov_path, 'bam')
                if not os.path.exists(bam_path):
                    os.makedirs(bam_path)
                for rd in range(100):
                    command = 'bwa mem {} {}/{}1.fq {}/{}2.fq > {}/{}.sam\n'.format(fasta_path, fastq_path, str(rd).zfill(2), fastq_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
                    sort_com = 'samtools view -Sb {}/{}.sam | samtools sort -o {}/{}.bam\n'.format(bam_path, str(rd).zfill(2), bam_path, str(rd).zfill(2))
                    index_com = 'samtools index {}/{}.bam\n'.format(bam_path, str(rd).zfill(2))
                    rmv_com = 'rm {}/{}.sam\n\n'.format(bam_path, str(rd).zfill(2))
                    to_print += command
                    to_print += sort_com
                    to_print += index_com
                    to_print += rmv_com
            with open(this_sh_path, 'w') as f:
                f.write(to_print)

        '/home/FCAM/mhosseini/HaplOrbit/extract_poly/build/extractHAIRS --bam /labs/Aguiar/pHapCompass/simulated_data_NEW/contig_100/ploidy_3/cov_10/bam/00.bam --vcf /labs/Aguiar/pHapCompass/simulated_data_NEW/contig_100/ploidy_3/merged_haplotypes.vcf --out /labs/Aguiar/pHapCompass/simulated_data_NEW/contig_100/ploidy_3/cov_10/frag/00.frag'



if __name__ == '__main__':
    generate_genomes_fasta()