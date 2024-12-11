import os
import random
import pandas as pd

reference_path = '/labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna'

def extract_contigs():
    sh_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/contigs/extract_contigs.sh'
    dataset = 'SRR942191'
    contigs = ['AHIQ01000' + str(cn).zfill(3) + '.1' for cn in range(1, 325)]
    to_print = ''
    for contig in contigs:
        command = 'samtools view -b {}_sorted.bam {} > {}_{}.bam\n'.format(dataset, contig, dataset, contig)
        to_print += command
    with open(sh_path, 'w') as f:
        f.write(to_print)

    # command = 'samtools view -b SRR942191_sorted.bam AHIQ01000001.1 > SRR942191_AHIQ01000001.1.bam'


def convert_to_vcf():
    sh_path = '/labs/Aguiar/pHapCompass/scripts/convert_2_vcf.sh'
    to_print = 'module load freebayes/1.3.4\n'
    bam_folder = '/labs/Aguiar/pHapCompass/datasets/SRR942191/contigs'
    vcf_folder = '/labs/Aguiar/pHapCompass/datasets/SRR942191/vcf_files'
    reference_path = '/labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna'
    dataset = 'SRR942191'
    contigs = ['AHIQ01000' + str(cn).zfill(3) + '.1' for cn in range(1, 325)]
    for contig in contigs:
        command = 'freebayes -f {} --ploidy 3 {}/{}_{}.bam > {}/{}_{}.vcf\n'.format(reference_path, bam_folder, dataset, contig, vcf_folder, dataset, contig)
        to_print += command

    with open(sh_path, 'w') as f:
        f.write(to_print)


def split_reference_fasta():
    sh_path = '/labs/Aguiar/pHapCompass/scripts/01_split_reference.sh'
    # to_print = ''
    reference_path = '/labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna'
    to_print = 'samtools faidx {}\n\n'.format(reference_path)
    reference_chr_path = '/labs/Aguiar/pHapCompass/references/AWRI1499/contigs'
    if not os.path.exists(reference_chr_path):
        os.makedirs(reference_chr_path)
    # split_command = 'samtools faidx input_reference.fa AHIQ01000001.1 > contig_AHIQ01000001.1.fa'
    contigs = ['AHIQ01000' + str(cn).zfill(3) + '.1' for cn in range(1, 325)]
    for contig in contigs:
        split_command = 'samtools faidx {} {} > {}/contig_{}.fa\n'.format(reference_path, contig, reference_chr_path, contig)
        index_command = 'samtools faidx {}/contig_{}.fa\n\n'.format(reference_chr_path, contig)
        to_print += split_command
        to_print += index_command
    with open(sh_path, 'w') as f:
        f.write(to_print)


def simulating_script():
    chr_reference_path = '/labs/Aguiar/pHapCompass/references/AWRI1499/contigs/contig_AHIQ01000001.1.fa'
    # command1 = 'art_454 -s -t -r 43 input_reference.fa -o output_single -f 5'
    # command2 = 'art_454 -s -t -r 43 input_reference.fa -o output_paired -f 5 -m 800 -s 150'
    to_print = ''

    simulated_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/sim'
    sh_path = '/labs/Aguiar/pHapCompass/scripts/02_simulate.sh'
    # random_seed = 43
    coverages = [5, 10, 15, 20, 25]
    mean_insert_lengths = [800]
    std_insert_lengths = [150]
    for mil in mean_insert_lengths:
        for sil in std_insert_lengths:
            for coverage in coverages:
                this_sim_path = os.path.join(simulated_path, 'cov_{}'.format(coverage))
                if not os.path.exists(this_sim_path):
                    os.makedirs(this_sim_path)
                for rd in range(100):
                    rn = random.randint(1, 2**32)
                    single_cmd = 'art_454 -s -t -r {} {} {}/{}_single {}\n'.format(rn, chr_reference_path, this_sim_path, str(rd).zfill(2), coverage)
                    paired_cmd = 'art_454 -s -t -r {} {} {}/{}_paired {} {} {}\n\n'.format(rn, chr_reference_path, this_sim_path, str(rd).zfill(2), coverage, mil, sil)
                    to_print += single_cmd
                    to_print += paired_cmd

    with open(sh_path, 'w') as f:
        f.write(to_print)


def merge_fastq_files():
    sh_path = '/labs/Aguiar/pHapCompass/scripts/03_merge_fastq.sh'
    to_print = ''
    simulated_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/contig1_k3'
    fastq_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq'
    coverages = [10, 20, 30, 40, 50]
    for coverage in coverages:
        this_sim_path = os.path.join(simulated_path, 'cov_{}'.format(coverage))
        this_fq_path = os.path.join(fastq_path, 'cov_{}'.format(coverage))
        # merged_bam = os.path.join(simulated_path, 'cov_{}'.format(coverage), 'merged')
        if not os.path.exists(this_fq_path):
            os.makedirs(this_fq_path)
        for rd in range(100):
            concat_fastqs = 'cat {}/{}_*.fq > {}/{}.fastq\n'.format(this_sim_path, str(rd).zfill(2), this_fq_path, str(rd).zfill(2))
            to_print += concat_fastqs

    with open(sh_path, 'w') as f:
        f.write(to_print)
            # single_sam = '{}/{}_single.sam'.format(this_sim_path, str(rd).zfill(2))
            # paired_sam = '{}/{}_paired.sam'.format(this_sim_path, str(rd).zfill(2))
            # convert_sam_bam = 'samtools view -b {} > {}/{}_single.bam\n'.format(single_sam, merged_bam, str(rd).zfill(2))


def align_fastq_files():
    reference_path = '/labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna'
    sh_path = '/labs/Aguiar/pHapCompass/scripts/04_align_fastq.sh'
    to_print = 'module load bwa-mem2/2.1\nmodule load bwa/0.7.17\n\n'
    simulated_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/contig1_k3'
    fastq_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq'
    bam_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam'
    coverages = [10, 20, 30, 40, 50]
    for coverage in coverages:
        this_fq_path = os.path.join(fastq_path, 'cov_{}'.format(coverage))
        this_bam_path = os.path.join(bam_path, 'cov_{}'.format(coverage))
        # merged_bam = os.path.join(simulated_path, 'cov_{}'.format(coverage), 'merged')
        if not os.path.exists(this_bam_path):
            os.makedirs(this_bam_path)
        for rd in range(100):
            align_com ='bwa mem {} {}/{}.fastq > {}/{}.sam\n'.format(reference_path, this_fq_path, str(rd).zfill(2), this_bam_path, str(rd).zfill(2))
            sort_com = 'samtools view -Sb {}/{}.sam | samtools sort -o {}/{}.bam\n'.format(this_bam_path, str(rd).zfill(2), this_bam_path, str(rd).zfill(2))
            index_com = 'samtools index {}/{}.bam\n'.format(this_bam_path, str(rd).zfill(2))
            rmv_com = 'rm {}/{}.sam\n\n'.format(this_bam_path, str(rd).zfill(2))
            to_print += align_com
            to_print += sort_com
            to_print += index_com
            to_print += rmv_com

    with open(sh_path, 'w') as f:
        f.write(to_print)

    

    def step1():
        """
        'sed 's/[YRWSMK]/N/g' /labs/Aguiar/pHapCompass/references/AWRI1499/contigs/contig_AHIQ01000001.1.fa > /labs/Aguiar/pHapCompass/references/AWRI1499/contigs_noamb/contig_AHIQ01000001.1.fa'
        'samtools faidx /labs/Aguiar/pHapCompass/references/AWRI1499/contigs_noamb/contig_AHIQ01000001.1.fa'        
        """

        chr_reference_path = '/labs/Aguiar/pHapCompass/references/AWRI1499/contigs_noamb/contig_AHIQ01000001.1.fa'
        chr_vcf_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/vcf_files/SRR942191_AHIQ01000001.1.vcf'
        chr_bam_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/contigs/SRR942191_AHIQ01000001.1.bam'
        new_genome_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/genomes'



        to_print = 'module load bcftools/1.20\nmodule load htslib/1.20\nmodule load freebayes/1.3.4\n\n'
        vcf_cmd = 'freebayes -f {} --ploidy 3 {} > {}\n'.format(chr_reference_path, chr_bam_path, chr_vcf_path)
        zip_com = 'bgzip {}\n'.format(chr_vcf_path)
        index_com1 = 'tabix -p vcf {}.gz\n'.format(chr_vcf_path)
        index_com2 ='bcftools index {}.gz\n\n'.format(chr_vcf_path)
        to_print += vcf_cmd
        to_print += zip_com
        to_print += index_com1
        to_print += index_com2
        # print(to_print)
        # new_genome_file = os.path.join(new_genome_path, 'genome_1.fa')
        if not os.path.exists(new_genome_path):
            os.makedirs(new_genome_path)
        ploidy = 3
        for i in range(ploidy):
            new_genome_file = os.path.join(new_genome_path, 'genome_{}.fa'.format(i+1))
            command = 'bcftools consensus -f {} -o {} {}.gz\n'.format(chr_reference_path, new_genome_file, chr_vcf_path)
            index_com = 'samtools faidx {}\n\n'.format(new_genome_file)
        # 'bcftools consensus -f /labs/Aguiar/pHapCompass/references/AWRI1499/contigs/contig_AHIQ01000001.1_noambig.fa -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/genomes/genome_1.fa /labs/Aguiar/pHapCompass/datasets/SRR942191/vcf_files/SRR942191_AHIQ01000001.1.vcf.gz'
            to_print += command
            to_print += index_com
        print(to_print)


def generate_genomes_fasta(input_vcf_path, contig_fasta, output_fasta_path, output_vcf, ploidy):
    # vcf_in = pysam.VariantFile(input_vcf_path)
    # output_vcf_prefix = 'vcf_out'
    ploidy = 3 
    input_vcf_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/vcf_files/SRR942191_AHIQ01000001.1.vcf'
    contig_fasta = '/labs/Aguiar/pHapCompass/references/AWRI1499/contigs_noamb/contig_AHIQ01000001.1.fa'
    output_fasta_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/genomes/'
    output_vcf = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/genomes/vcf_out.vcf'

    haplotypes = [[] for _ in range(ploidy)]  # Initialize ploidy haplotype sequences

    with open(contig_fasta) as f:
        contig_name = f.readline().strip()[1:]  # Get contig name without '>'
        contig_seq = f.read().replace("\n", "")  # Get sequence without newlines

    # Initialize ploidy genomes with the reference sequence
    genomes = [list(contig_seq) for _ in range(ploidy)]  # List of lists for mutable strings

    # Open the input VCF file for reading and the output VCF for writing manually
    vcf_in = pysam.VariantFile(input_vcf_path)
    with open(output_vcf, "w") as vcf_out:
        # Write the header from the input VCF to the output VCF
        vcf_out.write(str(vcf_in.header))

        # Iterate over each SNP in the VCF
        for record in vcf_in.fetch():
            pos = record.pos - 1  # Convert 1-based VCF position to 0-based for indexing
            ref = record.ref
            alts = record.alts[:1]  # Only use the first alternative allele

            # Generate a random heterozygous genotype for ploidy = 3
            if len(alts) > 0:
                # print(pos)
                
                g = random.randint(1, ploidy - 1)
                selected_genomes = random.sample(range(ploidy), g)
                genotype = tuple(1 if i in selected_genomes else 0 for i in range(ploidy))

                # Assign alleles in the genomes based on the selected genotype
                for i in range(ploidy):
                    genomes[i][pos] = alts[0] if genotype[i] == 1 else ref
                    haplotypes[i].append(genotype[i])

                # Format the genotype as VCF-style (e.g., 0/1/0)
                genotype_str = "/".join(map(str, genotype))
                qual = round(random.uniform(3, 75), 4)

                # Write the record line directly to the output VCF file
                vcf_out.write(f"{record.chrom}\t{record.pos}\t.\t{ref}\t{alts[0]}\t{qual}\tPASS\t*\tGT:GQ\t{genotype_str}:100\n")


    # Write each genome to a separate FASTA file
    for i in range(ploidy):
        # output_fasta = f"{contig_name}_genome_{i + 1}.fa"
        output_fasta = os.path.join(output_fasta_path, 'genome_{}.fa'.format(i+1))
        with open(output_fasta, "w") as f_out:
            f_out.write(f">{contig_name}_genome_{i + 1}\n")
            f_out.write("".join(genomes[i]) + "\n")

    # write haplotypes to a DataFrame
    haplotype_df = pd.DataFrame(haplotypes).T  # Transpose so each column represents a haplotype

    # Rename columns for clarity
    haplotype_df.columns = [f"haplotype_{i + 1}" for i in range(len(haplotypes))]

    # Save the DataFrame as a CSV if needed
    haplotype_df.to_csv(os.path.join(output_fasta_path, 'haplotypes.csv'), index=False)


# def simulating_script():
#     chr_reference_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/genomes/genome_1.fa'
#     # command1 = 'art_454 -s -t -r 43 input_reference.fa -o output_single -f 5'
#     # command2 = 'art_454 -s -t -r 43 input_reference.fa -o output_paired -f 5 -m 800 -s 150'
#     to_print = ''

#     simulated_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/contig1_k3'
#     if not os.path.exists(simulated_path):
#         os.makedirs(simulated_path)
#     sh_path = '/labs/Aguiar/pHapCompass/scripts/02_simulate.sh'
#     # random_seed = 43
#     coverages = [5, 10, 15, 20, 25]
#     coverages = [10, 20, 30, 40, 50]
#     mean_insert_lengths = [800]
#     std_insert_lengths = [150]
#     for mil in mean_insert_lengths:
#         for sil in std_insert_lengths:
#             for coverage in coverages:
#                 this_sim_path = os.path.join(simulated_path, 'cov_{}'.format(coverage))
#                 if not os.path.exists(this_sim_path):
#                     os.makedirs(this_sim_path)
#                 for rd in range(100):
#                     rn = random.randint(1, 2**32)
#                     for i in range(ploidy):
#                         this_chr_reference_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/genomes/genome_{}.fa'.format(i+1)
#                         coverage = ? 
#                         single_cmd = './art_454 -s -t -r {} {} {}/{}_single_{}_{} {}\n'.format(rn, this_chr_reference_path, this_sim_path, str(rd).zfill(2), i+1, coverage)
#                         paired_cmd = './art_454 -s -t -r {} {} {}/{}_paired_{}_{} {} {} {}\n\n'.format(rn, this_chr_reference_path, this_sim_path, str(rd).zfill(2), i+1, coverage, mil, sil)
#                         to_print += single_cmd
#                         to_print += paired_cmd
#                     # single_cmd = './art_454 -s -t -r {} {} {}/{}_single {}\n'.format(rn, chr_reference_path, this_sim_path, str(rd).zfill(2), coverage)
#                     # paired_cmd = './art_454 -s -t -r {} {} {}/{}_paired {} {} {}\n\n'.format(rn, chr_reference_path, this_sim_path, str(rd).zfill(2), coverage, mil, sil)
#                     # to_print += single_cmd
#                     # to_print += paired_cmd

#     with open(sh_path, 'w') as f:
#         f.write(to_print)


def simulate_for_coverage():
    to_print = ''
    simulated_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/contig1_k3'
    if not os.path.exists(simulated_path):
        os.makedirs(simulated_path)
    sh_path = '/labs/Aguiar/pHapCompass/scripts/02_simulate.sh'

    coverages = [10, 20, 30, 40, 50]
    mean_insert_lengths = [800]
    std_insert_lengths = [150]

    ploidy = 3

    for mil in mean_insert_lengths:
        for sil in std_insert_lengths:
            for target_coverage in coverages:
                # Set up each coverage according to the target and constraints
                a = target_coverage / 6  # Initial approximation
                b = a / 2
                c = target_coverage / 6
                d = c / 2
                e = target_coverage / 3 - c  # Remaining coverage to balance
                f = e / 2

                # Adjusting to ensure total sum is exactly the target
                # Calculate current sum and scale to match target
                current_sum = a + b + c + d + e + f
                scale_factor = target_coverage / current_sum
                a, b, c, d, e, f = [int(round(x * scale_factor)) for x in [a, b, c, d, e, f]]

                # Check the sum for validation
                print(a + b + c + d + e + f, target_coverage)
                # assert a + b + c + d + e + f == target_coverage, "Coverages do not sum to target"

                # Create output directories for this coverage
                this_sim_path = os.path.join(simulated_path, f'cov_{target_coverage}')
                if not os.path.exists(this_sim_path):
                    os.makedirs(this_sim_path)

                # Loop to simulate reads for each genome
                for rd in range(100):
                    rn = random.randint(1, 2**32)
                    for i in range(ploidy):
                        this_chr_reference_path = f'/labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/genomes/genome_{i + 1}.fa'
                        
                        # Assign coverages to single and paired commands
                        single_coverage = [a, c, e][i]
                        paired_coverage = [b, d, f][i]

                        # Construct command strings
                        single_cmd = f'./art_454 -s -t -r {rn} {this_chr_reference_path} {this_sim_path}/{str(rd).zfill(2)}_single_{i + 1} {single_coverage}\n'
                        paired_cmd = f'./art_454 -s -t -r {rn} {this_chr_reference_path} {this_sim_path}/{str(rd).zfill(2)}_paired_{i + 1} {paired_coverage} {mil} {sil}\n\n'
                        
                        # Add commands to output script
                        to_print += single_cmd
                        to_print += paired_cmd

    # Save commands to shell script
    with open(sh_path, "w") as sh_file:
        sh_file.write(to_print)



def inspect_snps_positions_dist():
    # vcf_path = '/labs/Aguiar/pHapCompass/datasets/SRR942191/vcf_files/SRR942191_AHIQ01000001.1.vcf'
    vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/SRR942191/vcf_files/SRR942191_AHIQ01000001.1.vcf'
    vcf_in = pysam.VariantFile(vcf_path)
    positions = [record.pos for record in vcf_in.fetch()]
    differences = [positions[i+1] - positions[i] for i in range(len(positions) - 1)]

    # plt.hist(differences, bins=6, edgecolor='black')
    # # Add labels and title
    # plt.xlabel('Value')
    # plt.ylabel('Frequency')
    # plt.title('Simple Histogram')
    # # Show the plot
    # plt.show()


    
