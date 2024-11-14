module load bwa-mem2/2.1
module load bwa/0.7.17

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/00.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/00.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/00.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/00.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/00.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/00.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/01.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/01.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/01.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/01.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/01.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/01.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/02.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/02.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/02.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/02.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/02.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/02.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/03.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/03.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/03.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/03.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/03.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/03.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/04.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/04.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/04.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/04.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/04.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/04.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/05.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/05.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/05.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/05.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/05.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/05.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/06.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/06.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/06.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/06.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/06.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/06.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/07.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/07.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/07.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/07.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/07.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/07.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/08.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/08.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/08.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/08.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/08.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/08.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/09.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/09.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/09.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/09.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/09.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/09.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/10.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/10.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/10.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/10.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/10.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/10.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/11.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/11.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/11.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/11.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/11.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/11.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/12.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/12.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/12.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/12.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/12.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/12.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/13.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/13.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/13.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/13.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/13.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/13.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/14.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/14.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/14.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/14.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/14.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/14.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/15.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/15.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/15.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/15.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/15.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/15.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/16.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/16.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/16.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/16.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/16.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/16.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/17.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/17.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/17.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/17.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/17.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/17.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/18.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/18.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/18.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/18.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/18.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/18.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/19.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/19.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/19.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/19.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/19.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/19.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/20.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/20.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/20.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/20.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/20.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/20.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/21.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/21.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/21.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/21.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/21.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/21.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/22.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/22.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/22.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/22.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/22.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/22.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/23.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/23.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/23.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/23.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/23.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/23.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/24.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/24.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/24.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/24.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/24.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/24.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/25.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/25.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/25.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/25.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/25.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/25.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/26.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/26.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/26.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/26.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/26.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/26.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/27.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/27.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/27.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/27.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/27.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/27.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/28.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/28.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/28.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/28.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/28.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/28.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/29.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/29.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/29.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/29.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/29.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/29.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/30.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/30.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/30.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/30.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/30.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/30.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/31.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/31.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/31.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/31.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/31.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/31.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/32.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/32.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/32.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/32.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/32.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/32.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/33.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/33.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/33.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/33.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/33.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/33.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/34.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/34.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/34.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/34.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/34.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/34.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/35.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/35.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/35.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/35.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/35.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/35.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/36.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/36.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/36.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/36.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/36.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/36.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/37.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/37.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/37.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/37.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/37.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/37.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/38.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/38.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/38.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/38.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/38.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/38.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/39.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/39.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/39.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/39.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/39.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/39.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/40.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/40.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/40.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/40.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/40.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/40.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/41.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/41.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/41.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/41.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/41.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/41.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/42.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/42.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/42.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/42.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/42.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/42.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/43.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/43.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/43.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/43.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/43.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/43.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/44.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/44.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/44.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/44.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/44.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/44.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/45.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/45.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/45.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/45.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/45.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/45.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/46.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/46.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/46.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/46.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/46.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/46.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/47.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/47.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/47.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/47.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/47.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/47.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/48.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/48.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/48.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/48.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/48.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/48.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/49.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/49.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/49.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/49.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/49.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/49.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/50.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/50.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/50.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/50.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/50.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/50.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/51.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/51.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/51.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/51.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/51.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/51.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/52.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/52.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/52.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/52.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/52.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/52.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/53.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/53.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/53.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/53.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/53.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/53.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/54.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/54.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/54.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/54.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/54.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/54.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/55.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/55.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/55.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/55.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/55.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/55.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/56.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/56.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/56.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/56.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/56.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/56.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/57.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/57.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/57.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/57.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/57.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/57.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/58.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/58.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/58.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/58.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/58.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/58.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/59.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/59.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/59.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/59.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/59.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/59.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/60.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/60.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/60.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/60.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/60.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/60.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/61.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/61.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/61.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/61.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/61.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/61.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/62.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/62.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/62.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/62.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/62.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/62.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/63.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/63.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/63.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/63.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/63.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/63.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/64.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/64.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/64.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/64.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/64.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/64.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/65.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/65.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/65.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/65.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/65.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/65.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/66.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/66.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/66.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/66.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/66.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/66.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/67.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/67.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/67.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/67.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/67.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/67.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/68.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/68.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/68.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/68.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/68.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/68.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/69.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/69.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/69.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/69.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/69.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/69.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/70.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/70.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/70.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/70.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/70.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/70.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/71.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/71.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/71.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/71.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/71.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/71.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/72.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/72.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/72.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/72.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/72.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/72.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/73.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/73.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/73.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/73.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/73.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/73.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/74.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/74.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/74.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/74.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/74.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/74.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/75.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/75.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/75.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/75.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/75.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/75.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/76.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/76.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/76.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/76.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/76.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/76.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/77.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/77.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/77.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/77.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/77.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/77.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/78.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/78.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/78.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/78.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/78.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/78.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/79.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/79.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/79.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/79.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/79.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/79.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/80.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/80.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/80.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/80.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/80.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/80.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/81.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/81.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/81.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/81.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/81.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/81.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/82.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/82.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/82.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/82.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/82.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/82.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/83.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/83.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/83.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/83.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/83.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/83.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/84.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/84.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/84.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/84.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/84.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/84.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/85.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/85.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/85.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/85.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/85.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/85.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/86.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/86.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/86.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/86.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/86.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/86.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/87.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/87.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/87.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/87.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/87.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/87.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/88.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/88.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/88.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/88.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/88.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/88.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/89.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/89.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/89.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/89.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/89.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/89.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/90.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/90.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/90.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/90.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/90.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/90.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/91.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/91.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/91.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/91.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/91.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/91.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/92.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/92.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/92.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/92.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/92.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/92.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/93.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/93.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/93.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/93.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/93.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/93.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/94.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/94.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/94.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/94.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/94.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/94.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/95.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/95.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/95.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/95.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/95.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/95.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/96.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/96.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/96.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/96.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/96.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/96.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/97.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/97.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/97.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/97.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/97.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/97.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/98.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/98.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/98.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/98.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/98.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/98.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_10/99.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/99.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/99.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/99.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/99.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_10/99.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/00.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/00.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/00.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/00.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/00.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/00.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/01.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/01.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/01.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/01.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/01.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/01.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/02.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/02.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/02.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/02.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/02.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/02.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/03.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/03.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/03.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/03.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/03.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/03.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/04.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/04.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/04.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/04.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/04.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/04.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/05.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/05.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/05.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/05.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/05.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/05.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/06.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/06.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/06.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/06.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/06.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/06.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/07.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/07.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/07.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/07.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/07.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/07.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/08.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/08.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/08.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/08.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/08.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/08.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/09.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/09.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/09.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/09.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/09.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/09.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/10.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/10.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/10.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/10.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/10.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/10.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/11.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/11.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/11.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/11.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/11.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/11.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/12.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/12.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/12.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/12.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/12.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/12.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/13.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/13.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/13.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/13.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/13.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/13.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/14.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/14.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/14.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/14.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/14.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/14.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/15.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/15.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/15.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/15.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/15.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/15.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/16.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/16.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/16.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/16.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/16.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/16.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/17.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/17.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/17.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/17.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/17.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/17.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/18.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/18.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/18.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/18.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/18.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/18.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/19.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/19.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/19.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/19.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/19.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/19.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/20.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/20.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/20.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/20.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/20.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/20.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/21.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/21.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/21.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/21.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/21.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/21.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/22.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/22.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/22.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/22.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/22.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/22.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/23.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/23.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/23.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/23.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/23.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/23.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/24.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/24.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/24.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/24.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/24.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/24.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/25.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/25.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/25.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/25.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/25.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/25.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/26.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/26.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/26.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/26.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/26.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/26.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/27.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/27.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/27.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/27.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/27.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/27.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/28.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/28.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/28.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/28.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/28.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/28.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/29.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/29.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/29.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/29.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/29.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/29.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/30.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/30.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/30.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/30.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/30.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/30.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/31.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/31.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/31.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/31.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/31.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/31.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/32.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/32.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/32.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/32.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/32.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/32.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/33.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/33.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/33.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/33.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/33.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/33.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/34.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/34.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/34.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/34.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/34.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/34.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/35.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/35.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/35.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/35.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/35.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/35.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/36.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/36.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/36.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/36.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/36.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/36.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/37.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/37.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/37.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/37.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/37.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/37.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/38.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/38.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/38.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/38.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/38.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/38.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/39.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/39.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/39.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/39.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/39.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/39.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/40.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/40.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/40.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/40.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/40.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/40.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/41.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/41.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/41.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/41.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/41.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/41.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/42.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/42.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/42.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/42.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/42.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/42.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/43.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/43.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/43.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/43.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/43.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/43.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/44.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/44.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/44.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/44.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/44.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/44.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/45.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/45.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/45.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/45.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/45.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/45.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/46.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/46.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/46.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/46.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/46.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/46.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/47.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/47.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/47.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/47.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/47.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/47.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/48.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/48.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/48.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/48.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/48.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/48.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/49.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/49.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/49.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/49.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/49.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/49.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/50.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/50.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/50.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/50.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/50.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/50.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/51.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/51.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/51.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/51.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/51.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/51.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/52.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/52.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/52.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/52.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/52.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/52.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/53.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/53.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/53.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/53.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/53.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/53.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/54.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/54.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/54.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/54.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/54.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/54.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/55.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/55.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/55.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/55.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/55.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/55.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/56.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/56.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/56.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/56.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/56.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/56.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/57.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/57.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/57.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/57.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/57.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/57.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/58.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/58.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/58.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/58.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/58.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/58.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/59.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/59.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/59.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/59.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/59.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/59.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/60.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/60.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/60.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/60.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/60.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/60.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/61.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/61.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/61.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/61.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/61.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/61.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/62.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/62.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/62.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/62.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/62.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/62.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/63.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/63.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/63.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/63.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/63.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/63.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/64.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/64.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/64.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/64.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/64.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/64.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/65.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/65.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/65.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/65.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/65.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/65.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/66.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/66.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/66.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/66.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/66.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/66.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/67.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/67.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/67.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/67.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/67.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/67.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/68.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/68.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/68.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/68.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/68.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/68.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/69.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/69.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/69.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/69.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/69.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/69.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/70.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/70.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/70.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/70.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/70.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/70.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/71.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/71.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/71.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/71.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/71.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/71.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/72.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/72.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/72.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/72.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/72.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/72.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/73.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/73.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/73.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/73.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/73.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/73.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/74.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/74.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/74.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/74.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/74.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/74.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/75.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/75.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/75.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/75.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/75.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/75.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/76.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/76.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/76.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/76.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/76.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/76.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/77.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/77.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/77.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/77.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/77.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/77.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/78.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/78.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/78.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/78.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/78.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/78.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/79.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/79.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/79.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/79.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/79.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/79.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/80.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/80.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/80.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/80.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/80.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/80.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/81.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/81.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/81.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/81.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/81.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/81.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/82.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/82.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/82.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/82.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/82.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/82.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/83.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/83.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/83.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/83.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/83.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/83.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/84.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/84.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/84.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/84.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/84.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/84.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/85.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/85.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/85.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/85.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/85.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/85.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/86.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/86.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/86.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/86.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/86.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/86.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/87.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/87.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/87.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/87.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/87.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/87.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/88.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/88.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/88.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/88.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/88.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/88.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/89.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/89.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/89.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/89.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/89.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/89.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/90.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/90.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/90.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/90.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/90.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/90.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/91.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/91.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/91.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/91.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/91.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/91.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/92.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/92.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/92.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/92.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/92.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/92.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/93.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/93.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/93.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/93.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/93.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/93.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/94.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/94.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/94.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/94.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/94.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/94.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/95.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/95.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/95.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/95.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/95.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/95.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/96.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/96.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/96.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/96.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/96.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/96.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/97.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/97.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/97.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/97.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/97.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/97.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/98.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/98.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/98.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/98.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/98.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/98.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_20/99.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/99.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/99.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/99.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/99.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_20/99.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/00.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/00.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/00.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/00.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/00.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/00.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/01.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/01.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/01.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/01.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/01.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/01.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/02.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/02.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/02.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/02.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/02.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/02.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/03.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/03.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/03.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/03.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/03.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/03.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/04.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/04.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/04.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/04.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/04.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/04.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/05.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/05.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/05.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/05.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/05.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/05.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/06.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/06.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/06.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/06.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/06.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/06.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/07.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/07.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/07.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/07.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/07.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/07.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/08.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/08.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/08.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/08.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/08.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/08.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/09.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/09.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/09.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/09.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/09.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/09.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/10.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/10.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/10.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/10.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/10.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/10.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/11.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/11.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/11.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/11.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/11.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/11.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/12.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/12.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/12.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/12.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/12.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/12.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/13.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/13.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/13.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/13.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/13.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/13.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/14.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/14.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/14.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/14.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/14.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/14.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/15.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/15.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/15.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/15.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/15.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/15.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/16.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/16.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/16.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/16.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/16.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/16.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/17.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/17.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/17.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/17.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/17.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/17.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/18.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/18.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/18.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/18.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/18.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/18.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/19.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/19.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/19.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/19.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/19.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/19.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/20.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/20.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/20.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/20.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/20.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/20.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/21.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/21.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/21.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/21.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/21.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/21.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/22.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/22.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/22.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/22.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/22.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/22.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/23.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/23.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/23.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/23.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/23.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/23.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/24.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/24.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/24.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/24.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/24.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/24.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/25.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/25.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/25.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/25.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/25.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/25.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/26.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/26.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/26.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/26.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/26.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/26.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/27.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/27.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/27.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/27.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/27.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/27.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/28.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/28.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/28.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/28.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/28.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/28.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/29.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/29.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/29.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/29.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/29.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/29.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/30.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/30.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/30.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/30.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/30.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/30.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/31.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/31.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/31.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/31.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/31.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/31.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/32.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/32.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/32.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/32.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/32.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/32.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/33.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/33.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/33.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/33.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/33.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/33.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/34.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/34.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/34.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/34.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/34.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/34.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/35.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/35.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/35.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/35.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/35.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/35.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/36.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/36.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/36.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/36.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/36.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/36.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/37.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/37.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/37.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/37.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/37.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/37.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/38.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/38.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/38.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/38.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/38.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/38.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/39.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/39.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/39.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/39.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/39.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/39.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/40.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/40.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/40.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/40.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/40.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/40.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/41.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/41.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/41.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/41.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/41.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/41.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/42.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/42.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/42.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/42.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/42.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/42.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/43.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/43.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/43.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/43.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/43.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/43.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/44.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/44.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/44.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/44.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/44.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/44.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/45.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/45.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/45.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/45.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/45.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/45.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/46.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/46.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/46.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/46.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/46.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/46.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/47.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/47.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/47.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/47.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/47.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/47.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/48.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/48.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/48.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/48.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/48.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/48.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/49.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/49.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/49.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/49.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/49.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/49.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/50.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/50.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/50.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/50.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/50.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/50.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/51.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/51.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/51.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/51.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/51.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/51.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/52.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/52.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/52.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/52.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/52.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/52.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/53.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/53.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/53.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/53.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/53.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/53.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/54.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/54.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/54.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/54.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/54.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/54.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/55.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/55.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/55.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/55.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/55.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/55.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/56.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/56.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/56.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/56.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/56.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/56.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/57.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/57.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/57.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/57.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/57.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/57.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/58.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/58.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/58.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/58.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/58.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/58.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/59.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/59.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/59.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/59.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/59.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/59.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/60.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/60.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/60.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/60.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/60.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/60.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/61.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/61.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/61.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/61.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/61.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/61.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/62.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/62.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/62.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/62.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/62.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/62.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/63.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/63.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/63.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/63.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/63.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/63.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/64.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/64.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/64.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/64.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/64.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/64.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/65.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/65.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/65.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/65.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/65.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/65.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/66.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/66.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/66.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/66.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/66.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/66.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/67.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/67.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/67.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/67.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/67.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/67.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/68.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/68.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/68.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/68.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/68.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/68.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/69.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/69.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/69.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/69.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/69.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/69.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/70.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/70.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/70.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/70.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/70.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/70.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/71.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/71.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/71.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/71.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/71.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/71.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/72.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/72.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/72.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/72.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/72.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/72.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/73.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/73.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/73.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/73.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/73.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/73.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/74.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/74.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/74.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/74.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/74.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/74.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/75.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/75.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/75.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/75.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/75.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/75.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/76.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/76.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/76.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/76.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/76.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/76.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/77.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/77.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/77.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/77.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/77.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/77.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/78.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/78.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/78.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/78.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/78.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/78.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/79.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/79.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/79.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/79.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/79.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/79.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/80.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/80.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/80.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/80.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/80.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/80.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/81.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/81.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/81.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/81.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/81.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/81.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/82.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/82.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/82.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/82.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/82.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/82.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/83.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/83.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/83.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/83.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/83.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/83.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/84.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/84.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/84.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/84.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/84.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/84.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/85.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/85.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/85.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/85.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/85.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/85.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/86.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/86.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/86.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/86.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/86.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/86.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/87.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/87.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/87.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/87.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/87.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/87.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/88.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/88.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/88.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/88.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/88.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/88.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/89.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/89.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/89.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/89.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/89.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/89.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/90.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/90.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/90.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/90.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/90.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/90.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/91.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/91.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/91.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/91.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/91.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/91.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/92.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/92.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/92.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/92.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/92.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/92.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/93.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/93.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/93.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/93.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/93.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/93.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/94.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/94.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/94.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/94.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/94.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/94.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/95.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/95.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/95.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/95.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/95.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/95.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/96.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/96.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/96.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/96.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/96.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/96.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/97.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/97.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/97.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/97.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/97.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/97.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/98.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/98.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/98.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/98.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/98.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/98.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_30/99.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/99.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/99.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/99.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/99.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_30/99.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/00.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/00.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/00.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/00.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/00.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/00.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/01.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/01.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/01.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/01.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/01.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/01.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/02.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/02.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/02.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/02.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/02.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/02.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/03.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/03.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/03.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/03.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/03.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/03.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/04.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/04.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/04.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/04.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/04.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/04.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/05.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/05.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/05.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/05.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/05.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/05.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/06.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/06.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/06.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/06.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/06.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/06.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/07.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/07.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/07.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/07.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/07.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/07.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/08.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/08.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/08.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/08.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/08.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/08.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/09.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/09.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/09.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/09.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/09.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/09.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/10.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/10.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/10.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/10.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/10.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/10.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/11.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/11.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/11.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/11.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/11.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/11.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/12.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/12.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/12.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/12.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/12.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/12.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/13.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/13.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/13.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/13.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/13.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/13.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/14.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/14.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/14.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/14.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/14.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/14.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/15.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/15.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/15.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/15.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/15.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/15.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/16.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/16.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/16.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/16.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/16.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/16.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/17.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/17.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/17.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/17.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/17.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/17.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/18.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/18.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/18.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/18.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/18.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/18.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/19.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/19.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/19.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/19.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/19.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/19.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/20.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/20.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/20.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/20.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/20.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/20.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/21.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/21.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/21.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/21.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/21.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/21.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/22.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/22.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/22.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/22.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/22.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/22.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/23.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/23.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/23.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/23.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/23.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/23.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/24.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/24.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/24.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/24.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/24.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/24.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/25.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/25.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/25.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/25.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/25.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/25.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/26.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/26.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/26.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/26.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/26.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/26.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/27.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/27.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/27.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/27.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/27.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/27.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/28.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/28.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/28.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/28.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/28.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/28.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/29.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/29.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/29.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/29.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/29.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/29.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/30.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/30.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/30.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/30.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/30.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/30.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/31.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/31.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/31.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/31.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/31.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/31.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/32.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/32.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/32.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/32.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/32.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/32.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/33.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/33.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/33.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/33.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/33.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/33.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/34.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/34.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/34.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/34.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/34.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/34.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/35.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/35.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/35.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/35.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/35.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/35.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/36.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/36.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/36.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/36.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/36.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/36.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/37.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/37.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/37.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/37.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/37.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/37.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/38.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/38.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/38.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/38.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/38.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/38.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/39.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/39.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/39.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/39.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/39.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/39.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/40.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/40.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/40.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/40.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/40.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/40.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/41.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/41.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/41.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/41.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/41.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/41.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/42.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/42.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/42.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/42.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/42.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/42.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/43.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/43.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/43.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/43.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/43.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/43.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/44.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/44.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/44.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/44.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/44.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/44.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/45.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/45.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/45.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/45.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/45.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/45.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/46.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/46.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/46.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/46.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/46.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/46.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/47.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/47.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/47.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/47.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/47.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/47.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/48.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/48.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/48.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/48.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/48.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/48.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/49.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/49.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/49.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/49.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/49.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/49.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/50.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/50.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/50.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/50.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/50.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/50.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/51.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/51.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/51.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/51.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/51.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/51.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/52.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/52.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/52.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/52.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/52.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/52.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/53.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/53.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/53.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/53.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/53.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/53.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/54.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/54.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/54.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/54.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/54.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/54.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/55.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/55.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/55.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/55.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/55.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/55.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/56.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/56.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/56.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/56.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/56.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/56.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/57.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/57.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/57.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/57.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/57.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/57.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/58.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/58.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/58.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/58.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/58.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/58.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/59.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/59.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/59.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/59.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/59.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/59.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/60.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/60.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/60.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/60.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/60.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/60.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/61.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/61.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/61.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/61.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/61.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/61.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/62.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/62.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/62.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/62.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/62.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/62.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/63.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/63.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/63.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/63.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/63.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/63.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/64.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/64.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/64.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/64.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/64.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/64.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/65.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/65.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/65.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/65.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/65.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/65.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/66.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/66.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/66.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/66.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/66.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/66.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/67.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/67.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/67.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/67.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/67.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/67.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/68.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/68.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/68.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/68.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/68.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/68.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/69.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/69.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/69.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/69.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/69.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/69.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/70.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/70.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/70.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/70.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/70.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/70.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/71.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/71.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/71.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/71.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/71.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/71.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/72.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/72.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/72.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/72.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/72.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/72.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/73.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/73.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/73.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/73.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/73.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/73.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/74.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/74.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/74.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/74.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/74.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/74.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/75.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/75.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/75.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/75.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/75.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/75.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/76.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/76.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/76.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/76.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/76.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/76.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/77.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/77.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/77.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/77.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/77.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/77.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/78.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/78.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/78.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/78.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/78.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/78.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/79.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/79.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/79.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/79.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/79.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/79.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/80.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/80.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/80.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/80.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/80.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/80.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/81.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/81.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/81.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/81.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/81.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/81.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/82.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/82.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/82.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/82.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/82.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/82.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/83.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/83.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/83.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/83.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/83.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/83.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/84.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/84.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/84.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/84.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/84.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/84.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/85.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/85.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/85.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/85.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/85.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/85.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/86.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/86.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/86.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/86.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/86.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/86.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/87.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/87.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/87.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/87.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/87.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/87.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/88.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/88.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/88.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/88.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/88.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/88.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/89.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/89.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/89.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/89.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/89.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/89.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/90.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/90.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/90.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/90.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/90.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/90.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/91.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/91.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/91.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/91.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/91.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/91.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/92.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/92.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/92.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/92.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/92.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/92.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/93.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/93.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/93.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/93.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/93.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/93.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/94.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/94.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/94.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/94.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/94.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/94.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/95.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/95.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/95.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/95.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/95.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/95.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/96.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/96.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/96.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/96.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/96.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/96.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/97.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/97.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/97.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/97.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/97.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/97.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/98.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/98.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/98.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/98.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/98.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/98.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_40/99.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/99.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/99.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/99.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/99.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_40/99.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/00.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/00.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/00.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/00.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/00.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/00.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/01.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/01.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/01.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/01.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/01.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/01.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/02.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/02.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/02.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/02.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/02.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/02.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/03.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/03.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/03.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/03.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/03.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/03.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/04.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/04.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/04.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/04.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/04.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/04.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/05.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/05.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/05.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/05.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/05.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/05.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/06.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/06.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/06.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/06.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/06.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/06.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/07.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/07.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/07.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/07.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/07.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/07.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/08.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/08.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/08.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/08.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/08.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/08.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/09.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/09.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/09.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/09.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/09.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/09.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/10.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/10.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/10.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/10.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/10.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/10.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/11.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/11.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/11.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/11.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/11.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/11.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/12.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/12.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/12.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/12.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/12.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/12.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/13.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/13.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/13.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/13.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/13.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/13.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/14.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/14.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/14.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/14.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/14.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/14.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/15.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/15.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/15.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/15.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/15.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/15.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/16.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/16.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/16.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/16.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/16.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/16.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/17.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/17.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/17.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/17.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/17.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/17.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/18.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/18.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/18.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/18.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/18.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/18.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/19.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/19.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/19.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/19.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/19.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/19.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/20.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/20.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/20.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/20.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/20.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/20.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/21.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/21.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/21.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/21.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/21.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/21.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/22.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/22.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/22.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/22.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/22.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/22.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/23.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/23.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/23.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/23.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/23.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/23.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/24.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/24.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/24.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/24.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/24.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/24.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/25.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/25.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/25.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/25.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/25.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/25.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/26.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/26.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/26.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/26.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/26.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/26.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/27.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/27.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/27.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/27.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/27.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/27.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/28.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/28.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/28.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/28.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/28.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/28.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/29.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/29.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/29.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/29.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/29.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/29.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/30.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/30.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/30.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/30.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/30.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/30.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/31.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/31.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/31.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/31.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/31.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/31.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/32.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/32.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/32.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/32.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/32.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/32.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/33.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/33.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/33.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/33.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/33.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/33.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/34.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/34.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/34.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/34.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/34.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/34.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/35.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/35.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/35.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/35.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/35.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/35.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/36.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/36.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/36.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/36.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/36.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/36.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/37.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/37.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/37.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/37.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/37.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/37.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/38.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/38.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/38.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/38.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/38.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/38.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/39.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/39.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/39.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/39.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/39.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/39.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/40.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/40.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/40.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/40.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/40.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/40.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/41.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/41.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/41.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/41.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/41.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/41.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/42.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/42.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/42.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/42.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/42.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/42.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/43.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/43.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/43.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/43.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/43.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/43.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/44.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/44.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/44.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/44.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/44.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/44.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/45.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/45.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/45.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/45.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/45.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/45.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/46.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/46.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/46.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/46.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/46.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/46.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/47.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/47.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/47.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/47.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/47.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/47.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/48.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/48.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/48.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/48.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/48.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/48.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/49.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/49.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/49.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/49.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/49.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/49.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/50.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/50.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/50.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/50.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/50.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/50.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/51.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/51.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/51.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/51.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/51.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/51.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/52.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/52.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/52.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/52.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/52.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/52.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/53.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/53.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/53.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/53.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/53.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/53.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/54.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/54.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/54.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/54.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/54.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/54.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/55.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/55.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/55.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/55.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/55.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/55.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/56.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/56.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/56.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/56.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/56.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/56.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/57.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/57.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/57.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/57.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/57.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/57.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/58.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/58.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/58.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/58.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/58.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/58.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/59.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/59.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/59.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/59.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/59.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/59.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/60.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/60.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/60.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/60.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/60.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/60.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/61.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/61.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/61.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/61.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/61.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/61.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/62.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/62.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/62.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/62.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/62.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/62.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/63.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/63.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/63.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/63.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/63.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/63.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/64.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/64.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/64.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/64.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/64.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/64.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/65.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/65.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/65.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/65.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/65.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/65.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/66.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/66.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/66.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/66.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/66.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/66.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/67.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/67.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/67.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/67.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/67.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/67.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/68.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/68.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/68.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/68.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/68.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/68.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/69.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/69.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/69.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/69.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/69.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/69.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/70.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/70.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/70.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/70.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/70.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/70.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/71.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/71.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/71.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/71.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/71.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/71.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/72.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/72.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/72.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/72.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/72.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/72.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/73.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/73.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/73.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/73.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/73.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/73.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/74.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/74.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/74.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/74.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/74.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/74.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/75.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/75.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/75.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/75.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/75.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/75.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/76.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/76.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/76.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/76.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/76.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/76.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/77.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/77.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/77.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/77.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/77.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/77.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/78.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/78.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/78.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/78.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/78.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/78.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/79.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/79.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/79.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/79.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/79.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/79.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/80.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/80.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/80.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/80.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/80.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/80.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/81.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/81.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/81.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/81.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/81.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/81.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/82.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/82.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/82.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/82.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/82.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/82.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/83.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/83.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/83.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/83.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/83.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/83.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/84.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/84.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/84.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/84.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/84.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/84.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/85.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/85.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/85.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/85.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/85.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/85.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/86.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/86.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/86.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/86.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/86.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/86.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/87.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/87.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/87.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/87.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/87.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/87.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/88.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/88.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/88.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/88.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/88.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/88.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/89.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/89.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/89.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/89.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/89.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/89.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/90.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/90.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/90.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/90.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/90.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/90.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/91.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/91.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/91.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/91.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/91.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/91.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/92.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/92.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/92.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/92.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/92.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/92.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/93.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/93.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/93.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/93.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/93.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/93.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/94.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/94.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/94.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/94.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/94.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/94.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/95.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/95.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/95.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/95.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/95.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/95.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/96.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/96.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/96.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/96.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/96.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/96.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/97.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/97.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/97.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/97.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/97.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/97.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/98.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/98.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/98.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/98.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/98.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/98.sam

bwa mem /labs/Aguiar/pHapCompass/references/AWRI1499/GCA_000259595.1/GCA_000259595.1_AWRI1499_v1.0_genomic.fna /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/fastq/cov_50/99.fastq > /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/99.sam
samtools view -Sb /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/99.sam | samtools sort -o /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/99.bam
samtools index /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/99.bam
rm /labs/Aguiar/pHapCompass/datasets/SRR942191/simulated/bam/cov_50/99.sam

