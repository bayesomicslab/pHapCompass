Extracting haplotype information from BAM and VCF file for polyploids
======

## About:
This is an edited version (under test) of [extracthairs](https://github.com/vibansal/HapCUT2) in which polyploid genomes are also allowed. The goal of this code is to generate fragment file needed for haplotyping algorithm like [hapcut](https://github.com/vibansal/hapcut2), [sdhap](https://sourceforge.net/projects/sdhap/), [althap](https://github.com/realabolfazl/AltHap), [Haptree v0.1](http://cb.csail.mit.edu/cb/haptree/), [HapMc](https://github.com/smajidian/HapMC), and [ H-popG](https://github.com/MinzhuXie/H-PoPG).



## To build:

```

git clone https://github.com/smajidian/extract_poly
cd extract_poly
make 
```

The makefile will attempt to build samtools and htslib as git submodules. The output of this step is a binary file `extractHAIRS` in folder `build`.





## Input:
It requires the following input:
- BAM file for an individual containing reads aligned to a reference genome
- VCF file containing only heterzygous SNVs . Complex SNVs and indels should be handled beforehand. 




## Run for Illumina dataset:


(1) Filtering VCF file (removing homozygous and non-SNP variants) for tetraploid

```
cd test_data

cat vars.vcf | grep -v "1/1/1/1" | grep -v "0/0/0/0" | grep -v "mnp" > vars_filtered.vcf
```



(2) Using extractHAIRS to convert BAM file to the compact fragment file format containing only haplotype-relevant information. 

```
../build/extractHAIRS  --bam reads_sorted.bam --VCF vars_filtered.vcf --out fragment_file
```


(3) If you need to use the fragment file for sdhap and althap, use

```
python2 ../FragmentPoly.py -f fragment_file  -o fragment_file_sdhap -x SDhaP 
```


or for Haptree v0.1. Note that haptree v1 is only for diploid.
```
python2 ../FragmentPoly.py -f fragment_file  -o fragment_file_haptree -x HapTree 
```




## Run for 10x dataset:

*important* Please run the pipeline for each chromosome separately.


(1) Filtering VCF file

```
cat variants.vcf | grep -v "0/0/0" | grep -v "1/1/1/1" | grep -v "0/0/0/0" | grep -v "mnp" > variants_filtered.vcf
```


(2) use extractHAIRS to convert BAM file to the compact fragment file format containing only haplotype-relevant information. 

```
./build/extractHAIRS --10X 1 --bam reads_sorted.bam --VCF variants_filtered.vcf --out unlinked_fragment_file
```

(3) Link fragments into barcode-specific fragment:
```
python3 utilities/LinkFragments_brcd_based.py  unlinked_fragment_file linked_fragment_file
```


(4) If you need to use the fragment file for sdhap and althap, use

```
python2 $fragpoly -f fragment_file  -o fragment_file_sdhap -x SDhaP 
```


or for 
```
python2 $fragpoly -f fragment_file  -o fragment_file_haptree -x HapTree 
```





NOTE: It is required that the BAM reads have the BX (corrected barcode) tag.






## Citation:

[Extracthairs](https://github.com/vibansal/HapCUT2)

[Haplosim](https://github.com/EhsanMotazedi/Haplosim)

[HapMc](https://github.com/smajidian/HapMC)

[Hap10] My paper under preparation.





