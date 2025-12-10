# pHapCompass: Probabilistic Polyploid Haplotype Assembly

pHapCompass is a unified probabilistic framework for **polyploid haplotype assembly** supporting both  
**short-read** and **long-read** sequencing data. It integrates a Bayesian/graph‑based formulation with  
alternative decoding strategies and optional solution sampling.  

Unlike earlier tools, **pHapCompass bundles a fully polyploid-compatible extractHAIRS implementation**,  
so users **do not need to install extractHAIRS separately**. All fragment extraction happens internally.

---

# 1. Installation

pHapCompass requires **Python 3.10**.

### **Install from GitHub (recommended)**

```bash
pip install "git+https://github.com/bayesomicslab/pHapCompass.git"
```

### **Install from local clone**

```bash
git clone https://github.com/bayesomicslab/pHapCompass.git
cd pHapCompass
pip install .
```

This automatically builds and uses the bundled C binary:

```
third_party/extract_poly/build/extractHAIRS
```

No external tools are required.

---

# 2. Input Requirements

To run pHapCompass, you need:

### **Required**
- **BAM file**: aligned reads from a single individual  
- **VCF file**: containing heterozygous SNPs (biallelic or multiallelic)

### **Optional**
- A pre‑computed `.frag` fragment file (if you prefer not to use internal extractHAIRS)

The tool infers ploidy automatically from the VCF unless specified.

---

# 3. Basic Usage

The standard and most common usage is to run pHapCompass **directly from BAM + VCF**, letting the internal  
polyploid extractHAIRS generate fragments.

## **Short‑read model**

### **From BAM + VCF (recommended)**

```bash
phapcompass   --data-type short   --bam-path sample.bam   --vcf-path sample.vcf.gz   --result-path output_short.vcf.gz 
```

### **Using a precomputed fragment file**

```bash
phapcompass   --data-type short   --frag-path sample.frag   --vcf-path sample.vcf.gz   --result-path output_short.vcf.gz
```

---

## **Long‑read model**

### **From BAM + VCF**

```bash
phapcompass   --data-type long   --bam-path sample.bam   --vcf-path sample.vcf.gz   --result-path output_long.vcf.gz 
```

### **Using a precomputed fragment file**

```bash
phapcompass   --data-type long   --frag-path sample.frag   --vcf-path sample.vcf.gz   --result-path output_long.vcf.gz
```

---

# 4. Command‑line Arguments

## **Core I/O**
| Argument | Description |
|---------|-------------|
| `--bam-path PATH` | BAM file; triggers internal extractHAIRS. |
| `--frag-path PATH` | Optional: use an existing fragment file. |
| `--vcf-path PATH` | Required. Input VCF containing heterozygous SNPs. |
| `--result-path PATH` | Required. Output VCF path. |
| `--ploidy INT` | Optional. If omitted, inferred from VCF. |

## **Model selection**
- `--data-type short`
- `--data-type long`

## **Short‑read model hyperparameters**
- `--mw` MEC weight  
- `--lw` likelihood weight  
- `--sw` FFBS sample weight  

## **Long‑read model hyperparameters**
- `--delta`  
- `--learning-rate`  

## **Other**
- `--epsilon` sequencing error rate  
- `--uncertainty [N]` enable sampling mode (N samples; default = 3)  
- `--verbose`  

---

# 5. Output Format (Updated VCF Specification)

pHapCompass outputs a single **phased polyploid VCF** with:

### FORMAT fields:
```
GT   Genotype (phased or unphased)
PS   Phase‑set identifier
```

If uncertainty mode is enabled, we also add **probability headers** (one per solution):

```
##phapcompass_solution=<ID=i,Probability=p_i>
```

### **GT formatting**
- Phased alleles use **pipes**: `0|1|0`
- Unphased alleles use **slashes**: `0/1/0`
- Values correspond to REF/ALT encodings in the VCF.

### **PS formatting**
- Integer block ID for phased SNPs  
- `.` for unphased positions  

### **Multisolution output (uncertainty mode)**

If `--uncertainty N` is used:

- GT fields for different solutions appear **separated by ':'**  
- PS fields also appear **separated by ':'**  
- Probabilities appear in VCF **header only**, not per‑SNP

Example:

```
GT:PS
0|0|1:3529 : 0|1|0:3529 : 1|0|0:3529
```

---

# 6. Example Output (Truncated)

```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
##phapcompass_solution=<ID=1,Probability=0.812345>
##phapcompass_solution=<ID=2,Probability=0.187655>
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
Chr1   3529 .  A   T   .    PASS   .    GT:PS   0|0|1:0
Chr1   3781 .  G   C   .    PASS   .    GT:PS   1|0|0:0
Chr1   5934 .  A   T   .    PASS   .    GT:PS   1|0|0:0
```

---

# 7. Availability of Simulated Datasets

A subset of our simulated polyploid benchmarking data is publicly available:

**Zenodo dataset:**  
https://zenodo.org/records/17667753

The remaining datasets will be released upon acceptance of the manuscript.

---

# 8. Citation

If you use pHapCompass, please cite our preprint:

**Hosseini et al.**  
*pHapCompass: Probabilistic Polyploid Haplotype Assembly*  
arXiv:2512.04393  
https://doi.org/10.48550/arXiv.2512.04393

BibTeX:

```
@article{hosseini2025phapcompass,
  title={pHapCompass: Probabilistic Polyploid Haplotype Assembly},
  author={Hosseini, Marjan and McConnell, Devin and Aguiar, Derek},
  journal={arXiv preprint arXiv:2512.04393},
  year={2025},
  doi={10.48550/arXiv.2512.04393}
}
```

---

# 9. License

pHapCompass is released under the BSD‑2 license.  
The bundled extract_poly code is licensed under the BSD‑2 license of the Bansal Lab.

---

# 10. Contact

For questions or issues, please open a GitHub issue on the project repository.

