# pHapCompass: Probabilistic Polyploid Haplotype Assembly

pHapCompass is a unified probabilistic framework for **polyploid haplotype assembly** supporting both  
**short-read** and **long-read** sequencing data. 


---

# 1. Installation

pHapCompass requires **Python 3.10+** and standard build tools (gcc, make).


### Option 1: Install from GitHub (Recommended)

```bash
# Clone the repository
git clone https://github.com/bayesomicslab/pHapCompass.git
cd pHapCompass

# Compile extractHAIRS (required for BAM input)
cd third_party/extract_poly
make
cd ../..

# Install pHapCompass
pip install -e .
```

### Option 2: Install with pre-computed fragment files

If you already have fragment files (`.frag`) or plan to use extractHAIRS separately:

```bash
pip install git+https://github.com/bayesomicslab/pHapCompass.git
```

Then either:
- Install extractHAIRS separately and add to PATH, OR
- Use `--frag-path` to provide pre-computed fragment files

### Requirements

- Python 3.10 or higher
- gcc and make (for compiling extractHAIRS)
- Python packages: pandas, scipy, opt_einsum, pysam (installed automatically)

### Verifying Installation

```bash
# Check that pHapCompass is installed
phapcompass --help

# Test with example data
phapcompass --data-type short \
  --bam-path test_data/short_data_example/0.bam \
  --vcf-path test_data/ref_example/Chr1_unphased.vcf \
  --result-path output.vcf.gz
```


---

# 2. Input Requirements

To run pHapCompass, you need:

### **Required**
- **BAM file**: aligned reads from a single individual  
- **VCF file**: containing heterozygous SNPs (biallelic or multiallelic)

### **Optional**
- A pre‑computed `.frag` fragment file. 

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

Optionally, you may adjust the hyperparameters of the short-read model:

- `--mw` : MEC weight (default: 10.0)
- `--lw` : likelihood weight (default: 1.0)
- `--sw` : FFBS sample weight (default: 1.0)
- `--epsilon` : sequencing error rate (default: 1e-5)
- `--uncertainty [N]` : enable N-sample FFBS solution sampling (default N=3 when no value is provided)

Example with custom parameters:

```bash
phapcompass  --data-type short --bam-path sample.bam  --vcf-path sample.vcf.gz  --result-path output_short.vcf.gz   --mw 8 --lw 2 --sw 0.5 
```
Note: The weights do not have to be between 0 and 1 as inputs.

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

Optionally, you may adjust the hyperparameters of the long-read model:

- `--delta` : transition penalty parameter (default: 5)
- `--learning-rate` : optimization learning rate (default: 0.02)
- `--epsilon` : sequencing error rate (default: 1e-5)
- `--uncertainty [N]` : enable N-sample solution sampling (default N=3 when no value is provided)

Example with custom parameters:

```bash
phapcompass --data-type long --bam-path sample.bam  --vcf-path sample.vcf.gz  --result-path output_long.vcf.gz  --delta 4 --learning-rate 0.01 --epsilon 0.00002 
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

# 5. Simulation Pipeline (Haplotype References and Optional Reads)

pHapCompass includes a simulator for generating polyploid haplotype references (and optionally reads) for benchmarking.
Read simulation depends on simulated haplotypes, so the pipeline is organized as:

1. **Haplotype simulation** (required)
2. **Read simulation** (optional; uses output of step 1)

### 5.1 Simulate haplotype references

Help:
```bash
phapcompass simulation haplotypes -h
```

Autopolyploidy example:
```bash
phapcompass simulation haplotypes \
  --reference_path reference/potato_tetra/He1_Chr1_only.fasta \
  --output_dir sim_out \
  --structure autopolyploidy \
  --num_samples 1 \
  --ploidies 4 \
  --mutation_rates 0.001
```

Allopolyploidy example:
```bash
phapcompass simulation haplotypes \
  --reference_path reference/potato_tetra/He1_Chr1_only.fasta \
  --output_dir sim_out \
  --structure allopolyploidy \
  --num_samples 1 \
  --sg_rates 0.0005 0.0001 \
  --mutation_rates 0.00005 0.0001
```

**Reference length note:** the haplotype simulator uses a fixed region window (`500000–1000000`) when `shifted=True`.
Ensure your reference contig is long enough, or adjust the window in `src/phapcompass/simulator/simulate_haplotypes.py`.

### 5.2 Simulate reads (planned)

Read simulation is under development and will be exposed through the same pipeline entry. It will always run **after** haplotype simulation.

---

# 6. Evaluation (VER, MEC, Geometric MEC)

pHapCompass includes utilities to evaluate predicted polyploid haplotypes against truth.

Help:
```bash
phapcompass eval -h
```

### 6.1 VER (Vector Error Rate)

```bash
phapcompass eval ver \
  --truth-vcf path/to/truth.vcf.gz \
  --pred-vcf path/to/pred.vcf.gz \
  --ploidy 4
```

### 6.2 MEC (Minimum Error Correction)

Using a fragment file:
```bash
phapcompass eval mec \
  --pred-vcf path/to/pred.vcf.gz \
  --ploidy 4 \
  --frag path/to/reads.frag
```

Using a BAM file (runs extractHAIRS; requires the VCF used for extraction):
```bash
phapcompass eval mec \
  --pred-vcf path/to/pred.vcf.gz \
  --ploidy 4 \
  --bam path/to/reads.bam \
  --vcf path/to/input_unphased.vcf.gz
```

### 6.3 Geometric MEC

```bash
phapcompass eval geom-mec \
  --pred-vcf path/to/pred.vcf.gz \
  --ploidy 4 \
  --frag path/to/reads.frag
```

---

# 7. Output Format (Updated VCF Specification)

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

# 8. Example Output (Truncated)

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

# 9. Availability of Simulated Datasets

A subset of our simulated polyploid benchmarking data is publicly available:

**Zenodo dataset:**  
https://zenodo.org/records/17667753

The remaining datasets will be released upon acceptance of the manuscript.

---

# 10. Citation

If you use pHapCompass, please cite our preprint:

**Hosseini et al.**  
*pHapCompass: Probabilistic Assembly and Uncertainty Quantification of Polyploid Haplotype Phase*  
arXiv:2512.04393  
https://doi.org/10.48550/arXiv.2512.04393

BibTeX:

```
@article{hosseini2025phapcompass,
  title={pHapCompass: Probabilistic Assembly and Uncertainty Quantification of Polyploid Haplotype Phase},
  author={Hosseini, Marjan and Veiner, Ella and Bergendahl, Thomas and Yasenpoor, Tala and Smith, Zane and Staton, Margaret and Aguiar, Derek},
  journal={arXiv preprint arXiv:2512.04393},
  year={2025}
}
```

---

# 11. Contact

For questions or issues, please open a GitHub issue on the project repository.
