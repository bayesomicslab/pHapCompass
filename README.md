# pHapCompass: Probabilistic Polyploid Haplotype Assembly

pHapCompass is a unified probabilistic framework for **polyploid haplotype assembly** supporting both  
**short-read** and **long-read** sequencing data. 


---

# 1. Installation

## Prerequisites

Before installing pHapCompass, ensure you have:

- **Python 3.10+**
- **C compiler** (gcc or cc)
- **make** build tool
- **git** (for downloading submodules)
- **zlib development headers** (required for compiling extractHAIRS)

### Installing Prerequisites by OS

> **Note:** pHapCompass has been tested on Ubuntu, Debian, Fedora, and RHEL/Rocky Linux. Windows and macOS are not officially supported.

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install build-essential git python3-pip zlib1g-dev libbz2-dev liblzma-dev
```

**Fedora:**
```bash
sudo dnf install -y git gcc make zlib-devel bzip2-devel xz-devel python3-pip
```

**RHEL/Rocky Linux/CentOS:**
```bash
sudo dnf install -y git gcc make zlib-devel bzip2-devel xz-devel epel-release
# RHEL default Python may be 3.9 - install Python 3.11 explicitly:
sudo dnf install -y python3.11 python3.11-pip
# Then use python3.11 instead of pip when installing pHapCompass:
python3.11 -m pip install -e .
```

## Installation Options

### Option 1: Install from GitHub (Recommended)

This option automatically compiles extractHAIRS during installation:

```bash
# Clone the repository
git clone --recursive https://github.com/bayesomicslab/pHapCompass.git
cd pHapCompass

# Install pHapCompass (this will compile extractHAIRS automatically)
pip install -e .
```

**Note:** The `--recursive` flag is important as it downloads required submodules (htslib and samtools).

### Option 2: Manual Compilation (If automatic compilation fails)

If the automatic compilation fails, you can compile extractHAIRS manually:

```bash
# Clone the repository
git clone --recursive https://github.com/bayesomicslab/pHapCompass.git
cd pHapCompass

# Manually compile extractHAIRS
cd third_party/extract_poly
make
cd ../..

# Install pHapCompass
pip install -e .
```

### Option 3: Install without extractHAIRS

If you plan to use pre-computed fragment files (`.frag`) or install extractHAIRS separately:

```bash
pip install git+https://github.com/bayesomicslab/pHapCompass.git
```

Then either:
- Install extractHAIRS separately and add to PATH, OR
- Use `--frag-path` to provide pre-computed fragment files

## Troubleshooting Installation

### Common Issues

**1. "gcc: command not found" or "make: command not found"**
- Install build tools using the commands in the Prerequisites section above

**2. "git submodule" errors**
- Ensure you cloned with `--recursive` flag
- Or manually initialize submodules:
  ```bash
  git submodule update --init --recursive
  ```

**3. Compilation errors related to htslib or samtools**
- Make sure git submodules are initialized (see above)
- Try cleaning and rebuilding:
  ```bash
  cd third_party/extract_poly
  make clean
  make
  ```

**4. "extractHAIRS binary not found" during runtime**
- Check if extractHAIRS was compiled:
  ```bash
  ls -lh src/phapcompass/bin/extractHAIRS
  ```
- If not present, manually compile as shown in Option 2

### Verifying Installation

```bash
# Check that pHapCompass is installed
phapcompass --help

# Check that extractHAIRS is available
phapcompass --data-type short --help

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
- `--uncertainty [N]` enable sampling mode (N samples; default = 3)  
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
