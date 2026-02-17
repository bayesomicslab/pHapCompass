# pHapCompass: Probabilistic Polyploid Haplotype Assembly

pHapCompass is a unified probabilistic framework for **polyploid haplotype assembly** supporting both  
**short-read** and **long-read** sequencing data. 

---

## Table of Contents

- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Installation Options](#installation-options)
  - [Troubleshooting Installation](#troubleshooting-installation)
  - [Verifying Installation](#verifying-installation)
- [Input Requirements](#input-requirements)
- [Basic Usage](#basic-usage)
  - [Short-read model](#short-read-model)
  - [Long-read model](#long-read-model)
- [Command-line Arguments](#command-line-arguments)
- [Simulation Pipeline](#simulation-pipeline)
- [Evaluation](#evaluation)
- [Output Format](#output-format)
- [Example Output](#example-output)
- [Availability of Simulated Datasets](#availability-of-simulated-datasets)
- [Citation](#citation)
- [Contact](#contact)

---

## Installation

### Prerequisites

Before installing pHapCompass, ensure you have:

- **Python 3.10+**
- **C compiler** (gcc or cc)
- **make** build tool
- **git** (for downloading submodules)
- **zlib development headers** (required for compiling extractHAIRS)

> **Note:** pHapCompass has been tested on Ubuntu, Debian, Fedora, and RHEL/Rocky Linux.

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install build-essential git python3-pip zlib1g-dev libbz2-dev liblzma-dev
```

**Fedora:**
```bash
sudo dnf install -y git gcc make zlib-devel bzip2-devel xz-devel python3-pip
```

**RHEL/Rocky Linux:**
```bash
sudo dnf install -y git gcc make zlib-devel bzip2-devel xz-devel epel-release
# RHEL ships with Python 3.9 by default - install Python 3.11 explicitly:
sudo dnf install -y python3.11 python3.11-pip
# Then use python3.11 instead of pip when installing pHapCompass (see below)
```

---

### Installation Options

**Option 1: Install from GitHub (Recommended)**

This automatically compiles extractHAIRS during installation:

```bash
git clone --recursive https://github.com/bayesomicslab/pHapCompass.git
cd pHapCompass
pip install -e .
```

> The `--recursive` flag is important — it downloads required submodules (htslib and samtools).

On **RHEL/Rocky Linux**, use Python 3.11 explicitly:
```bash
python3.11 -m pip install -e .
```

**Option 2: Manual Compilation (if automatic compilation fails)**

```bash
git clone --recursive https://github.com/bayesomicslab/pHapCompass.git
cd pHapCompass

# Manually compile extractHAIRS
cd third_party/extract_poly
make
cd ../..

# Install pHapCompass
pip install -e .
```

**Option 3: Install without extractHAIRS**

If you plan to use pre-computed fragment files (`.frag`):

```bash
pip install git+https://github.com/bayesomicslab/pHapCompass.git
```

Then use `--frag-path` to provide pre-computed fragment files.

---

### Troubleshooting Installation

**"gcc: command not found" or "make: command not found"**
- Install build tools using the commands in the Prerequisites section above.

**"zlib.h: No such file or directory" during compilation**
- Install zlib development headers:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install zlib1g-dev

  # Fedora/RHEL
  sudo dnf install zlib-devel
  ```

**"Package requires a different Python: 3.9.x not in >=3.10" on RHEL/Rocky Linux**
- Install Python 3.11 and use it explicitly:
  ```bash
  sudo dnf install -y epel-release python3.11 python3.11-pip
  python3.11 -m pip install -e .
  ```
- Or use conda to manage the Python version:
  ```bash
  conda create -n phapcompass python=3.10 -y
  conda activate phapcompass
  pip install -e .
  ```

**Git submodule errors**
- Ensure you cloned with `--recursive`, or initialize manually:
  ```bash
  git submodule update --init --recursive
  ```

**"extractHAIRS binary not found" during runtime**
- Check if the binary was compiled:
  ```bash
  ls -lh src/phapcompass/bin/extractHAIRS
  ```
- If not present, manually compile as shown in Option 2 above.

For more detailed troubleshooting, see [TROUBLESHOOTING.md](TROUBLESHOOTING.md).

---

### Verifying Installation

```bash
# Check that pHapCompass is installed
phapcompass --help

# Test with example data (this also exercises extractHAIRS)
phapcompass --data-type short \
  --bam-path test_data/short_data_example/0.bam \
  --vcf-path test_data/ref_example/Chr1_unphased.vcf \
  --result-path output.vcf.gz
```

---

## Input Requirements

To run pHapCompass, you need:

**Required**
- **BAM file**: aligned reads from a single individual  
- **VCF file**: containing heterozygous SNPs (biallelic or multiallelic)

**Optional**
- A pre-computed `.frag` fragment file

The tool infers ploidy automatically from the VCF unless specified.

---

## Basic Usage

The standard usage is to run pHapCompass directly from BAM + VCF, letting the internal polyploid extractHAIRS generate fragments.

### Short-read model

**From BAM + VCF (recommended)**

```bash
phapcompass --data-type short \
  --bam-path sample.bam \
  --vcf-path sample.vcf.gz \
  --result-path output_short.vcf.gz
```

Optional hyperparameters:

- `--mw` : MEC weight (default: 10.0)
- `--lw` : likelihood weight (default: 1.0)
- `--sw` : FFBS sample weight (default: 1.0)
- `--epsilon` : sequencing error rate (default: 1e-5)
- `--uncertainty [N]` : enable N-sample FFBS solution sampling (default N=3)

Example with custom parameters:

```bash
phapcompass --data-type short \
  --bam-path sample.bam \
  --vcf-path sample.vcf.gz \
  --result-path output_short.vcf.gz \
  --mw 8 --lw 2 --sw 0.5
```

Note: The weights do not have to be between 0 and 1.

**Using a precomputed fragment file**

```bash
phapcompass --data-type short \
  --frag-path sample.frag \
  --vcf-path sample.vcf.gz \
  --result-path output_short.vcf.gz
```

---

### Long-read model

**From BAM + VCF**

```bash
phapcompass --data-type long \
  --bam-path sample.bam \
  --vcf-path sample.vcf.gz \
  --result-path output_long.vcf.gz
```

Optional hyperparameters:

- `--delta` : transition penalty parameter (default: 5)
- `--learning-rate` : optimization learning rate (default: 0.02)
- `--epsilon` : sequencing error rate (default: 1e-5)
- `--uncertainty [N]` : enable N-sample solution sampling (default N=3)

Example with custom parameters:

```bash
phapcompass --data-type long \
  --bam-path sample.bam \
  --vcf-path sample.vcf.gz \
  --result-path output_long.vcf.gz \
  --delta 4 --learning-rate 0.01 --epsilon 0.00002
```

**Using a precomputed fragment file**

```bash
phapcompass --data-type long \
  --frag-path sample.frag \
  --vcf-path sample.vcf.gz \
  --result-path output_long.vcf.gz
```

---

## Command-line Arguments

**Core I/O**
| Argument | Description |
|---------|-------------|
| `--bam-path PATH` | BAM file; triggers internal extractHAIRS. |
| `--frag-path PATH` | Optional: use an existing fragment file. |
| `--vcf-path PATH` | Required. Input VCF containing heterozygous SNPs. |
| `--result-path PATH` | Required. Output VCF path. |
| `--ploidy INT` | Optional. If omitted, inferred from VCF. |

**Model selection**
- `--data-type short`
- `--data-type long`

**Short-read model hyperparameters**
- `--mw` MEC weight  
- `--lw` likelihood weight  
- `--sw` FFBS sample weight  

**Long-read model hyperparameters**
- `--delta`  
- `--learning-rate`  

**Other**
- `--epsilon` sequencing error rate  
- `--uncertainty [N]` enable sampling mode (N samples; default = 3)  
- `--verbose`  

---

## Simulation Pipeline

pHapCompass includes a simulator for generating polyploid haplotype references (and optionally reads) for benchmarking. The pipeline is organized as:

1. **Haplotype simulation** (required)
2. **Read simulation** (optional; uses output of step 1)

**Simulate haplotype references**

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

> **Note:** the haplotype simulator uses a fixed region window (`500000–1000000`) when `shifted=True`. Ensure your reference contig is long enough, or adjust the window in `src/phapcompass/simulator/simulate_haplotypes.py`.

**Simulate reads** (planned)

Read simulation is under development and will be exposed through the same pipeline entry.

---

## Evaluation

pHapCompass includes utilities to evaluate predicted polyploid haplotypes against truth.

```bash
phapcompass eval -h
```

**VER (Vector Error Rate)**

```bash
phapcompass eval ver \
  --truth-vcf path/to/truth.vcf.gz \
  --pred-vcf path/to/pred.vcf.gz \
  --ploidy 4
```

**MEC (Minimum Error Correction)**

Using a fragment file:
```bash
phapcompass eval mec \
  --pred-vcf path/to/pred.vcf.gz \
  --ploidy 4 \
  --frag path/to/reads.frag
```

Using a BAM file:
```bash
phapcompass eval mec \
  --pred-vcf path/to/pred.vcf.gz \
  --ploidy 4 \
  --bam path/to/reads.bam \
  --vcf path/to/input_unphased.vcf.gz
```

**Geometric MEC**

```bash
phapcompass eval geom-mec \
  --pred-vcf path/to/pred.vcf.gz \
  --ploidy 4 \
  --frag path/to/reads.frag
```

---

## Output Format

pHapCompass outputs a single **phased polyploid VCF** with the following FORMAT fields:

```
GT   Genotype (phased or unphased)
PS   Phase-set identifier
```

If uncertainty mode is enabled, probability headers are added (one per solution):

```
##phapcompass_solution=<ID=i,Probability=p_i>
```

**GT formatting**
- Phased alleles use pipes: `0|1|0`
- Unphased alleles use slashes: `0/1/0`

**PS formatting**
- Integer block ID for phased SNPs  
- `.` for unphased positions  

**Multisolution output (uncertainty mode)**

If `--uncertainty N` is used, GT and PS fields for different solutions appear separated by `:`, and probabilities appear in the VCF header only:

```
GT:PS
0|0|1:3529 : 0|1|0:3529 : 1|0|0:3529
```

---

## Example Output

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

## Availability of Simulated Datasets

A subset of our simulated polyploid benchmarking data is publicly available:

**Zenodo dataset:**  
https://zenodo.org/records/17667753

The remaining datasets will be released upon acceptance of the manuscript.

---

## Citation

If you use pHapCompass, please cite our preprint:

**Hosseini et al.**  
*pHapCompass: Probabilistic Assembly and Uncertainty Quantification of Polyploid Haplotype Phase*  
arXiv:2512.04393  
https://doi.org/10.48550/arXiv.2512.04393

```
@article{hosseini2025phapcompass,
  title={pHapCompass: Probabilistic Assembly and Uncertainty Quantification of Polyploid Haplotype Phase},
  author={Hosseini, Marjan and Veiner, Ella and Bergendahl, Thomas and Yasenpoor, Tala and Smith, Zane and Staton, Margaret and Aguiar, Derek},
  journal={arXiv preprint arXiv:2512.04393},
  year={2025}
}
```

---

## Contact

For questions or issues, please open a GitHub issue on the project repository.
