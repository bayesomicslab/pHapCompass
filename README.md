# pHapCompass

pHapCompass is a probabilistic haplotype assembly framework for polyploid genomes.  
It supports both **short-read** and **long-read** data and optionally performs **uncertainty quantification** across haplotype solutions.

---

## Installation

pHapCompass is a Python package targeting **Python 3.10**.

### Option 1 — Using Conda + Pip (recommended)

Create and activate a dedicated environment:

```bash
conda create -n phapcompass python=3.10 -y
conda activate phapcompass
```

Then install pHapCompass:

```bash
pip install "git+https://github.com/bayesomicslab/pHapCompass.git"
```

### Option 2 — From a local clone

```bash
git clone https://github.com/bayesomicslab/pHapCompass.git
cd pHapCompass

conda create -n phapcompass python=3.10 -y
conda activate phapcompass

pip install .
```

### Option 3 — Using the provided environment / requirements files

If you prefer using the provided files:

```bash
conda env create -f environment.yml
conda activate hap
pip install .
```

or, if you already have a suitable Python 3.10 environment:

```bash
pip install -r requirements.txt
pip install .
```

> **Note:** If starting from BAM + VCF, you must install **extractHAIRS** and ensure it is available on your `$PATH` (or provide `--extracthairs-bin`).

> **Conda package:** A pure `conda install phapcompass` (e.g., via bioconda) is not yet available. Once a bioconda recipe is merged, installation will be possible via:
> ```bash
> conda install -c bioconda phapcompass
> ```

---

## Quick Start

Below are examples for **both models**, with **both input modes** (fragment file OR BAM+VCF).  
Replace paths with your own data.

> In all examples below, `phapcompass` refers to the installed CLI.  

### 🔹 Short-read model (data-type = short)

#### **Option A — Start from BAM + VCF**

```bash
phapcompass --data-type short  --bam-path sample.bam --vcf-path sample.vcf.gz --result-path outputs/sample_short.vcf.gz --extracthairs-bin extractHAIRS
```

#### **Option B — Start from a precomputed fragment file**

```bash
phapcompass --data-type short --frag-path sample.frag --vcf-path sample.vcf.gz --result-path outputs/sample_short.vcf.gz
```

---

### 🔹 Long-read model (data-type = long)


#### **Option A — Start from BAM + VCF**

```bash
phapcompass   --data-type long   --bam-path sample.bam   --vcf-path sample.vcf.gz   --result-path outputs/sample_long.vcf.gz   --delta 0.1   --learning-rate 0.001   --epsilon 0.01   --mbq 13   --uncertainty 4
```

#### **Option B — Start from a precomputed fragment file**

```bash
phapcompass   --data-type long   --frag-path sample.frag   --vcf-path sample.vcf.gz   --result-path outputs/sample_long.vcf.gz   --delta 0.1   --learning-rate 0.001   --uncertainty
```

---

## Command-line Arguments

### Core I/O

- `--bam-path PATH`  
- `--frag-path PATH`  
- `--vcf-path PATH` **(required)**
- `--result-path PATH` **(required)**
- `--ploidy INT` (optional; inferred from VCF if not given)
- `--epsilon FLOAT`

### Model type

- `--data-type {short,long}` **(required)**

### Short-read hyperparameters

- `--mw FLOAT`
- `--lw FLOAT`
- `--sw FLOAT`

### Long-read hyperparameters

- `--delta FLOAT`
- `--learning-rate FLOAT`

### Fragment extraction (extractHAIRS)

- `--mbq INT` (default = 13)
- `--extracthairs-bin PATH`

### Uncertainty quantification

- `--uncertainty [N]`  
  - no flag → disabled  
  - `--uncertainty` → enabled, **3 samples**  
  - `--uncertainty N` → enabled, N samples  

### Logging

- `--log-level LEVEL`
- `--verbose`

---

# Output Format (VCF)

pHapCompass outputs a **single VCF file** (typically ending in `.vcf.gz`) that contains the **phased haplotypes**.

It preserves all header lines from the input VCF and adds:

```text
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
##FORMAT=<ID=LK,Number=1,Type=Float,Description="Phasing likelihood/confidence">
```

---

## GT — Genotype (Phased or Unphased)

- **Phased** uses `|`  
  Example: `0|0|1`
- **Unphased** uses `/`  
  Example: `0/0/1`
- Values:  
  - `0` = reference allele  
  - `1` = alternate allele  
  - `.` = missing

---

## PS — Phase Set

- Shared integer ID for SNPs belonging to the same phase block
- Typically the coordinate or an index corresponding to the first SNP in the block
- Omitted (set to `.`) for unphased positions

Example:

```text
FORMAT  GT:PS
SAMPLE  0|0|1:3529
```

All sites with `PS=3529` are phased relative to each other.

---

## LK — Phasing Likelihood (Optional)

- Added only when `--uncertainty` is enabled and likelihoods are provided.
- Higher values indicate greater confidence for that phasing solution.
- For multiple solutions, LK can be written as colon-separated values (`LK=0.87:0.82:0.90`).

---

# Uncertainty Quantification Mode

When enabled (`--uncertainty N`), the VCF FORMAT fields can contain *multiple phasing solutions* for each variant.

Conceptually:

```text
GT:PS:LK
0|0|1:3529:0.876543
0|1|0:3529:0.823456
1|0|0:3529:0.901234
```

Interpretation:

- Multiple GT strings represent different sampled phasing configurations.
- PS values indicate the phase block IDs per solution.
- LK values (if present) represent relative model confidence per configuration.

This allows users to:

- inspect phasing uncertainty,
- distinguish stable vs. ambiguous regions,
- perform downstream consensus or Bayesian aggregation.

---

# Example Output (Truncated)

```text
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
##FORMAT=<ID=LK,Number=1,Type=Float,Description="Phasing confidence">
#CHROM  POS     ID   REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
Chr1    2       .    C    G    .     PASS    .     GT      0/0/1
Chr1    3529    .    A    T    .     PASS    .     GT:PS   0|0|1:0
Chr1    3533    .    A    T    .     PASS    .     GT:PS   0|0|1:0
Chr1    3780    .    A    T    .     PASS    .     GT:PS   0|0|1:1
Chr1    3781    .    G    C    .     PASS    .     GT:PS   1|0|0:1
Chr1    5934    .    A    T    .     PASS    .     GT:PS   1|0|0:1
```

The output is compatible with:

- `bcftools`
- `vcftools`
- polyploid haplotype analysis pipelines
- visualization tools

The VCF can be compressed and indexed with:

```bash
bgzip output.vcf
tabix -p vcf output.vcf.gz
```

---

# Citation

(Add citation once available.)
