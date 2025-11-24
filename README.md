# pHapCompass

pHapCompass is a probabilistic haplotype assembly framework for polyploid genomes.  
It supports both **short-read** and **long-read** data and optionally performs **uncertainty quantification** across haplotype solutions.

---

## Installation

### Option 1 — Conda (recommended)

```bash
conda env create -f environment.yaml
conda activate hap
```

### Option 2 — Pip

```bash
pip install -r requirements.txt
```

> Note: If starting from BAM + VCF, you must install **extractHAIRS** and ensure it is available on your PATH (or provide `--extracthairs-bin`).

---

## Quick Start

Below are examples for **both models**, with **both input modes** (fragment file OR BAM+VCF).

---

### 🔹 Short-read model (data-type = short)

#### **Option A — Start from BAM + VCF**
```bash
python main.py   --data-type short   --bam-path sample.bam   --vcf-path sample.vcf.gz   --result-path outputs/sample_short.vcf.gz   --mw 1.0 --lw 1.0 --sw 1.0   --epsilon 0.01   --mbq 13   --uncertainty 5   --verbose
```

#### **Option B — Start from a precomputed fragment file**
```bash
python main.py   --data-type short   --fragment-file-path sample.frag   --vcf-path sample.vcf.gz   --result-path outputs/sample_short.vcf.gz   --mw 1.0 --lw 1.0 --sw 1.0   --epsilon 0.01   --uncertainty
```

---

### 🔹 Long-read model (data-type = long)

#### **Option A — Start from BAM + VCF**
```bash
python main.py   --data-type long   --bam-path sample.bam   --vcf-path sample.vcf.gz   --result-path outputs/sample_long.vcf.gz   --delta 0.1   --learning-rate 0.001   --epsilon 0.01   --mbq 13   --uncertainty 4
```

#### **Option B — Start from a precomputed fragment file**
```bash
python main.py   --data-type long   --fragment-file-path sample.frag   --vcf-path sample.vcf.gz   --result-path outputs/sample_long.vcf.gz   --delta 0.1   --learning-rate 0.001   --uncertainty
```

---

## Command-line Arguments

### Core I/O
- `--bam-path PATH`  
- `--fragment-file-path PATH`  
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

```
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

---

## PS — Phase Set

- Shared integer ID for SNPs belonging to the same phase block
- Typically the coordinate of the first SNP in that block
- Omitted for unphased positions

Example:
```
GT:PS
0|0|1:3529
```

---

## LK — Phasing Likelihood (Optional)

- Added only when `--uncertainty` is enabled
- Higher values → higher confidence

---

# Uncertainty Quantification Mode

When enabled (`--uncertainty N`), the VCF FORMAT fields contain *multiple phasing solutions*:

Example:

```
GT:PS:LK
0|0|1:3529:0.876543
0|1|0:3529:0.823456
1|0|0:3529:0.901234
```

Interpretation:

- Multiple GTs correspond to N sampled phasing solutions
- PS values repeated for each sample
- LK gives the model confidence for each sampled haplotype configuration

This allows users to:

- inspect phasing uncertainty
- identify stable vs. variable positions
- downstream weighting or consensus phasing

---

# Example Output (Truncated)

```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
##FORMAT=<ID=LK,Number=1,Type=Float,Description="Phasing likelihood">
#CHROM  POS     ID   REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
Chr1    2       .    C    G    .     PASS    .     GT      0/0/1
Chr1    3529    .    A    T    .     PASS    .     GT:PS   0|0|1:3529
Chr1    3533    .    A    T    .     PASS    .     GT:PS   0|0|1:3529
Chr1    3780    .    A    T    .     PASS    .     GT:PS   0|0|1:3529
Chr1    3781    .    G    C    .     PASS    .     GT:PS   1|0|0:3781
Chr1    5934    .    A    T    .     PASS    .     GT:PS   1|0|0:3781
```

The output is compatible with:

- `bcftools`
- `vcftools`
- polyploid haplotype analysis pipelines
- visualization tools

The VCF can be compressed with:

```bash
bgzip output.vcf
tabix -p vcf output.vcf.gz
```

---

# Citation

(Add citation once available.)
