# pHapCompass

pHapCompass is a probabilistic haplotype assembly framework for polyploid genomes.  
It supports both **short-read** and **long-read** data and can optionally perform **uncertainty quantification** over haplotypes.

---

## Installation

### Option 1: Using Conda (recommended)

```bash
conda env create -f environment.yaml
conda activate hap
```

### Option 2: Using pip

```bash
pip install -r requirements.txt
```

Note: If using BAM + VCF input, you must install **extractHAIRS** separately or provide its path via `--extracthairs-bin`.

---

## Quick Start

### Short-read example (BAM + VCF)

```bash
python main.py   --data-type short   --bam-path sample.bam   --vcf-path sample.vcf.gz   --result-path outputs/sample_short.txt   --mw 1.0 --lw 1.0 --sw 1.0   --epsilon 0.01   --mbq 13   --uncertainty 5   --verbose
```

### Long-read example (fragment file)

```bash
python main.py   --data-type long   --fragment-file-path sample.frag   --vcf-path sample.vcf.gz   --result-path outputs/sample_long.txt   --delta 0.1   --learning-rate 0.001   --uncertainty
```

---

## Arguments

### Core I/O

- `--bam-path PATH`  
- `--fragment-file-path PATH`  
- `--vcf-path PATH` (required)  
- `--ploidy INT`  
- `--result-path PATH` (required)

### Model Type

- `--data-type {short,long}`  
- `--epsilon FLOAT`

### Short-read hyperparameters

- `--mw FLOAT`  
- `--lw FLOAT`  
- `--sw FLOAT`

### Long-read hyperparameters

- `--delta FLOAT`  
- `--learning-rate FLOAT`

### Fragment extraction

- `--mbq INT`  
- `--extracthairs-bin PATH`

### Uncertainty quantification

- `--uncertainty [N]`  
  - no flag → disabled  
  - `--uncertainty` → 3 samples  
  - `--uncertainty N` → N samples  

### Logging

- `--log-level LEVEL`  
- `--verbose`

---

## Outputs

- `genotype.csv` (auto-generated)  
- Fragment file(s) (*.frag)  
- Main haplotype output (`--result-path`)  

---

## Citation


