# GWAS Pipeline

## Overview

This folder contains scripts used for genome-wide association studies (GWAS), meta-analysis, SNP pruning, and downstream preprocessing for functional annotation.

The pipeline was developed for cardiac radiomic feature GWAS analyses and downstream genomic interpretation.

---

## Contents

| File / Folder | Description |
|---|---|
| `FUMApreproessing&lifeover.R` | Preprocessing and liftover script for preparing GWAS results for FUMA |
| `GWAS_regenie` | REGENIE-based GWAS analysis scripts |
| `Meta-analysis` | Scripts for GWAS meta-analysis across datasets/cohorts |
| `prune` | SNP pruning and LD filtering scripts for GSA-53K and MEG-12K datasets |

---

## Main Analyses


### 1. SNP Pruning
LD pruning was performed to:
- Remove highly correlated variants
- Generate approximately independent SNP sets
- Support downstream analyses including PCA and clumping

### 2. GWAS Using REGENIE
Genome-wide association analyses were performed using REGENIE.

Pipeline includes:
- Step 1 polygenic prediction model fitting
- Step 2 association testing
- Covariate adjustment
- Chromosome-wise analysis

---

### 3. Meta-analysis
GWAS summary statistics from multiple datasets were combined using meta-analysis approaches.

Outputs include:
- Combined effect sizes
- Meta-analysis P-values
- Harmonized summary statistics

---

### 4. FUMA Preprocessing and Liftover
Scripts prepare GWAS summary statistics for:
- FUMA annotation
- hg19/hg38 coordinate conversion
- SNP formatting and filtering

---

## Software Requirements

### Main Tools
- REGENIE
- PLINK
- R

### Common R Packages
```r
library(data.table)
library(dplyr)
