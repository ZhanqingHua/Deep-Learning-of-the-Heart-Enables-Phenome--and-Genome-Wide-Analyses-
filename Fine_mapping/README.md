# Fine Mapping Pipeline

## Overview 

This script performs statistical fine-mapping using SuSiE (`susieR`) on GWAS summary statistics within a target genomic region.

The pipeline:
1. Loads GWAS summary statistics
2. Extracts SNPs within ±500 kb of a target locus
3. Matches SNPs with genotype BIM files
4. Computes LD matrices using PLINK
5. Runs SuSiE fine-mapping using summary statistics and LD information

---

## Main Steps

### 1. Load GWAS Summary Statistics
- Import GWAS meta-analysis results
- Parse chromosome and position information from `MarkerName`
- Generate harmonized SNP IDs

---

### 2. Extract Regional SNPs
Select SNPs within:
```text
target position ±100 kb
