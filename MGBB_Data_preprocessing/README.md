# MGBB CT Full Cohort Data Cleaning

## Overview 

This R Markdown script cleans the MGBB CT full cohort phenotype dataset and prepares heart and aortic phenotype files for downstream GWAS analysis.

The pipeline:
1. Loads the full MGBB CT cohort phenotype data
2. Removes repeated scans by keeping the earliest scan per patient
3. Matches EMPI to Biobank Subject ID
4. Creates phenotype datasets for heart and aortic traits
5. Saves cleaned phenotype files in CSV and TSV formats

---

## Input Files

| File | Description |
|---|---|
| `MGBB_CT_Full_Cohort_03032025.csv` | Full CT cohort phenotype dataset |
| `matched_linker_with_full.csv` | Linker file matching EMPI to Biobank Subject ID |

---

## Main Steps

### 1. Remove Repeated CT Scans

For patients with multiple CT scans, only the earliest scan is retained based on:

```r
Image_Date_Numeric
