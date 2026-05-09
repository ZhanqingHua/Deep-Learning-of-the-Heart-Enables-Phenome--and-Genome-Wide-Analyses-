# Partial Correlation Analysis

## Overview

This R script performs partial correlation analysis between radiomic features and clinical biomarkers, adjusting for height.

The pipeline:
1. Merges phenotype/radiomic data with average height data
2. Calculates partial correlations between radiomic traits and biomarkers
3. Adjusts correlations for height
4. Generates a heatmap of Pearson partial correlations

---

## Input Files

| File | Description |
|---|---|
| `heart_merged_data.csv` | Heart phenotype/radiomic dataset |
| `Height_EMPI_avg.csv` | Average height data matched by EMPI |
| `aortic_Height_merged.tsv` | Aortic phenotype/radiomic dataset merged with height |

---

## Main Steps

### 1. Merge Height Data

The script merges the main dataset with average height by:

```r
EMPI
