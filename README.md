# Deep Learning of the Heart Enables Phenome- and Genome-Wide Analyses

## Overview

This repository contains the analysis pipelines, scripts, and visualization code used for the manuscript:

**“Deep Learning of the Heart Enables Phenome- and Genome-Wide Analyses”**

The project integrates deep learning–derived cardiac radiomic phenotypes with large-scale genomic and phenomic analyses to investigate the genetic architecture of cardiac morphology and its relationship with cardiovascular disease.

---

## Main Components

### 1. Phenotype Processing
Scripts for:
- Cleaning and preprocessing MGBB CT cohort data
- Removing duplicate scans
- Generating heart and aortic phenotype datasets
- Matching imaging data with biobank identifiers

---

### 2. Radiomic Feature Analysis
Analyses include:
- Correlation analysis with clinical biomarkers
- Partial correlation adjusted for height
- Heatmap visualization of radiomic-biomarker associations

---

### 3. Genome-Wide Association Studies (GWAS)
GWAS analyses were performed using:
- REGENIE
- Meta-analysis pipelines
- SNP pruning and LD filtering

The repository includes scripts for:
- GWAS preprocessing
- Meta-analysis
- FUMA formatting and liftover
- Independent SNP identification

---

### 4. Fine-Mapping and Functional Annotation
Downstream genomic analyses include:
- LD clumping
- SuSiE fine-mapping
- Colocalization analysis
- Nearest gene mapping
- Functional prioritization

---

### 5. Visualization
Scripts for generating:
- Manhattan plots
- Circular/circos plots
- Correlation heatmaps
- Summary tables

---

## Software Requirements

### Main Tools
- R
- PLINK
- REGENIE
- bash utilities
