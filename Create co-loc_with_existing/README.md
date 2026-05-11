Cardiac Radiomics GWAS Project

This repository contains the analysis pipeline and supporting files for the manuscript:

“Integrating Deep Learning-Based Radiomics and Genomics to Define the Genetic Architecture of Cardiac Morphology and its Link to Cardiovascular Disease”

Repository Structure
overlap&clumping&annotation/

Contains scripts and intermediate files used for:

SNP overlap analysis between radiomic GWAS and published cardiovascular disease (CVD) GWAS
LD clumping of significant variants
SNP-to-gene annotation
nearest gene mapping
downstream functional interpretation
Co-loc/

Contains scripts and workflows for colocalization analysis between radiomic GWAS loci and cardiovascular disease GWAS signals.

This section is currently under active development and will be updated following ongoing analyses and discussions.

README.md

Overview of repository structure and analysis workflow.

Analysis Overview

The overall workflow includes:

Deep learning–based extraction of cardiac radiomic features from CT imaging
Genome-wide association studies (GWAS) of imaging-derived phenotypes
Meta-analysis across cohorts
Overlap analysis with published cardiovascular disease GWAS
LD clumping and independent lead SNP identification
Functional annotation and nearest-gene mapping
Colocalization analysis
Downstream biological interpretation
Software and Tools

Major tools used in this project include:

R
PLINK
METAL
SuSiE
FUMA
PoPS

Key R packages include:

data.table
tidyverse
ggplot2
Notes
Some analyses require access to protected genotype and phenotype data.
Intermediate files may not be included due to size limitations and data-sharing restrictions.
Additional scripts and documentation will be continuously updated.
