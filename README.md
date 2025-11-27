# Candidate-SNP GWAS Pipeline for Age-related Macular Degeneration (AMD)

## Overview
This repository contains the R source code utilized for the genetic association analysis presented in the manuscript regarding **AMD-associated genetic variants in the Jeju cohort**.

The pipeline performs a comprehensive analysis of candidate SNPs, including quality control (QC), logistic regression with Firth's penalization fallback, subgroup comparisons (Drusen-driven vs. Pachychoroid-driven AMD), and visualization of results.

> **Note**: Due to privacy regulations and IRB requirements regarding human subject data, the raw genotype and phenotype datasets (`*.csv`) are **not** included in this repository.

## Features
* **Quality Control**:
    * SNP Call Rate filtering (≥ 0.98)
    * Minor Allele Frequency (MAF) filtering (Controls ≥ 0.01)
    * Hardy-Weinberg Equilibrium (HWE) filtering (Controls p ≥ 1e-6)
    * Biallelic SNP checks
* **Association Analysis**:
    * Logistic regression models adjusted for **Age, Sex, and Smoking (PAQ8)**.
    * Automatic fallback to **Firth's penalized likelihood** (`logistf`) when standard GLM fails (e.g., due to complete separation).
    * Correction for multiple testing (Benjamini-Hochberg FDR & Bonferroni).
* **Subgroup Analysis**:
    * Comparisons between phenotypic subgroups: **Drusen-driven** vs. **Pachychoroid-driven** AMD.
* **Visualization**:
    * Manhattan-like plots for candidate regions.
    * Volcano plots and QQ plots.
    * MAF vs. |Delta MAF| scatter plots.

## Project Structure
```text
.
├── R/
│   └── utils_gwas.R        # Helper functions (Genotype parsing, HWE, MAF calculation)
├── scripts/
│   └── run_analysis.R      # Main execution script
├── AMD_GWAS_out/           # Output directory (Results, Tables, Figures)
├── README.md               # Project documentation