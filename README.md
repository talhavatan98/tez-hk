# IDH Mutation Detection via Multi-Omics Data Integration (GBMLGG)

Machine learning pipeline for predicting IDH1/IDH2 mutation status in 
glioma patients (TCGA GBMLGG cohort) by integrating transcriptomic 
(mRNA-seq) and epigenomic (DNA methylation 450k) data.

## Overview

This repository contains the full analysis pipeline comparing single-omics 
and multi-omics machine learning approaches for IDH mutation classification:

- **sPLS-DA** (mRNA only)
- **sPLS-DA** (Methylation only)
- **DIABLO** (multi-omics integration, mixOmics)
- **Concatenated Random Forest** (mRNA + Methylation features combined)

Model performance is evaluated using ROC-AUC, PR-AUC, F1, MCC, bootstrap 
95% confidence intervals, and DeLong pairwise tests. Kaplan-Meier survival 
analysis is included for clinical validation.

## Data Source

All data is retrieved programmatically from **UCSC Xena** (TCGA GBMLGG cohort):
- Gene expression: `HiSeqV2` (IlluminaHiSeq RNA-seq, log2 RSEM-normalized)
- DNA methylation: `HumanMethylation450` (Illumina 450k array, beta values)
- Somatic mutations: `mc3_gene_level` (used to define IDH1/IDH2 mutation status)
- Clinical data: `GBMLGG_clinicalMatrix`

No local data files are required — the script downloads data directly 
into memory at runtime.

## Environment

- **Platform:** Google Colab (R kernel)
- **Language:** R
- **Key packages:** `mixOmics`, `sva`, `impute`, `randomForest`, `caret`, 
  `pROC`, `PRROC`, `survival`, `survminer`, 
  `IlluminaHumanMethylation450kanno.ilmn12.hg19` (CRAN + Bioconductor)

All required packages are automatically installed at the start of the script.

## How to Run

1. Open the script in Google Colab with an R runtime.
2. Run all cells sequentially (the script handles package installation, 
   data download, preprocessing, modeling, and evaluation automatically).
3. Set `FAST_MODE <- TRUE` for quick validation (lower computational cost) 
   or `FAST_MODE <- FALSE` for full analysis with extended cross-validation 
   and permutation testing.
4. All outputs (figures, tables, selected feature lists) are saved to 
   `./IDH_Pipeline_Results/`.

## Outputs

- Model comparison tables (ROC-AUC, PR-AUC, F1, MCC)
- ROC curves (individual and combined)
- Confusion matrices
- Bootstrap 95% CI and DeLong pairwise significance tests
- Selected gene/probe lists per model (Enrichr/g:Profiler-ready)
- Demographic summary table and Kaplan-Meier survival curves

## Reproducibility

- Fixed random seed (`SEED_VAL = 42`) across all modeling steps
- Train/test split performed before scaling to prevent data leakage
- Permutation testing (label shuffling) for empirical p-values

## License

This project is licensed under the MIT License.

## Citation

If you use this pipeline, please cite the associated publication 
(citation details to be added upon publication).