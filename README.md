# HCA-scRNAseq-RP-genes
This repository contains a collection of R scripts used for preprocessing, cell-type-specific merging, and downstream analysis of single-cell RNA sequencing (scRNA-seq) data across multiple tissues and cell types.

Repository Structure
| File                          | Description                                                                                 |
|-------------------------------|---------------------------------------------------------------------------------------------|
| `Data_input.R`                | Loads and preprocesses raw data for downstream analyses.                                    |
| `seurat_preprocessing.R`      | Preprocessing pipeline using Seurat (QC, normalisation, integration) to generate RDS files  |
| `Analysis_final.R`            | Main script for final data analysis of tissues.                                             |
| `Multitissue_var_final.R`     | Final analysis of gene expression mean and variability across multiple tissues.             |
| `GSE159929_DEG_volcano.R`     | Generates volcano plots for DEG analysis of tissues from GSE159929 dataset.                 |
| `testis_pseudotime.R`         | Pseudotime analysis for testis cell populations.                                            |
| `ery_pseudotime.R`            | Pseudotime analysis for erythroid lineage.                                                  |
| `DPGP_input.R`                | Prepares input data for DPGP (Dirichlet Process Gaussian Process) analysis.                 |
| `DEG.R`                       | Differential expression gene (DEG) analysis for male infertility using DEseq2.              |
| `B_merging.R`                 | Merging and processing of B cell data.                                                      |
| `CD8T_NKG7_merging.R`         | Merging of CD8 T cells expressing NKG7.                                                     |
| `CD8T_merging.R`              | Merging and processing of CD8 T cell data.                                                  |
| `Epithelial_merging.R`        | Merging and processing of epithelial cell data.                                             |
| `Fibro_merging.R`             | Merging and processing of fibroblast cell data.                                             |
| `NKT_merging.R`               | Merging and processing of NKT cell data.                                                    |
| `Plasma_merging.R`            | Merging and processing of plasma cell data.                                                 |
| `SMC_merging.R`               | Merging and processing of smooth muscle cell (SMC) data.                                    |
| `T_merged.R`                  | Merging and processing of T cell data.                                                      |
| `Tgd_merging.R`               | Merging and processing of gamma delta T cell data.                                          |
| `mono_merging.R`              | Merging and processing of monocyte data.                                                    |
| `celltype_analysis.R`         | Analysis of common cell types obtained from merging tissues                                 |
| `LICENSE`                     | License information.                                                                        |
| `README.md`                   | This file.                                                                                  |


🧾 Miscellaneous
README.md: This file.

Usage
Each script is written in R and assumes input data formatted for use with the Seurat package. 

Requirements
R (>= 4.3.1)
Seurat V4.3.3
Packages: Seurat, dplyr, ggplot2, Matrix, EnhancedVolcano, DPGP, and others, depending on the script.


### Usage

1. **Clone the repository:**
   
2. **Prepare your data:**  
   Place your raw or preprocessed data files in the appropriate directory as expected by `Data_input.R`.

3. **Run preprocessing:**  
   Start with  `Data_input.R` and `seurat_preprocessing.R` to prepare your data.

4. **Downstream analyses:**  
   - basic analysis: `Analysis_final.R` 
   - DEG analysis: `DEG.R`, `GSE159929_DEG_volcano.R`
   - Pseudotime: `ery_pseudotime.R`, `testis_pseudotime.R`
   - Intertissue variability: `Multitissue_var_final.R`

5. **Cell type analysis and merging:**  
   Use the relevant `*_merging.R` and `celltype_analysis.R` scripts for your cell types of interest.







