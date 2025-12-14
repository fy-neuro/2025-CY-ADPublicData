# AD Public Single-Cell Data Analysis

This repository contains R-based analyses of Alzheimer's disease (AD) single-cell RNA-seq datasets from public repositories (GEO). The project focuses on cell type annotation, differential expression analysis, and extracellular matrix (ECM) gene expression profiling in AD.

## Project Overview

- **Domain**: Alzheimer's disease research using public single-cell datasets
- **Primary Analysis**: Single-cell RNA-seq analysis, cell type annotation, differential expression, AUC analysis of ECM genes
- **Key Datasets**: GSE138852 (mouse), GSE174367 (human snRNA-seq + snATAC-seq), GSE188545 (human MTG)
- **Tools**: Seurat v5, Harmony, presto, ggplot2, dplyr, clusterProfiler
- **Language**: R with Chinese comments, English variable/function names

## Repository Structure

```
├── src/                    # Reusable R functions and utilities
│   ├── 1stAnnotation.R    # Cell type marker definitions
│   ├── getdir.R           # Directory path utilities
│   ├── recluster.R        # Reclustering functions
│   └── plot_function/     # Plotting utilities
├── data/                  # Raw and processed datasets (excluded from Git)
│   ├── GSE138852/        # Mouse AD single-cell data
│   ├── GSE174367/        # Human AD snRNA-seq + snATAC-seq multiome
│   ├── GSE188545/        # Human AD middle temporal gyrus data
│   └── ECM_related_genes/ # ECM gene lists
├── task/                  # Analysis tasks (organized by date/dataset)
│   ├── YYYYMMDD-Description/
│   │   ├── code/         # Task-specific R scripts/RMarkdown
│   │   ├── plot/         # Generated plots (excluded from Git)
│   │   └── Documents/    # Additional documentation/outputs (excluded)
└── report/               # Reports and presentations
```

**Note**: Large data files (`.rds`, `.h5`, `.csv.gz`, plots, documents) are excluded via `.gitignore` to keep repository size manageable.

## Key Analyses

### 1. GSE188545 Human AD Analysis (Latest)
- **Dataset**: 62,263 cells from human middle temporal gyrus (AD vs Healthy Control)
- **Analyses Completed**:
  - Cell type annotation (8 major CNS types)
  - Differential expression analysis (AD vs HC per cell type)
  - ECM gene AUC analysis across cell types
- **Output**: Annotated Seurat objects, volcano plots, AUC heatmaps

### 2. GSE174367 Human Multiome Analysis
- **Dataset**: snRNA-seq + snATAC-seq from human AD brain
- **Analyses Completed**:
  - AUC analysis of ECM-related genes (communication, structure, all ECM)
  - Cell type-specific differential expression
- **Output**: AUC tables, heatmap visualizations

### 3. GSE138852 Mouse AD Analysis
- **Dataset**: Mouse Alzheimer's model single-cell data
- **Analyses Completed**:
  - Basic data loading and preprocessing
  - Initial visualization and quality control

## Installation and Setup

### R Package Dependencies
```r
install.packages(c("Seurat", "dplyr", "ggplot2", "here", "harmony", "presto"))
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db"))
```

### Project Setup
```r
# Set working directory to project root
project_root <- "/path/to/20251120-CY-ADpublic_data"
setwd(project_root)

# Source utility functions
source(here::here("src/getdir.R"))
```

## Usage Examples

### Running Analysis Scripts
```bash
# From project root directory
Rscript task/20251214-10-GSE188545_Annotation/code/fast_annotation.R
```

### Knitting RMarkdown Documents
```r
rmarkdown::render("task/20251212-6-GSE174367_Load_Data+UMAP+Markerplot/code/1.0Load Data.Rmd")
```

## Code Conventions

- **Function/variable names**: snake_case (e.g., `plot_proportion_dodge`)
- **Indentation**: 2 spaces (tidyverse style)
- **Comments**: Chinese with `#` for single-line, `##` for section headers
- **Path handling**: Use `here::here()` for project-relative paths

## File Formats

- **Count matrices**: `.csv.gz` (compressed CSV)
- **Seurat objects**: `.rds` (R serialized format)
- **10X data**: `.h5` (H5 files) or sample directories
- **Gene lists**: `.xlsx` (Excel files)
- **Plots**: `.pdf`, `.png`, `.eps`

## Recent Updates (December 2025)

- **2025-12-14**: GSE188545 cell type annotation and differential expression completed
- **2025-12-14**: GSE188545 ECM gene AUC analysis completed
- **2025-12-12**: GSE174367 AUC analysis for ECM genes completed
- **2025-12-09**: GSE188545 data loading, QC, and integration completed

## License

This project contains analysis code for public datasets. The code is provided for research purposes. Dataset usage should comply with GEO data use agreements.

## Contact

For questions about the analysis code, please refer to the project documentation in `AGENTS.md`.

---

*Last updated: December 14, 2025*