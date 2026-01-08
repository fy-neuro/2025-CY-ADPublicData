# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R-based single-cell RNA-seq analysis project for Alzheimer's disease (AD) research using public datasets from GEO. Primary focus: cell type annotation, differential expression analysis, and extracellular matrix (ECM) gene expression profiling.

**Key datasets**: GSE138852 (mouse AD), GSE174367 (human snRNA-seq + snATAC-seq multiome), GSE188545 (human MTG, 62K cells).

## Running R Scripts

From project root directory:
```bash
Rscript task/20251214-10-GSE188545_Annotation/code/fast_annotation.R
```

Interactive R session:
```r
setwd("/media/bigbaby/Zhang-1号/Projects/Others/20251120-CY-ADpublic_data")
```

Knit RMarkdown:
```r
rmarkdown::render("task/20251212-6-GSE174367_Load_Data+UMAP+Markerplot/code/1.0Load Data.Rmd")
```

## Code Architecture

### Task-Based Workflow
- Each analysis isolated in `task/YYYYMMDD-N-Description/` directories
- Subdirectories: `code/` (R scripts/RMarkdown), `plot/` (outputs), `Documents/` (results)
- Tasks share data from `data/` directory and utilities from `src/`

### Reusable Functions (`src/`)

**Core utilities**:
- `src/getdir.R`: Path handling using `here::here()` pattern for project-relative paths
- `src/1stAnnotation.R`: Cell type marker gene definitions (mouse/human) - defines `Markers_CNS` for CNS cell types (glutamatergic/GABAergic neurons, OPC, oligodendrocytes, astrocytes, T cells, microglia, vascular cells)
- `src/recluster.R`: `recluster()` function - performs SCTransform normalization, PCA, Harmony batch correction (using `orig.ident`), clustering, and UMAP

**Plotting functions** (`src/plot_function/`):
- `Dotplot.R`: `dotplot_fy()` - cell type marker dot plots with cols4all colors
- `volcano.R`: Volcano plots for DE results
- `Vlnplot.R`: Violin plots
- `Proportion.R`: Cell proportion plots
- `GO Enrichment.R`: GO enrichment analysis (uses mouse `org.Mm.eg.db`)

**Specialized utilities**:
- `src/create_loupe.R`: Convert Seurat objects to Loupe Browser (.cloupe) files

### Data Structure Patterns

**Count matrices**: Compressed CSV (`.csv.gz`) with genes as rows, cells as columns

**Seurat objects**: Saved as `.rds` files, can be multi-layered in Seurat v5

**10X data**: H5 files (`.h5`) or per-sample directories with separate `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` files (e.g., `data/GSE188545/GSE188545_GEM/`)

**ECM gene lists**: Excel files (`.xlsx`) in `data/ECM_related_genes/`:
- `ECM_communication_genesV3.xlsx` (most recent, 49 genes)
- `ECM_communication_genesV2.xlsx`, `ECM_communication_genesV4.xlsx`
- `ECM_structure_genes.xlsx`
- `ECM_Genes_Grouped_Final.xlsx`

**Metadata**: CSV files (`.csv.gz`) for sample covariates

### Path Handling Conventions

Use `here::here()` for project-relative paths:
```r
here("src/getdir.R") %>% source()
```

Many scripts set working directory to project root:
```r
project_root <- file.path(script_dir, "..", "..", "..")
setwd(project_root)
```

**Critical**: Project path contains Chinese characters (`Zhang-1号`) - ensure UTF-8 encoding.

## Seurat v5 Compatibility

**Key differences from v4**:
- Multi-layered data structure: layers include `counts`, `data`, `scale.data`
- **Call `JoinLayers()` before** `FindAllMarkers()` and expression scoring
- Use `layer = "data"` instead of deprecated `slot = "data"` in `GetAssayData()`
- Check layers with `Layers()`
- Memory warning: `JoinLayers()` can increase object size 2-3x

**Load plyr before dplyr** to avoid namespace conflicts:
```r
library(plyr)
library(dplyr)
```

## Standard Analysis Pipeline

1. **Data loading**: Read count matrices → Create Seurat object
2. **Preprocessing**: `SCTransform()` normalization, PCA, Harmony batch correction
3. **Clustering**: `FindNeighbors()`, `FindClusters()`, `RunUMAP()`
4. **Annotation**: Score-based using markers from `src/1stAnnotation.R`
5. **Differential expression**: presto AUC analysis or Seurat `FindMarkers()`
6. **Visualization**: Functions from `src/plot_function/`

### Typical Annotation Workflow

```r
# Score-based annotation (fast method)
markers <- list(
  Glutamatergic_Neurons = c("RALYL","LDB2","NELL2"),
  GABAergic_Neurons = c("GAD1","GAD2","SLC6A1"),
  OPC = c("PDGFRA","VCAN","OLIG1"),
  # ... from src/1stAnnotation.R Markers_CNS
)
# Score each cell type and assign max score
```

## Differential Expression Methods

### Wilcoxon AUC (presto)
Preferred for cell-level DE across conditions:
```r
presto::wilcoxauc(sobj, features = gene_list, group_by = "condition")
```

### Pseudobulk (DESeq2/edgeR)
Challenges encountered with Seurat v5 data structure:
- `AggregateExpression()` replaces underscores with dashes in sample IDs
- Need to create sample metadata mapping between original and modified IDs
- Requires ≥2 samples per condition per cell type
- See `task/20260104-15-ECM_comparison_wilcoxauc_pseudobulk/` for debugging attempts

## Naming Conventions

- **Functions/variables**: snake_case (`plot_proportion_dodge`, `run_pseudobulk_final`)
- **File names**: `DescriptiveName.R` (spaces allowed but discouraged)
- **Directories**: `YYYYMMDD-N-Description` with hyphens
- **Code style**: 2-space indentation, Chinese comments with `#` or `##`

## Important Gotchas

### Working Directory
Many scripts assume running from project root after `setwd(project_root)`. Plotting functions create `../plot` and `../Documents` relative to current working directory.

### File Name Spaces
Some files contain spaces (e.g., `1.0Load Data.R`). Use quotes or escape spaces in shell commands.

### GO Enrichment Database
GO functions in `src/plot_function/GO Enrichment.R` default to mouse database (`org.Mm.eg.db`). For human datasets, switch to `org.Hs.eg.db`.

### Memory Usage
Large datasets (60K+ cells) require 8-16GB RAM. Use sparse matrices and Seurat parameters conservatively.

### Package Versions
No formal package management (no renv/packrat). Critical versions observed: R 4.5.2, Seurat 5.3.1.

## Core Library Stack

**Single-cell**: Seurat, harmony (batch correction), presto (fast Wilcoxon AUC)
**Data manipulation**: dplyr, tidyr, tibble, plyr (load before dplyr)
**Visualization**: ggplot2, cols4all, patchwork, ggpubr, pheatmap, viridis
**Annotation**: AnnotationDbi, org.Mm.eg.db, org.Hs.eg.db
**Enrichment**: clusterProfiler, enrichplot, stringr, purrr
**Data I/O**: openxlsx (Excel), Matrix (sparse matrices)

## Current Technical Challenges

### Pseudobulk Analysis Compatibility
Seurat v5's layered data structure and `AggregateExpression()` ID transformation cause issues with DESeq2/edgeR metadata alignment. Manual pseudobulk aggregation implementation in progress (see debug scripts in `task/20260104-15-ECM_comparison_wilcoxauc_pseudobulk/code/`).

### Large-Scale Visualization
60K+ cells require efficient visualization strategies; consider downsampling for interactive plots.
