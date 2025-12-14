# Agent Guide for AD Public Single-Cell Data Analysis Project

This document provides essential information for AI agents working in this R-based single-cell RNA-seq analysis project.

## Project Overview

- **Type**: R project (RStudio) for single-cell RNA-seq data analysis
- **Domain**: Alzheimer's disease (AD) research, public single-cell datasets
- **Primary Tools**: Seurat, Harmony, ggplot2, dplyr, here package
- **Language**: R with Chinese comments (code uses English variable/function names)

## Directory Structure

```
├── 20251120-CY-ADpublic_data.Rproj    # RStudio project file
├── src/                               # Reusable R functions and utilities
│   ├── 1stAnnotation.R               # Cell type marker definitions
│   ├── getdir.R                      # Directory path utilities
│   ├── recluster.R                   # Reclustering functions
│   └── plot_function/                # Plotting utilities
│       ├── Dotplot.R
│       ├── GO Enrichment.R
│       ├── Proportion.R
│       ├── Vlnplot.R
│       └── volcano.R
├── data/                              # Raw and processed datasets
│   ├── GSE138852/                    # Dataset-specific files
│   ├── GSE174367/
│   ├── GSE188545/
│   └── ECM_related_genes/            # ECM gene lists
├── task/                              # Analysis tasks (organized by date/dataset)
│   ├── YYYYMMDD-Description/         # Each task has its own directory
│   │   ├── code/                     # Task-specific R scripts/RMarkdown
│   │   ├── plot/                     # Generated plots
│   │   └── Documents/                # Additional documentation/outputs
└── report/                            # Reports and presentations
```

## Data Formats

- **Count matrices**: Compressed CSV (`.csv.gz`) with genes as rows, cells as columns
- **Seurat objects**: Saved as `.rds` files (R serialized format)
- **10X Genomics data**: H5 files (`.h5`) for filtered feature-bc matrices
- **10X sample directories**: Some datasets contain per-sample directories with separate barcodes/genes/matrix files (e.g., `sample_dirs/` in task directories)
- **Gene lists**: Excel files (`.xlsx`) for ECM-related genes
- **Metadata**: CSV files (`.csv.gz`) for sample covariates

## Essential Commands

### Running R Scripts
```bash
# Execute an R script from project root
Rscript task/20251126-1-GSE138852/code/1.0Load\ Data.R

# Execute with relative path (from within project root)
Rscript src/recluster.R
```

### Knitting RMarkdown Documents
```bash
# Render RMarkdown to HTML
Rscript -e "rmarkdown::render('task/20251212-6-GSE174367_Load_Data+UMAP+Markerplot/code/1.0Load Data.Rmd')"
```

### Interactive R Session
```bash
# Start R in project directory
R

# Within R, set working directory to project root
setwd("/media/bigbaby/Zhang-1号/Projects/Others/20251120-CY-ADpublic_data")
```

### Package Management
No formal package management (no renv, packrat). Dependencies are installed globally:
```r
install.packages(c("Seurat", "dplyr", "ggplot2", "here", "harmony"))
```

## Code Organization and Patterns

### Task-Based Workflow
- Each analysis task is isolated in its own directory under `task/`
- Task directories follow naming pattern: `YYYYMMDD-N-Description` (N is optional sequence number)
- Within each task: `code/` contains R scripts/RMarkdown files, `plot/` for output plots, `Documents/` for additional outputs
- Tasks are generally independent but may share data from `data/` directory

### Reusable Functions
- `src/` contains reusable functions for annotation, clustering, plotting
- Functions are sourced using `here()` package pattern:
  ```r
  here("src/getdir.R") %>% source()
  ```

### Path Handling
- Use `here::here()` for project-relative paths
- Some scripts set working directory relative to script location:
  ```r
  project_root <- file.path(script_dir, "..", "..", "..")
  setwd(project_root)
  ```
- Data files are referenced relative to project root (e.g., `data/GSE138852_counts.csv.gz`)

### Directory Creation
- Scripts often create output directories using `dir.create()`
- Typical pattern: `dir.create("../plot", showWarnings = FALSE, recursive = TRUE)`
- Directories created relative to working directory (usually project root)

### Analysis Pipeline Patterns
1. **Data Loading**: Read count matrices (CSV.gz or 10X H5) → Create Seurat object
2. **Preprocessing**: Normalization (SCTransform), PCA, batch correction (Harmony)
3. **Clustering**: Find neighbors, find clusters, UMAP/t-SNE visualization
4. **Annotation**: Use marker genes from `src/1stAnnotation.R`
5. **Plotting**: Use functions from `src/plot_function/`

## Naming Conventions and Style

### Functions and Variables
- **snake_case** for function and variable names (`plot_proportion_dodge`, `perform_GO_analysis`)
- Descriptive names in English despite Chinese comments

### Code Style
- **Indentation**: 2 spaces (tidyverse style)
- **Line length**: No strict limit observed
- **Spacing**: Spaces around operators (`<-`, `+`, `=`), after commas
- **Comments**: Chinese comments with `#` for single-line, `##` for section headers
- **Function definitions**: No space before parentheses, opening brace on same line:
  ```r
  function_name <- function(arg1, arg2) {
    # body
  }
  ```

### File Names
- R scripts: `DescriptiveName.R` (spaces allowed but not recommended)
- RMarkdown: `Description.Rmd` (spaces allowed)
- Directories: `Descriptive-Name` with hyphens

## Testing

No formal test suite (no testthat, no tests/ directory). Validation is done through:
- Visual inspection of plots
- Checking Seurat object dimensions and metadata
- Manual review of analysis outputs

## Important Gotchas

### Working Directory Assumptions
- Many scripts assume they are run from project root after setting `setwd(project_root)`
- Some plotting functions create `../plot` and `../Documents` directories relative to current working directory
- Always verify working directory before running scripts

### Path Encoding
- Project path contains Chinese characters (`Zhang-1号`)
- Ensure file paths are properly encoded in R (UTF-8)

### Memory Usage
- Single-cell datasets can be large (GBs of memory)
- Use sparse matrices and appropriate Seurat parameters

### Package Versions
- No version locking; ensure compatible versions of Seurat, Harmony, etc.
- Some functions may depend on specific Seurat API (v3/v4)
- **Current versions observed**: R 4.5.2, Seurat 5.3.1 (check session_info.txt in task directories)

### System Locale
- System uses Chinese locale (`zh_CN.UTF-8`)
- May affect file reading/writing with special characters
- UTF-8 encoding used throughout

### File Name Spaces
- Some file names contain spaces (e.g., `1.0Load Data.R`)
- Escape spaces or use quotes in shell commands

### Undefined Functions in Sourced Scripts
- Some older scripts call `get_script_path()` which is not defined in `src/getdir.R`
- Use `here::here()` pattern for project-relative paths instead
- Newer scripts use `project_root <- here::here(); setwd(project_root)`

### Mouse vs Human Annotation Databases
- GO enrichment functions in `src/plot_function/GO Enrichment.R` use mouse database (`org.Mm.eg.db`)
- For human datasets, ensure appropriate organism database (e.g., `org.Hs.eg.db`)
- Cell type markers in `src/1stAnnotation.R` include both mouse and human genes

## Project-Specific Context

### Key Datasets
- **GSE138852**: Mouse Alzheimer's model single-cell data
- **GSE174367**: Human AD snRNA-seq + snATAC-seq multiome data
- **GSE188545**: Additional AD single-cell dataset

### Analysis Focus
- Cell type annotation using canonical markers
- Reclustering with Harmony for batch correction
- ECM (Extracellular matrix) related gene expression analysis
- Cell proportion analysis across conditions

### Common Libraries

#### Core Libraries (used in most scripts)
- **Seurat**: Single-cell analysis toolkit
- **Harmony**: Batch effect correction
- **ggplot2**: Plotting
- **dplyr/tidyr**: Data manipulation
- **here**: Project-relative paths
- **Matrix**: Sparse matrix support

#### Visualization Libraries
- **cols4all**: Color palettes
- **patchwork**: Plot arrangement
- **ggpubr**: Publication-ready plots
- **ggrepel**: Label placement in plots
- **ggvolcano**: Volcano plot visualization

#### Specialized Analysis Libraries
- **AnnotationDbi**, **org.Mm.eg.db**: Gene annotation (mouse-specific)
- **clusterProfiler**, **enrichplot**: GO enrichment analysis
- **stringr**, **purrr**: String manipulation and functional programming

## Quick Start for New Tasks

1. Create new task directory: `task/YYYYMMDD-N-Description/` (N is optional sequence number, e.g., `20251212-6-GSE174367_Load_Data+UMAP+Markerplot`)
2. Create subdirectories: `code/`, `plot/`, `Documents/`
3. Write R scripts in `code/` using existing patterns:
   - Load libraries (Seurat, dplyr, here)
   - Source `src/getdir.R` for path utilities
   - Set working directory to project root
   - Load data from `data/` directory
   - Use functions from `src/` for common operations
   - Save outputs to task directory
4. Run scripts from project root with `Rscript`

## Recent Task Summary: GSE188545 Cell Type Annotation (2025-12-14)

**Task**: Perform comprehensive cell type annotation for GSE188545 dataset (62,263 cells from human AD middle temporal gyrus).

**Status**: 
- ✅ Basic cell type annotation completed using score-based method (`fast_annotation.R`)
- ✅ Annotated Seurat object saved: `task/20251214-10-GSE188545_Annotation/GSE188545_sobj_annotated_fast.rds` (1.75 GB)
- ✅ Visualizations generated: UMAP colored by cell type, dot plot of key markers
- ✅ Differential expression analysis (AD vs HC) completed using `run_de_analysis.R` (DE results and volcano plots generated)
- ❌ GO enrichment analysis pending due to `clusterProfiler` installation issues

**Key Details**:
- **Cell types assigned**: 8 major CNS types (Oligodendrocytes, Astrocytes, Glutamatergic_Neurons, GABAergic_Neurons, OPC, Microglia, Endothelial, Unknown)
- **Method**: Average expression scoring of canonical marker sets from `src/1stAnnotation.R`
- **Dataset**: GSE188545 human Alzheimer's disease single-cell RNA-seq (62,263 cells, 33 clusters at resolution 0.5)
- **Seurat version**: v5.3.1 (used `JoinLayers()` before marker detection, `layer = "data"`)

**Blockers**: 
- `clusterProfiler` installation requires system dependency `libcairo2-dev` (not installed)
- `FindAllMarkers()` on joined data layers may hang due to memory/timeout

**Next Steps**:
1. Install `libcairo2-dev` system package: `sudo apt-get install libcairo2-dev`
2. Install `clusterProfiler`: `BiocManager::install('clusterProfiler')`
3. Run GO enrichment analysis for significant DE genes using `annotation_pipeline.R` (step 8)
4. Consider reclustering major cell types for subpopulation analysis

**Files Created**:
- `task/20251214-10-GSE188545_Annotation/plot/annotation_fast/`
- `task/20251214-10-GSE188545_Annotation/Documents/annotation_fast/`
- `task/20251214-10-GSE188545_Annotation/Documents/annotation_summary.txt`

## Related Task Summary: GSE174367 ECM Gene AUC Analysis (2025-12-12)

**Task**: Perform AUC (Area Under the ROC Curve) analysis of ECM-related genes across cell types in GSE174367 dataset to identify genes differentially expressed between AD and Control groups.

**Analysis Completed**:
- ✅ AUC analysis for three ECM gene categories: all ECM genes, communication genes, structure genes
- ✅ Generated CSV tables with AUC values per gene per cell subclass (7 subclasses total)
- ✅ Created heatmap visualizations for each gene category
- ✅ Summary document created: `task/20251212-7-GSE174367_AUC/Documents/analysis_summary.md`

**Key Details**:
- **Dataset**: GSE174367 human Alzheimer's disease single-nucleus RNA-seq + ATAC-seq multiome
- **Method**: `presto::wilcoxauc` with comparison `c("AD", "Control")`
- **Cell types**: 7 subclasses (ASC/astrocytes, EX/excitatory neurons, INH/inhibitory neurons, MG/microglia, ODC/oligodendrocytes, OPC/OPC, PER.END/pericytes-endothelial)
- **Gene categories**: All ECM (4,397 gene-subclass pairs), Communication (281 pairs), Structure (435 pairs)

**AUC Interpretation**:
- AUC = 0.5: No discriminatory power (equal expression AD/Control)
- AUC > 0.5: Higher expression in AD group
- AUC < 0.5: Higher expression in Control group

**Output Files**:
- `task/20251212-7-GSE174367_AUC/Documents/` - CSV tables with AUC values
- `task/20251212-7-GSE174367_AUC/plot/` - Heatmap PDF visualizations
- `task/20251212-7-GSE174367_AUC/Documents/analysis_summary.md` - Comprehensive summary

**Note**: Additional file with typo (`Sturcture_gene_list_AUC_per_subclass.csv`) appears to be duplicate of structure gene results.

## New Task Created: GSE188545 AUC Analysis (2025-12-14)

**Task**: `20251214-11-GSE188545_AUC`

**Status**: Analysis completed (December 14, 2025).

**Analysis Completed**:
- ✅ ECM gene AUC analysis for four gene lists: communication, communication V2, all ECM genes, structure genes
- ✅ AUC calculation using `presto::wilcoxauc` with AD vs HC comparison across 8 cell types
- ✅ Generated CSV tables with AUC values per gene per cell type
- ✅ Created heatmap visualizations for each gene category
- ✅ Summary statistics and analysis report generated

**Key Results**:
- **Gene coverage**: 299/314 ECM genes present in dataset
- **Mean AUC**: ~0.497 across all gene lists (close to 0.5, indicating minimal differential expression)
- **Genes with AUC > 0.7**: 3 genes in all ECM list, 1 in structure list
- **Cell types analyzed**: 8 major CNS cell types (excluding Unknown due to insufficient cells)

**Output Files**:
- `task/20251214-11-GSE188545_AUC/Documents/` – CSV tables with AUC results and summary statistics
- `task/20251214-11-GSE188545_AUC/plot/` – Heatmap PDF visualizations
- `task/20251214-11-GSE188545_AUC/Documents/analysis_summary.md` – Comprehensive analysis summary

**Technical Details**:
- **Script**: `code/ecm_auc_analysis.R` (adapted with `JoinLayers()` for Seurat v5 compatibility)
- **Packages**: Seurat v5.3.1, presto, dplyr, tidyr, tibble, pheatmap, viridis
- **Runtime**: ~30 minutes
- **Memory**: Peak ~6 GB

---

*Last updated: December 14, 2025*

*Last updated: December 14, 2025*