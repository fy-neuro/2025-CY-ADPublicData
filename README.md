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
│   └── ECM_related_genes/ # ECM gene lists (V1-V4 versions)
├── task/                  # Analysis tasks (organized by date/dataset)
│   ├── YYYYMMDD-Description/
│   │   ├── code/         # Task-specific R scripts/RMarkdown
│   │   ├── plot/         # Generated plots (excluded from Git)
│   │   └── Documents/    # Additional documentation/outputs (excluded)
└── report/               # Analysis reports and organized results
    ├── 20260108/         # ECM V4 analysis results (Jan 2026)
    │   ├── csv/         # Wilcoxauc data files
    │   ├── figures/     # Wilcoxauc visualizations
    │   ├── docs/        # Analysis documentation
    │   ├── pseudobulk/  # Pseudobulk results (task18, task19)
    │   └── *.md         # Navigation guides and summaries
    └── YYYYMMDD/        # Other dated reports
```

**Note**: Large data files (`.rds`, `.h5`, `.csv.gz`, plots, documents) are excluded via `.gitignore` to keep repository size manageable.

## Key Analyses

### 1. ECM Communication V4 Gene Analysis (January 2026)
**Latest comprehensive analysis using two complementary methods**

- **Gene List**: ECM Communication V4 genes (20 genes involved in ECM-cell communication)
- **Datasets**:
  - GSE174367 (Human snRNA-seq): 7 cell types
  - GSE188545 (Human MTG): 8 cell types
  - GSE138852 (Mouse): Excluded from wilcoxauc (SCT data incompatibility)

**Methods Applied**:

1. **Wilcoxauc AUC Analysis** (Single-cell level)
   - Tool: `presto::wilcoxauc`
   - Output: AUC scores (0.5-1.0) for each gene-celltype combination
   - Results: Mean AUC = 0.5 (ECM V4 genes stable in AD)
   - Total combinations analyzed: 600 (15 celltypes × 20 genes)

2. **Pseudobulk DESeq2 Analysis** (Sample-level)
   - Method: Aggregate by sample → DESeq2 differential expression
   - Output: log2 fold change, adjusted p-values
   - Volcano plots for each cell type
   - Heatmaps of log2FC across genes and celltypes

**Results Location**: `report/20260108/`
- Navigation guides: README.md, QUICK_REFERENCE.md
- Wilcoxauc results: `csv/`, `figures/`, `docs/`
- Pseudobulk results: `pseudobulk/task18_GSE188545/`, `pseudobulk/task19_3datasets/`
- Complete file index: `COMPLETE_FILE_INDEX.txt`

**Key Findings**:
- ECM Communication V4 genes show high stability in AD (mean AUC = 0.5)
- Only 2/600 combinations showed AUC > 0.6 (GSE188545)
- No genes showed AUC > 0.7 (strong effect)
- Pseudobulk analysis provides more conservative, sample-level validation

### 2. ECM Gene V1-V3 Analysis (December 2025)
- Previous versions of ECM gene lists analyzed
- Methods: Wilcoxauc AUC analysis
- Datasets: GSE174367, GSE188545
- Results available in respective task directories

### 3. GSE188545 Human AD Analysis
- **Dataset**: 62,263 cells from human middle temporal gyrus (AD vs Healthy Control)
- **Cell Types**: 8 major CNS types (Glutamatergic/GABAergic neurons, Astrocytes, Microglia, Oligodendrocytes, OPC, Endothelial, Unknown)
- **Analyses Completed**:
  - Data loading and quality control
  - Harmony batch correction
  - Cell type annotation using marker genes
  - Differential expression analysis (AD vs HC per cell type)
  - ECM gene AUC analysis across cell types (V1-V4)
- **Output**: Annotated Seurat objects, UMAP plots, volcano plots, AUC heatmaps
- **Metadata Variables**: `condition` (AD/Control), `celltype`, `sample`

### 4. GSE174367 Human Multiome Analysis
- **Dataset**: snRNA-seq + snATAC-seq from human AD brain
- **Analyses Completed**:
  - Data integration of RNA and ATAC modalities
  - AUC analysis of ECM-related genes (V1-V4: communication, structure, all ECM)
  - Cell type-specific differential expression
  - Pseudobulk aggregation by sample
- **Output**: AUC tables, heatmap visualizations, DEG results
- **Metadata Variables**: `Diagnosis` (condition), `Cell.Type` (celltype), `SampleID` (sample)

### 5. GSE138852 Mouse AD Analysis
- **Dataset**: Mouse Alzheimer's model single-cell data
- **Analyses Completed**:
  - Basic data loading and preprocessing
  - SCTransform normalization
  - Initial visualization and quality control
  - Note: Uses SCT-transformed data, incompatible with wilcoxauc (requires log-normalized data)
- **Metadata Variables**: `oupSample.batchCond` (condition), `oupSample.cellType` (celltype)

## Installation and Setup

### R Package Dependencies
```r
# Core packages
install.packages(c("Seurat", "dplyr", "ggplot2", "here", "harmony", "presto"))

# Differential expression (pseudobulk)
BiocManager::install("DESeq2")

# Enrichment analysis
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db"))

# Visualization
install.packages(c("pheatmap", "viridis", "patchwork", "ComplexHeatmap"))
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

## Analysis Methods

### Wilcoxauc AUC Analysis (Single-cell level)
**Purpose**: Detect genes that discriminate between conditions (e.g., AD vs Control) at single-cell resolution

**Method**: Wilcoxon rank sum test with Area Under Curve (AUC) calculation
```r
library(presto)

# Run wilcoxauc
auc_results <- presto::wilcoxauc(
  seurat_object,
  group_by = "condition",
  assay = "data"
)
```

**Output Interpretation**:
- **AUC = 0.5**: No discriminatory power (genes expressed similarly in both conditions)
- **AUC > 0.6**: Moderate effect
- **AUC > 0.7**: Strong effect
- **AUC > 0.8**: Very strong effect

**Advantages**:
- Fast computation on single-cell level
- Sensitive to subtle changes
- No aggregation required

**Limitations**:
- May detect cell-level heterogeneity not consistent at sample level
- Higher false positive rate
- Requires log-normalized data (not compatible with SCTransform)

### Pseudobulk DESeq2 Analysis (Sample-level)
**Purpose**: Validate differential expression at sample level with more conservative statistics

**Method**: Aggregate cells by sample → Bulk RNA-seq differential expression
```r
library(Seurat)
library(DESeq2)

# Step 1: Aggregate by sample
sobj_agg <- AggregateExpression(
  sobj,
  group.by = "sample",
  assays = "RNA"
)

# Step 2: Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "AD", "Control"))
```

**Output Interpretation**:
- **log2FoldChange**: Magnitude of expression change (positive = up in AD, negative = down)
- **padj**: Adjusted p-value (FDR-controlled)
- **Significance**: Typically padj < 0.05 AND |log2FC| > 1

**Advantages**:
- Controls for sample-level variability
- More conservative, fewer false positives
- Compatible with count-based data
- Standard method in bulk RNA-seq

**Limitations**:
- Requires sufficient samples per group (≥3 recommended)
- Less sensitive to subtle changes
- Loses single-cell resolution

### Choosing the Right Method

| Use Case | Recommended Method | Rationale |
|----------|-------------------|-----------|
| Exploratory screening | Wilcoxauc | High sensitivity, detects subtle patterns |
| Publication validation | Pseudobulk DESeq2 | More conservative, sample-level rigor |
| Limited samples (<3 per group) | Wilcoxauc | DESeq2 underpowered with few samples |
| SCTransform data | Pseudobulk | Wilcoxauc incompatible with SCT assays |
| Final high-confidence genes | **Both** | Genes significant in both methods |

### Combined Analysis Workflow
1. **Screen with wilcoxauc** → Identify candidate genes
2. **Validate with pseudobulk** → Filter to high-confidence genes
3. **Compare results** → Genes significant in both methods are most reliable

## File Formats

- **Count matrices**: `.csv.gz` (compressed CSV)
- **Seurat objects**: `.rds` (R serialized format)
- **10X data**: `.h5` (H5 files) or sample directories
- **Gene lists**: `.xlsx` (Excel files)
- **Plots**: `.pdf`, `.png`, `.eps`

## Accessing Latest Results

The most recent comprehensive analysis (ECM Communication V4, January 2026) is available in `report/20260108/`:

**Quick Start**:
```bash
# Navigate to results directory
cd report/20260108/

# Read navigation guide
cat README.md

# View quick reference (30-second summary)
cat QUICK_REFERENCE.md

# Browse complete file catalog
cat COMPLETE_FILE_INDEX.txt
```

**Key Result Files**:
- **Wilcoxauc summary**: `csv/summary_statistics_per_dataset.csv`
- **Wilcoxauc heatmap**: `figures/All_Datasets_ECM_communicationV4_AUC_heatmap.pdf`
- **Pseudobulk volcanoes**: `pseudobulk/task18_GSE188545/volcano_plots/`
- **Pseudobulk combined results**: `pseudobulk/task19_3datasets/data/All_Datasets_ECM_communicationV4_pseudobulk_DEG_combined.csv`
- **Daily summary (Chinese)**: `20260108_task_summary.md`

**For detailed methodology**:
- Pseudobulk guide: `PSEUDOBULK_RESULTS.md`
- Wilcoxauc file catalog: `FILE_INDEX.txt`

## Recent Updates

### January 2026
- **2026-01-08**: ECM Communication V4 gene analysis completed using two methods:
  - Wilcoxauc AUC analysis across 2 human datasets (600 gene-celltype combinations)
  - Pseudobulk DESeq2 analysis with volcano plots and heatmaps
  - Results organized in `report/20260108/` with comprehensive navigation guides
  - Key finding: ECM V4 genes are highly stable in AD (mean AUC = 0.5)
- **2026-01-08**: Created 4 analysis tasks (tasks 16-19) for ECM V4 analysis
- **2026-01-08**: Fixed multiple technical issues (JoinLayers compatibility, presto column detection, SCT data handling)

### December 2025
- **2025-12-14**: GSE188545 cell type annotation and differential expression completed
- **2025-12-14**: GSE188545 ECM gene AUC analysis completed
- **2025-12-12**: GSE174367 AUC analysis for ECM genes completed
- **2025-12-09**: GSE188545 data loading, QC, and integration completed

## License

This project contains analysis code for public datasets. The code is provided for research purposes. Dataset usage should comply with GEO data use agreements.

## Contact

For questions about the analysis code, please refer to the project documentation in `AGENTS.md`.

---

*Last updated: January 8, 2026*