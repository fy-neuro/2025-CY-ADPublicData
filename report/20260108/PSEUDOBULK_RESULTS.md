# Pseudobulk DESeq2 Analysis Results

## 📁 Directory Structure

```
pseudobulk/
├── task18_GSE188545/          ← Single dataset pseudobulk analysis
│   ├── volcano_plots/         ← Volcano plots for each cell type (7 plots)
│   ├── heatmaps/              ← Log2FC heatmaps
│   ├── summary_plots/         ← Summary statistics visualizations
│   └── data/                  ← CSV data files (5 files)
│
└── task19_3datasets/          ← Multi-dataset pseudobulk analysis
    ├── heatmaps/              ← Log2FC heatmaps per dataset
    ├── summary_plots/         ← Summary statistics
    └── data/                  ← CSV data files (8 files)
```

---

## 📊 Analysis Overview

### Method
- **Pseudobulk aggregation**: Aggregate expression by sample
- **Differential expression**: DESeq2 (bulk RNA-seq method)
- **Level**: Sample-level (more conservative than single-cell)

### Datasets Analyzed
- **Task 18**: GSE188545 only (8 cell types)
- **Task 19**: GSE174367 + GSE188545 (GSE138852 skipped)

---

## 🔍 Task 18: GSE188545 Pseudobulk Analysis

### Key Files

#### 📈 Visualizations
**Volcano Plots** (`volcano_plots/`):
- `volcano_Astrocytes.pdf`
- `volcano_Endothelial.pdf`
- `volcano_GABAergic_Neurons.pdf`
- `volcano_Glutamatergic_Neurons.pdf`
- `volcano_Microglia.pdf`
- `volcano_Oligodendrocytes.pdf`
- `volcano_OPC.pdf`

Each plot shows:
- **X-axis**: log2 fold change (AD vs Control)
- **Y-axis**: -log10(padj)
- **Red points**: Significant genes (padj < 0.05, |log2FC| > 1)

**Heatmaps** (`heatmaps/`):
- `log2_fold_change_heatmap.pdf` - Log2FC for all genes across cell types

**Summary Plots** (`summary_plots/`):
- `summary_genes_tested.pdf` - Number of genes tested per cell type
- `summary_significant_genes.pdf` - Number of significant genes per cell type

#### 📊 Data Files (`data/`)

| File | Size | Description |
|------|------|-------------|
| `ECM_communicationV4_pseudobulk_DEG_results.csv` | 21 KB | Complete DESeq2 results for all genes and cell types |
| `overall_summary.csv` | 239 B | Quick overview statistics |
| `significant_genes_padj_0.05.csv` | 1.5 KB | Genes with padj < 0.05 |
| `significant_genes_padj_0.05_abs_log2FC_gt_1.csv` | 768 B | Genes with padj < 0.05 AND |log2FC| > 1 |
| `summary_statistics_per_celltype.csv` | 618 B | Per-celltype statistics |

---

## 🌐 Task 19: Multi-Dataset Pseudobulk Analysis

### Key Files

#### 📈 Visualizations
**Heatmaps** (`heatmaps/`):
- `All_Datasets_log2_fold_change_heatmap.pdf` - Combined heatmap
- `GSE174367_Human_log2_fold_change_heatmap.pdf` - GSE174367 only
- `GSE188545_Human_log2_fold_change_heatmap.pdf` - GSE188545 only

**Summary Plots** (`summary_plots/`):
- `Summary_statistics_per_dataset.pdf` - Per-dataset statistics visualization

#### 📊 Data Files (`data/`)

| File | Size | Description |
|------|------|-------------|
| `All_Datasets_ECM_communicationV4_pseudobulk_DEG_combined.csv` | 49 KB | Combined results from both datasets |
| `GSE174367_ECM_communicationV4_pseudobulk_DEG.csv` | 22 KB | GSE174367-specific results |
| `GSE188545_ECM_communicationV4_pseudobulk_DEG.csv` | 28 KB | GSE188545-specific results |
| `All_Datasets_ECM_communicationV4_pseudobulk_significant_padj_0.05.csv` | 2.2 KB | Significant genes (padj < 0.05) across datasets |
| `All_Datasets_ECM_communicationV4_pseudobulk_significant_padj_0.05_abs_log2FC_gt_1.csv` | 1.1 KB | Significant genes (padj < 0.05, |log2FC| > 1) |
| `Genes_significant_in_multiple_datasets.csv` | 205 B | Genes significant in both datasets |
| `summary_statistics_per_dataset.csv` | 325 B | Overview statistics per dataset |
| `summary_statistics_per_celltype.csv` | 1.1 KB | Per-celltype statistics |

---

## 🎯 Key Differences: Wilcoxauc vs Pseudobulk

| Aspect | Wilcoxauc | Pseudobulk DESeq2 |
|--------|-----------|-------------------|
| **Level** | Single-cell | Sample (aggregated) |
| **Method** | Wilcoxon rank sum test | Negative binomial GLM |
| **Sensitivity** | Higher (detects subtle changes) | Lower (more conservative) |
| **Specificity** | Lower (more false positives) | Higher (fewer false positives) |
| **Best for** | Exploratory analysis | Confirmatory analysis |
| **Output** | AUC (0.5-1.0) | log2FC, padj |

### When to Use Which?
- **Wilcoxauc**: Quick screening, hypothesis generation
- **Pseudobulk**: Validation, publication-quality results
- **Combined**: High-confidence genes = significant in BOTH methods

---

## 📈 How to Interpret Results

### Pseudobulk Output Columns
- **gene**: Gene symbol
- **celltype**: Cell type
- **log2FoldChange**: Log2 fold change (positive = up in AD, negative = down in AD)
- **lfcSE**: Standard error of log2FC
- **stat**: Wald statistic
- **pvalue**: Raw p-value
- **padj**: Adjusted p-value (Benjamini-Hochberg FDR)

### Significance Thresholds
- **padj < 0.05**: Statistically significant (FDR < 5%)
- **|log2FC| > 1**: Biologically meaningful (2-fold change)
- **padj < 0.05 AND |log2FC| > 1**: Both statistical and biological significance

### Example Interpretation
```
Gene: COL1A1
Cell type: Microglia
log2FC: 1.5
padj: 0.01
→ COL1A1 is upregulated 2.8-fold in AD microglia (statistically significant)
```

---

## 🔬 Biological Interpretation

### What Do the Results Mean?

**Stable Genes** (not significant in either method):
- Maintain baseline ECM structure
- May be "housekeeping" ECM genes
- Not involved in AD pathogenesis

**Dysregulated Genes** (significant in pseudobulk):
- Actively respond to AD pathology
- May be therapeutic targets
- Warrant further investigation

**Method-Specific Signals**:
- **Wilcoxauc only**: Cell-level heterogeneity, not detectable at sample level
- **Pseudobulk only**: Consistent changes across samples, robust signal

---

## 🚀 Analysis Workflow

```r
# Pseudobulk workflow (Task 18 & 19)
# 1. Aggregate by sample
sobj_agg <- AggregateExpression(
  sobj,
  group.by = "sample",
  assays = "RNA"
)

# 2. Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# 3. Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "AD", "Control"))

# 4. Filter by cell type
res_celltype <- res[gene %in% ECM_V4_genes]
```

---

## 📝 Next Steps

1. **Compare methods**: Identify genes significant in BOTH wilcoxauc and pseudobulk
2. **Cross-dataset validation**: Find consistent signals across GSE174367 and GSE188545
3. **Literature review**: Investigate known AD associations
4. **Pathway analysis**: GO/KEGG enrichment on significant genes
5. **Version comparison**: Compare ECM V1, V2, V3, V4 results

---

## 🔧 Technical Notes

### Sample ID Mapping
- `AggregateExpression()` converts underscores to dashes
- Created mapping table for DESeq2 metadata
- Critical for correct sample-condition assignment

### Cell Type Filtering
- DESeq2 run separately for each cell type
- More conservative than multi-celltype models
- Prevents confounding by celltype-specific effects

### Metadata Variables
- **GSE174367**: `SampleID` (aggregation), `Diagnosis` (condition)
- **GSE188545**: `sample` (aggregation), `condition` (condition)

---

**Analysis Date**: January 8, 2026
**Project**: AD Public Single-Cell Data Analysis
**Analyst**: Claude Code AI Assistant
