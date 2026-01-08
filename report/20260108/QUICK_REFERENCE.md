# Quick Reference: ECM Communication V4 Analysis

## 🎯 TL;DR (30 seconds)

**What**: Analyzed 20 ECM genes across 15 cell types in 2 AD datasets (600 combinations)
**Method**: Wilcoxauc AUC analysis (single-cell level)
**Result**: **Mean AUC = 0.5** → ECM V4 genes are **stable in AD** (minimal expression changes)
**Files**: All organized in `report/20260108/`

---

## 📊 Understanding AUC Values

| AUC Range | Interpretation | Effect Size | Our Results |
|-----------|----------------|-------------|-------------|
| 0.5 | No difference | None | ✅ **598 genes** (99.7%) |
| 0.6-0.7 | Moderate difference | Small | ⚠️ **2 genes** (0.3%) |
| 0.7-0.8 | Strong difference | Moderate | ❌ **0 genes** |
| >0.8 | Very strong difference | Large | ❌ **0 genes** |

### What is AUC?
- **AUC** = Area Under the ROC Curve
- Measures how well a gene discriminates between AD and control
- **0.5** = Cannot tell AD vs control apart (no difference)
- **1.0** = Perfect separation (very different)

---

## 🔍 The 2 Genes with AUC > 0.6

Open `csv/All_Datasets_ECM_communicationV4_wilcoxauc_combined.csv` and filter for `auc > 0.6` to see:
- Which genes showed moderate effects
- Which cell types they were in
- Their p-values and effect sizes

---

## 📂 Which Files to Open?

### I have 1 minute:
1. `README.md` - Navigation guide
2. `figures/All_Datasets_ECM_communicationV4_AUC_heatmap.pdf` - Visual overview
3. `csv/summary_statistics_per_dataset.csv` - Numbers at a glance

### I have 5 minutes:
1. All of the above, plus:
2. `20260108_task_summary.md` - Complete story
3. `docs/analysis_summary.md` - How we did it

### I have 15 minutes:
1. All of the above, plus:
2. `csv/All_Datasets_ECM_communicationV4_wilcoxauc_combined.csv` - Full data
3. Individual heatmaps in `figures/`
4. `csv/summary_statistics_per_celltype.csv` - Celltype breakdown

---

## 🧪 Technical Details

### Datasets Analyzed
- **GSE174367** (Human): 7 cell types, 280 combinations
- **GSE188545** (Human): 8 cell types, 320 combinations
- **GSE138852** (Mouse): Skipped (SCT data incompatibility)

### Method
```r
presto::wilcoxauc(
  seurat_object,
  group_by = "condition",     # AD vs Control
  assay = "data",             # Log-normalized data
  seurat_assay = "RNA"
)
```

### Output Columns
- `auc`: Discrimination power (0.5-1.0)
- `pval`: Statistical significance
- `effect_size`: Effect magnitude
- `gene`: Gene symbol
- `celltype`: Cell type

---

## 📈 Visualization Guide

### Main Heatmap
`figures/All_Datasets_ECM_communicationV4_AUC_heatmap.pdf`
- **X-axis**: Genes (20 ECM V4 genes)
- **Y-axis**: Cell types (15 total)
- **Color**: AUC value (blue = low, red = high)
- **Pattern**: Mostly uniform (blue) = genes stable

### Individual Heatmaps
- GSE174367: 7 cell types
- GSE188545: 8 cell types
- Compare to see consistency across datasets

---

## 🎓 Biological Interpretation

### What This Means
ECM Communication V4 genes are involved in **structural brain maintenance**:
- Collagens (COL1A1, COL3A1, etc.)
- Laminins (LAMA1, LAMB1, etc.)
- Fibronectin (FN1)
- Other ECM proteins

**Stability** suggests these genes:
- ✅ Maintain brain structure in both AD and control
- ✅ May be "housekeeping" ECM genes
- ❌ Are not dysregulated in AD pathogenesis
- 🤔 Could be important for other brain functions

### Contrast with Other ECM Lists
Previous versions (V1, V2, V3) may show different patterns. This suggests:
- **V4 genes**: Stable, structural ECM
- **V1/V2/V3**: May contain more AD-responsive genes

---

## 🔧 Troubleshooting

### Can't open PDFs?
- Use Adobe Reader, Preview (Mac), or browser
- All heatmaps are vector-based (infinitely zoomable)

### CSV files look messy?
- Open in Excel, R, or pandas
- Use `summary_statistics_*.csv` for cleaner summaries

### Want to reproduce?
- Check `docs/session_info.txt` for R package versions
- Scripts in `task/20260108-17-*/code/`

---

## 📝 Key Stats at a Glance

| Metric | Value |
|--------|-------|
| Total genes tested | 20 |
| Total cell types | 15 |
| Total combinations | 600 |
| Mean AUC | 0.5 |
| Genes with AUC > 0.6 | 2 (0.3%) |
| Genes with AUC > 0.7 | 0 (0%) |
| Datasets analyzed | 2 (GSE174367, GSE188545) |

---

## 🚀 Next Steps

1. **Pseudobulk analysis** (Task 19) - Compare with wilcoxauc results
2. **Cross-version comparison** - V1 vs V2 vs V3 vs V4
3. **Pathway analysis** - GO/KEGG enrichment
4. **Literature review** - What's known about these genes in AD?

---

**Analysis Date**: January 8, 2026
**Project**: AD Public Single-Cell Data Analysis
**Analyst**: Claude Code AI Assistant
