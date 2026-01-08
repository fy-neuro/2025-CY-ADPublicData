# ECM Communication V4 Wilcoxauc Analysis Summary

**Task Directory**: `20260108-17-ECM_communicationV4_wilcoxauc_3datasets`
**Date**:  Thu Jan  8 11:49:51 2026 

## Overview

This task performed AUC (Area Under the ROC Curve) analysis of ECM communication V4 genes across three Alzheimer's disease single-cell datasets to identify genes differentially expressed between disease and control conditions across different brain cell types.

## Datasets Analyzed

1. **GSE138852** (Mouse AD model)
   - Used oupSample.batchCond as condition
   - Used oupSample.cellType as celltype
2. **GSE174367** (Human snRNA-seq + snATAC-seq multiome)
   - Used Diagnosis as condition
   - Used Cell.Type as celltype
3. **GSE188545** (Human MTG single-cell RNA-seq)
   - Used condition as AD vs HC
   - Used celltype

## Gene List Analyzed

- **ECM communication genes V4** (`ECM_communication_genesV4.xlsx`)
-  22  genes in list

## Methodology

### Analysis Method
- Used `presto::wilcoxauc()` for AUC calculation
- Comparison: Disease vs Control conditions
- Minimum 10 cells per condition required per cell type
- Seurat v5 compatibility: Used `JoinLayers()` before analysis

### Key Results

- Datasets analyzed: 2 
- Total gene-cell type combinations: 600 
- Mean AUC across all datasets: 0.5 
- Genes with AUC > 0.7: 0 

### Per Dataset Statistics

#### GSE174367_Human 
- Cell types analyzed: 7 
- Genes tested: 20 
- Mean AUC: 0.5 
- Genes with AUC > 0.7: 0 

#### GSE188545_Human 
- Cell types analyzed: 8 
- Genes tested: 20 
- Mean AUC: 0.5 
- Genes with AUC > 0.7: 0 

## Files Generated

### Individual Dataset Results
1. `Documents/GSE138852_ECM_communicationV4_wilcoxauc.csv` - GSE138852 results
2. `Documents/GSE174367_ECM_communicationV4_wilcoxauc.csv` - GSE174367 results
3. `Documents/GSE188545_ECM_communicationV4_wilcoxauc.csv` - GSE188545 results

### Combined Results
4. `Documents/All_Datasets_ECM_communicationV4_wilcoxauc_combined.csv` - All datasets combined
5. `Documents/summary_statistics_per_dataset.csv` - Summary statistics per dataset
6. `Documents/summary_statistics_per_celltype.csv` - Summary statistics per cell type
7. `Documents/All_Datasets_ECM_communicationV4_high_auc_genes.csv` - Genes with AUC > 0.7

### Visualizations
8. `plot/GSE138852_ECM_communicationV4_AUC_heatmap.pdf` - GSE138852 heatmap
9. `plot/GSE174367_ECM_communicationV4_AUC_heatmap.pdf` - GSE174367 heatmap
10. `plot/GSE188545_ECM_communicationV4_AUC_heatmap.pdf` - GSE188545 heatmap
11. `plot/All_Datasets_ECM_communicationV4_AUC_heatmap.pdf` - Combined heatmap
12. `plot/Summary_statistics_per_dataset.pdf` - Summary bar plots

## AUC Interpretation

- **AUC = 0.5**: No discriminatory power (equal expression)
- **AUC > 0.5**: Higher expression in disease group
- **AUC < 0.5**: Higher expression in control group
- **AUC ≥ 0.7**: Moderate discriminatory power
- **AUC ≥ 0.8**: Strong discriminatory power

