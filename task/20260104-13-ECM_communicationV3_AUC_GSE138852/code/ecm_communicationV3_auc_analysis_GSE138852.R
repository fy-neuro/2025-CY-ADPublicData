# AUC Analysis for ECM_communication_genesV3 (GSE138852 Dataset)
# Updated script to handle data loading and preprocessing

# Load required libraries
library(Seurat)
library(presto)      # For wilcoxauc function (AUC calculation)
library(dplyr)
library(tidyr)
library(tibble)      # For column_to_rownames function
library(openxlsx)    # For reading Excel gene lists
library(pheatmap)    # For heatmap visualization
library(viridis)     # Color palette
library(here)        # Path management

# Set working directory to project root for consistent paths
project_root <- here::here()
setwd(project_root)

# Create output directories
output_dir <- "task/20260104-13-ECM_communicationV3_AUC_GSE138852"
dir.create(file.path(output_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plot"), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. LOAD GSE138852 RECLUSTERED DATA (with author annotations)
# ---------------------------------------------------------------------------
cat("Loading GSE138852 reclustered data with author annotations...\n")

# Load the reclustered object with author annotations
sobj_path <- "data/GSE138852/GSE138852_sobj_reclustered.rds"
if (!file.exists(sobj_path)) {
  stop("Reclustered GSE138852 object not found at: ", sobj_path)
}

sobj <- readRDS(sobj_path)
cat("Loaded reclustered GSE138852 object:", ncol(sobj), "cells,", nrow(sobj), "genes\n")

# Check available metadata columns
cat("Available metadata columns:", paste(colnames(sobj@meta.data), collapse = ", "), "\n")

# Find condition column (look for AD/WT or disease status)
condition_col <- NULL
possible_cols <- c("genotype", "condition", "group", "disease", "diagnosis", "status")

# First check for explicit condition columns
for (col in possible_cols) {
  if (col %in% colnames(sobj@meta.data)) {
    condition_col <- col
    break
  }
}

# If not found, look for columns with AD/WT values
if (is.null(condition_col)) {
  for (col in colnames(sobj@meta.data)) {
    unique_vals <- unique(sobj@meta.data[[col]])
    if (any(c("AD", "WT", "Alzheimer", "APP", "PS1", "PSAPP", "Control", "CTRL") %in% unique_vals)) {
      condition_col <- col
      cat("Found potential condition column:", col, "\n")
      cat("Unique values:", paste(unique_vals, collapse = ", "), "\n")
      break
    }
  }
}

if (is.null(condition_col)) {
  stop("No condition column found in metadata. Cannot perform AUC analysis.")
}

# Rename to standard 'condition' column
sobj$condition <- sobj@meta.data[[condition_col]]
# Convert "ct" to "WT" for consistency
sobj$condition[sobj$condition == "ct"] <- "WT"
cat("Using condition column:", condition_col, "\n")
cat("Conditions:", toString(unique(sobj$condition)), "\n")

# Check cell type annotation
celltype_col <- NULL
possible_celltype_cols <- c("celltype", "cell_type", "cell.type", "annotation", "cluster", "seurat_clusters", 
                          "oupSample.cellType", "oupSample.subclustID")

for (col in possible_celltype_cols) {
  if (col %in% colnames(sobj@meta.data)) {
    celltype_col <- col
    break
  }
}

if (is.null(celltype_col)) {
  stop("No cell type annotation found in metadata. Cannot perform AUC analysis.")
}

# Rename to standard 'celltype' column
sobj$celltype <- sobj@meta.data[[celltype_col]]
cat("Using cell type column:", celltype_col, "\n")
cat("Cell types:", toString(unique(sobj$celltype)), "\n")

# Check cell counts per condition and cell type
cat("\nCell counts per condition and cell type:\n")
cell_counts <- table(sobj$condition, sobj$celltype)
print(cell_counts)

# ---------------------------------------------------------------------------
# 2. LOAD ECM COMMUNICATION V3 GENE LIST
# ---------------------------------------------------------------------------
cat("Loading ECM communication V3 gene list...\n")
ecm_file <- "data/ECM_related_genes/ECM_communication_genesV3.xlsx"

if (!file.exists(ecm_file)) {
  stop("ECM communication V3 gene list not found at: ", ecm_file)
}

# Load gene list
gene_data <- read.xlsx(ecm_file, sheet = 1)
gene_list <- gene_data[[1]]  # Assume first column contains gene symbols
gene_list <- unique(gene_list[!is.na(gene_list)])  # Remove NAs and duplicates

cat("Loaded", length(gene_list), "unique genes from ECM_communication_genesV3\n")

# ---------------------------------------------------------------------------
# 3. AUC ANALYSIS FUNCTION
# ---------------------------------------------------------------------------
#' Calculate AUC values for genes across cell types
#'
#' @param sobj Seurat object with 'celltype' and 'condition' metadata
#' @param gene_list Vector of gene symbols to analyze
#' @param min_cells_per_group Minimum cells per condition to include cell type
#'
#' @return Data frame with columns: gene, auc, group, subclass, comparison
calculate_auc_by_celltype <- function(sobj, gene_list, min_cells_per_group = 10) {
  
  cell_types <- unique(sobj$celltype)
  auc_results <- list()
  
  for (ct in cell_types) {
    cat("  Processing cell type:", ct, "\n")
    
    # Subset to current cell type
    sobj_ct <- subset(sobj, subset = celltype == ct)
    
    # Check if both conditions present with sufficient cells
    cell_counts <- table(sobj_ct$condition)
    if (length(cell_counts) < 2) {
      cat("    Skipping: only one condition present\n")
      next
    }
    if (any(cell_counts < min_cells_per_group)) {
      cat("    Skipping: insufficient cells in one condition (", 
          paste(names(cell_counts), cell_counts, sep = "=", collapse = ", "), ")\n")
      next
    }
    
    # Find genes present in this subset
    available_genes <- rownames(sobj_ct)
    genes_to_test <- intersect(gene_list, available_genes)
    
    if (length(genes_to_test) == 0) {
      cat("    Skipping: no genes from list present in this cell type\n")
      next
    }
    
    # For Seurat v5 with SCT assay, extract data directly 
    # We'll get expression data using GetAssayData with layer="data"
    
    # Calculate AUC using presto::wilcoxauc
    Idents(sobj_ct) <- "condition"
    
    # Adjust comparison based on actual condition values
    condition_values <- unique(sobj_ct$condition)
    comparison <- condition_values  # Will use the two actual condition values
    
    auc_df <- tryCatch({
      # For SCT assay, we need to get data in the right format
      # Use seurat_assay = "SCT" explicitly
      presto::wilcoxauc(
        sobj_ct,
        comparison = comparison,
        seurat_assay = "SCT"
      )
    }, error = function(e) {
      cat("    Error in wilcoxauc:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(auc_df)) next
    
    # Filter for genes of interest and format results
    ct_results <- auc_df %>%
      filter(feature %in% genes_to_test) %>%
      select(feature, auc, group) %>%
      mutate(
        subclass = ct,
        comparison = paste(comparison[1], "vs", comparison[2])
      ) %>%
      rename(gene = feature)
    
    auc_results[[ct]] <- ct_results
  }
  
  # Combine all results
  if (length(auc_results) == 0) {
    return(data.frame())
  }
  
  do.call(rbind, auc_results)
}

# ---------------------------------------------------------------------------
# 4. RUN AUC ANALYSIS
# ---------------------------------------------------------------------------
cat("\n=== Analyzing ECM_communication_genesV3 (GSE138852) ===\n")
cat("Testing", length(gene_list), "genes\n")

# Calculate AUC
auc_results <- calculate_auc_by_celltype(sobj, gene_list)

if (nrow(auc_results) == 0) {
  cat("No AUC results obtained. This could be due to insufficient cells or missing cell types/conditions.\n")
  # Create a placeholder summary for demonstration
  stats <- data.frame(
    Gene_List = "ECM_communication_genesV3",
    Dataset = "GSE138852",
    Cell_Types = 0,
    Genes_Tested = 0,
    Mean_AUC = NA,
    SD_AUC = NA,
    Genes_AUC_gt_0.7 = 0,
    Genes_AUC_gt_0.8 = 0
  )
  
  cat("\nNote: Actual analysis failed. This is likely because:\n")
  cat("1. Cell type annotation is missing or incomplete\n")
  cat("2. Condition information is not properly loaded\n")
  cat("3. Insufficient cells per cell type/condition for statistical analysis\n")
  cat("\nFor full analysis, ensure proper cell type annotation and condition metadata.\n")
} else {
  # Save results
  output_file <- file.path(output_dir, "Documents", 
                           "ECM_communication_genesV3_AUC_per_celltype_GSE138852.csv")
  write.csv(auc_results, output_file, row.names = FALSE)
  cat("Saved results to:", output_file, "\n")
  
  # ---------------------------------------------------------------------------
  # 5. CREATE HEATMAP
  # ---------------------------------------------------------------------------
  cat("\n=== Creating heatmap ===\n")
  # Prepare matrix for heatmap (rows = cell types, columns = genes)
  auc_matrix <- auc_results %>%
    filter(group == auc_results$group[1]) %>%  # Use first group
    select(subclass, gene, auc) %>%
    pivot_wider(
      names_from = gene,
      values_from = auc,
      values_fill = 0  # Fill missing values with 0
    ) %>%
    column_to_rownames("subclass") %>%
    as.matrix()
  
  if (nrow(auc_matrix) > 1 && ncol(auc_matrix) > 1) {
    # Create heatmap
    heatmap_file <- file.path(output_dir, "plot", 
                              "ECM_communication_genesV3_AUC_heatmap_GSE138852.pdf")
    
    pdf(heatmap_file, width = max(10, ncol(auc_matrix) * 0.3), 
        height = max(6, nrow(auc_matrix) * 0.5))
    
    pheatmap(
      auc_matrix,
      main = "AUC of ECM Communication Genes V3 Across Cell Types (GSE138852)",
      color = viridis(100),
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      treeheight_row = 10,
      treeheight_col = 10,
      fontsize_row = 10,
      fontsize_col = 8,
      angle_col = 45,
      cellwidth = 12,
      cellheight = 20,
      border_color = NA
    )
    
    dev.off()
    cat("Created heatmap:", heatmap_file, "\n")
  } else {
    cat("Insufficient data for heatmap (need at least 2 cell types and 2 genes)\n")
  }
  
  # ---------------------------------------------------------------------------
  # 6. CREATE SUMMARY STATISTICS
  # ---------------------------------------------------------------------------
  cat("\n=== Creating summary statistics ===\n")
  stats <- auc_results %>%
    filter(group == auc_results$group[1]) %>%  # Use first group
    summarise(
      Gene_List = "ECM_communication_genesV3",
      Dataset = "GSE138852",
      Cell_Types = n_distinct(subclass),
      Genes_Tested = n_distinct(gene),
      Mean_AUC = mean(auc, na.rm = TRUE),
      SD_AUC = sd(auc, na.rm = TRUE),
      Genes_AUC_gt_0.7 = sum(auc > 0.7, na.rm = TRUE),
      Genes_AUC_gt_0.8 = sum(auc > 0.8, na.rm = TRUE)
    )
  
  # ---------------------------------------------------------------------------
  # 7. IDENTIFY GENES WITH HIGH AUC VALUES
  # ---------------------------------------------------------------------------
  cat("\n=== Genes with high AUC values (>0.7) ===\n")
  high_auc_genes <- auc_results %>%
    filter(group == auc_results$group[1], auc > 0.7) %>%
    arrange(desc(auc))
  
  if (nrow(high_auc_genes) > 0) {
    print(high_auc_genes)
    
    # Save high AUC genes list
    high_auc_file <- file.path(output_dir, "Documents", 
                              "ECM_communication_genesV3_high_auc_genes_GSE138852.csv")
    write.csv(high_auc_genes, high_auc_file, row.names = FALSE)
    cat("High AUC genes saved to:", high_auc_file, "\n")
  } else {
    cat("No genes with AUC > 0.7 found\n")
  }
}

# ---------------------------------------------------------------------------
# 8. SAVE SUMMARY STATISTICS
# ---------------------------------------------------------------------------
summary_file <- file.path(output_dir, "Documents", "auc_analysis_summary_GSE138852.csv")
write.csv(stats, summary_file, row.names = FALSE)
cat("Summary statistics saved to:", summary_file, "\n")

# Print summary to console
print(stats)

# ---------------------------------------------------------------------------
# 9. SESSION INFO
# ---------------------------------------------------------------------------
session_file <- file.path(output_dir, "Documents", "session_info_GSE138852.txt")
sink(session_file)
cat("ECM Communication V3 AUC Analysis Session Information (GSE138852)\n")
cat("==========================================================\n\n")
cat("Analysis date:", date(), "\n\n")
print(sessionInfo())
sink()

# ---------------------------------------------------------------------------
# 10. CREATE ANALYSIS SUMMARY DOCUMENT
# ---------------------------------------------------------------------------
summary_file <- file.path(output_dir, "Documents", "analysis_summary_GSE138852.md")
sink(summary_file)
cat("# Analysis Summary: ECM Communication V3 Gene AUC Analysis (GSE138852)\n\n")
cat("**Task Directory**: `20260104-13-ECM_communicationV3_AUC_GSE138852`\n")
cat("**Date**: ", date(), "\n")
cat("**Dataset**: GSE138852 mouse Alzheimer's disease single-cell RNA-seq\n")
cat("**Analysis Type**: AUC (Area Under the ROC Curve) analysis of ECM communication V3 genes across annotated cell types\n\n")
cat("## Overview\n\n")
cat("This task performed AUC analysis to identify ECM communication V3 genes that show differential expression between conditions across different brain cell types in GSE138852 mouse dataset. The analysis uses area under the receiver operating characteristic (ROC) curve to quantify each gene's ability to distinguish between conditions within each cell type.\n\n")
cat("## Methodology\n\n")
cat("### Gene List Analyzed\n")
cat("- **ECM communication genes V3** (`ECM_communication_genesV3.xlsx`)\n")
cat("- ", length(gene_list), " genes in list\n\n")
cat("### Cell Types Analyzed\n")
cat("- Cell types as annotated in dataset\n")
cat("- AUC calculated for condition comparison\n\n")
cat("### Key Results\n")
cat("- Mean AUC: ", ifelse(is.na(stats$Mean_AUC), "N/A", round(stats$Mean_AUC, 3)), "\n")
cat("- SD AUC: ", ifelse(is.na(stats$SD_AUC), "N/A", round(stats$SD_AUC, 3)), "\n")
cat("- Genes with AUC > 0.7: ", stats$Genes_AUC_gt_0.7, "\n")
cat("- Genes with AUC > 0.8: ", stats$Genes_AUC_gt_0.8, "\n\n")
cat("## Notes\n\n")
cat("For full analysis, ensure proper cell type annotation and condition metadata are available in the dataset. The current GSE138852 dataset may require additional preprocessing steps before AUC analysis can be performed.\n\n")
cat("## Files Generated\n\n")
cat("1. `Documents/auc_analysis_summary_GSE138852.csv` - Summary statistics\n")
cat("2. `Documents/session_info_GSE138852.txt` - R session information\n")
if (file.exists(file.path(output_dir, "Documents", "ECM_communication_genesV3_AUC_per_celltype_GSE138852.csv"))) {
  cat("3. `Documents/ECM_communication_genesV3_AUC_per_celltype_GSE138852.csv` - Complete AUC results\n")
  cat("4. `plot/ECM_communication_genesV3_AUC_heatmap_GSE138852.pdf` - Heatmap visualization\n")
  if (file.exists(file.path(output_dir, "Documents", "ECM_communication_genesV3_high_auc_genes_GSE138852.csv"))) {
    cat("5. `Documents/ECM_communication_genesV3_high_auc_genes_GSE138852.csv` - Genes with AUC > 0.7\n")
  }
}
cat("\n")
sink()

cat("\n=== Analysis complete ===\n")
cat("Outputs saved to:", output_dir, "\n")