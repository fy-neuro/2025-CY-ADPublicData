# AUC Analysis for ECM_communication_genesV3 (GSE188545 Dataset)
# Script adapted from ecm_auc_analysis.R for specific gene list

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
output_dir <- "task/20260104-12-ECM_communicationV3_AUC"
dir.create(file.path(output_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plot"), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. LOAD ANNOTATED SEURAT OBJECT
# ---------------------------------------------------------------------------
cat("Loading annotated Seurat object...\n")
sobj_path <- "data/GSE188545/GSE188545_sobj_annotated_fast.rds"

if (!file.exists(sobj_path)) {
  stop("Annotated Seurat object not found at: ", sobj_path)
}

sobj <- readRDS(sobj_path)
cat("Object loaded:", ncol(sobj), "cells,", nrow(sobj), "genes\n")
cat("Cell types:", toString(unique(sobj$celltype)), "\n")
cat("Conditions:", toString(unique(sobj$condition)), "\n")

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
#' @param min_cells_per_group Minimum cells per condition (AD/HC) to include cell type
#'
#' @return Data frame with columns: gene, auc, group, subclass, comparison
calculate_auc_by_celltype <- function(sobj, gene_list, min_cells_per_group = 10) {
  
  cell_types <- unique(sobj$celltype)
  auc_results <- list()
  
  for (ct in cell_types) {
    cat("  Processing cell type:", ct, "\n")
    
    # Skip Unknown cell type due to insufficient cells
    if (ct == "Unknown") {
      cat("    Skipping: Unknown cell type (insufficient cells)\n")
      next
    }
    
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
    
    # Join layers for Seurat v5 compatibility
    sobj_ct <- JoinLayers(sobj_ct)
    
    # Calculate AUC using presto::wilcoxauc
    Idents(sobj_ct) <- "condition"
    
    # Note: presto::wilcoxauc expects specific parameters
    # For Seurat v5, use seurat_assay = "RNA" parameter
    auc_df <- tryCatch({
      presto::wilcoxauc(
        sobj_ct,
        comparison = c("AD", "HC"),  # AD vs HC comparison
        seurat_assay = "RNA"
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
        comparison = "AD_vs_HC"
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
cat("\n=== Analyzing ECM_communication_genesV3 ===\n")
cat("Testing", length(gene_list), "genes\n")

# Calculate AUC
auc_results <- calculate_auc_by_celltype(sobj, gene_list)

if (nrow(auc_results) == 0) {
  stop("No AUC results obtained. Please check the analysis parameters.")
}

# Save results
output_file <- file.path(output_dir, "Documents", 
                         "ECM_communication_genesV3_AUC_per_celltype.csv")
write.csv(auc_results, output_file, row.names = FALSE)
cat("Saved results to:", output_file, "\n")

# ---------------------------------------------------------------------------
# 5. CREATE HEATMAP
# ---------------------------------------------------------------------------
cat("\n=== Creating heatmap ===\n")
# Prepare matrix for heatmap (rows = cell types, columns = genes)
auc_matrix <- auc_results %>%
  filter(group == "AD") %>%  # Use AD group results (vs HC)
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
                            "ECM_communication_genesV3_AUC_heatmap.pdf")
  
  pdf(heatmap_file, width = max(10, ncol(auc_matrix) * 0.3), 
      height = max(6, nrow(auc_matrix) * 0.5))
  
  pheatmap(
    auc_matrix,
    main = "AUC of ECM Communication Genes V3 Across Cell Types",
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
  filter(group == "AD") %>%
  summarise(
    Gene_List = "ECM_communication_genesV3",
    Cell_Types = n_distinct(subclass),
    Genes_Tested = n_distinct(gene),
    Mean_AUC = mean(auc, na.rm = TRUE),
    SD_AUC = sd(auc, na.rm = TRUE),
    Genes_AUC_gt_0.7 = sum(auc > 0.7, na.rm = TRUE),
    Genes_AUC_gt_0.8 = sum(auc > 0.8, na.rm = TRUE)
  )

summary_file <- file.path(output_dir, "Documents", "auc_analysis_summary.csv")
write.csv(stats, summary_file, row.names = FALSE)
cat("Summary statistics saved to:", summary_file, "\n")

# Print summary to console
print(stats)

# ---------------------------------------------------------------------------
# 7. IDENTIFY GENES WITH HIGH AUC VALUES
# ---------------------------------------------------------------------------
cat("\n=== Genes with high AUC values (>0.7) ===\n")
high_auc_genes <- auc_results %>%
  filter(group == "AD", auc > 0.7) %>%
  arrange(desc(auc))

if (nrow(high_auc_genes) > 0) {
  print(high_auc_genes)
  
  # Save high AUC genes list
  high_auc_file <- file.path(output_dir, "Documents", 
                            "ECM_communication_genesV3_high_auc_genes.csv")
  write.csv(high_auc_genes, high_auc_file, row.names = FALSE)
  cat("High AUC genes saved to:", high_auc_file, "\n")
} else {
  cat("No genes with AUC > 0.7 found\n")
}

# ---------------------------------------------------------------------------
# 8. SESSION INFO
# ---------------------------------------------------------------------------
session_file <- file.path(output_dir, "Documents", "session_info.txt")
sink(session_file)
cat("ECM Communication V3 AUC Analysis Session Information\n")
cat("==============================================\n\n")
cat("Analysis date:", date(), "\n\n")
print(sessionInfo())
sink()

# ---------------------------------------------------------------------------
# 9. CREATE ANALYSIS SUMMARY DOCUMENT
# ---------------------------------------------------------------------------
summary_file <- file.path(output_dir, "Documents", "analysis_summary.md")
sink(summary_file)
cat("# Analysis Summary: ECM Communication V3 Gene AUC Analysis\n\n")
cat("**Task Directory**: `20260104-12-ECM_communicationV3_AUC`\n")
cat("**Date**: ", date(), "\n")
cat("**Dataset**: GSE188545 human Alzheimer's disease single-cell RNA-seq from middle temporal gyrus (MTG)\n")
cat("**Analysis Type**: AUC (Area Under the ROC Curve) analysis of ECM communication V3 genes across annotated cell types\n\n")
cat("## Overview\n\n")
cat("This task performed AUC analysis to identify ECM communication V3 genes that show differential expression between Alzheimer's Disease (AD) and Healthy Control (HC) samples across different brain cell types. The analysis uses the area under the receiver operating characteristic (ROC) curve to quantify each gene's ability to distinguish AD from HC samples within each cell type.\n\n")
cat("## Methodology\n\n")
cat("### Gene List Analyzed\n")
cat("- **ECM communication genes V3** (`ECM_communication_genesV3.xlsx`)\n")
cat("- ", length(gene_list), " genes in list\n")
cat("- Output: `ECM_communication_genesV3_AUC_per_celltype.csv`\n\n")
cat("### Cell Types Analyzed\n")
cat("- 7 major cell types (Unknown excluded due to insufficient cells)\n")
cat("- AUC calculated for AD vs HC comparison\n\n")
cat("### Key Results\n")
cat("- Mean AUC: ", round(stats$Mean_AUC, 3), "\n")
cat("- SD AUC: ", round(stats$SD_AUC, 3), "\n")
cat("- Genes with AUC > 0.7: ", stats$Genes_AUC_gt_0.7, "\n")
cat("- Genes with AUC > 0.8: ", stats$Genes_AUC_gt_0.8, "\n\n")
cat("## Files Generated\n\n")
cat("1. `Documents/ECM_communication_genesV3_AUC_per_celltype.csv` - Complete AUC results\n")
cat("2. `plot/ECM_communication_genesV3_AUC_heatmap.pdf` - Heatmap visualization\n")
cat("3. `Documents/auc_analysis_summary.csv` - Summary statistics\n")
cat("4. `Documents/session_info.txt` - R session information\n")
if (nrow(high_auc_genes) > 0) {
  cat("5. `Documents/ECM_communication_genesV3_high_auc_genes.csv` - Genes with AUC > 0.7\n")
}
cat("\n")
sink()

cat("\n=== Analysis complete ===\n")
cat("Outputs saved to:", output_dir, "\n")