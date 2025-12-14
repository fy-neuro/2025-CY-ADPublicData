# AUC Analysis for GSE188545 ECM-related Genes
# Template script - adapt as needed for specific analysis

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
output_dir <- "task/20251214-11-GSE188545_AUC"
dir.create(file.path(output_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plot"), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. LOAD ANNOTATED SEURAT OBJECT
# ---------------------------------------------------------------------------
cat("Loading annotated Seurat object...\n")
sobj_path <- "task/20251214-10-GSE188545_Annotation/GSE188545_sobj_annotated_fast.rds"

if (!file.exists(sobj_path)) {
  stop("Annotated Seurat object not found at: ", sobj_path)
}

sobj <- readRDS(sobj_path)
cat("Object loaded:", ncol(sobj), "cells,", nrow(sobj), "genes\n")
cat("Cell types:", toString(unique(sobj$celltype)), "\n")
cat("Conditions:", toString(unique(sobj$condition)), "\n")

# Optional: Join layers if using Seurat v5 layered data
# sobj <- JoinLayers(sobj)

# ---------------------------------------------------------------------------
# 2. LOAD ECM GENE LISTS
# ---------------------------------------------------------------------------
cat("Loading ECM gene lists...\n")
ecm_data_dir <- "data/ECM_related_genes"
gene_list_files <- list.files(ecm_data_dir, pattern = "\\.xlsx$", full.names = TRUE)

if (length(gene_list_files) == 0) {
  stop("No ECM gene list files found in:", ecm_data_dir)
}

# Load all gene lists
gene_lists <- list()
for (file in gene_list_files) {
  list_name <- gsub("\\.xlsx$", "", basename(file))
  gene_data <- read.xlsx(file, sheet = 1)
  # Assume first column contains gene symbols
  gene_lists[[list_name]] <- gene_data[[1]]
  cat("Loaded", length(gene_lists[[list_name]]), "genes from", list_name, "\n")
}

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
# 4. RUN ANALYSIS FOR EACH GENE LIST
# ---------------------------------------------------------------------------
all_results <- list()

for (list_name in names(gene_lists)) {
  cat("\n=== Analyzing gene list:", list_name, "===\n")
  
  genes <- gene_lists[[list_name]]
  cat("Testing", length(genes), "genes\n")
  
  # Calculate AUC
  auc_results <- calculate_auc_by_celltype(sobj, genes)
  
  if (nrow(auc_results) == 0) {
    cat("No results for this gene list\n")
    next
  }
  
  # Save results
  output_file <- file.path(output_dir, "Documents", 
                           paste0(list_name, "_AUC_per_celltype.csv"))
  write.csv(auc_results, output_file, row.names = FALSE)
  cat("Saved results to:", output_file, "\n")
  
  all_results[[list_name]] <- auc_results
  
  # ---------------------------------------------------------------------------
  # 5. CREATE HEATMAP
  # ---------------------------------------------------------------------------
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
                              paste0(list_name, "_AUC_heatmap.pdf"))
    
    pdf(heatmap_file, width = max(10, ncol(auc_matrix) * 0.3), 
        height = max(6, nrow(auc_matrix) * 0.5))
    
    pheatmap(
      auc_matrix,
      main = paste("AUC of", list_name, "Genes Across Cell Types"),
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
}

# ---------------------------------------------------------------------------
# 6. CREATE SUMMARY REPORT
# ---------------------------------------------------------------------------
cat("\n=== Creating summary report ===\n")

summary_stats <- data.frame()
for (list_name in names(all_results)) {
  if (nrow(all_results[[list_name]]) > 0) {
    stats <- all_results[[list_name]] %>%
      filter(group == "AD") %>%
      summarise(
        Gene_List = list_name,
        Cell_Types = n_distinct(subclass),
        Genes_Tested = n_distinct(gene),
        Mean_AUC = mean(auc, na.rm = TRUE),
        SD_AUC = sd(auc, na.rm = TRUE),
        Genes_AUC_gt_0.7 = sum(auc > 0.7, na.rm = TRUE),
        Genes_AUC_gt_0.8 = sum(auc > 0.8, na.rm = TRUE)
      )
    summary_stats <- bind_rows(summary_stats, stats)
  }
}

if (nrow(summary_stats) > 0) {
  summary_file <- file.path(output_dir, "Documents", "auc_analysis_summary.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("Summary statistics saved to:", summary_file, "\n")
  
  # Print summary to console
  print(summary_stats)
}

# ---------------------------------------------------------------------------
# 7. SESSION INFO
# ---------------------------------------------------------------------------
session_file <- file.path(output_dir, "Documents", "session_info_auc.txt")
sink(session_file)
cat("AUC Analysis Session Information\n")
cat("================================\n\n")
cat("Analysis date:", date(), "\n\n")
print(sessionInfo())
sink()

cat("\n=== Analysis complete ===\n")
cat("Outputs saved to:", output_dir, "\n")