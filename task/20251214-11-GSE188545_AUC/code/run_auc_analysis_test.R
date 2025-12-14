# Test run of AUC analysis with limited genes
library(Seurat)
library(presto)
library(dplyr)
library(tidyr)
library(openxlsx)
library(pheatmap)
library(viridis)
library(here)

# Set working directory to project root
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

# ---------------------------------------------------------------------------
# 2. LOAD ECM GENE LISTS (only first file)
# ---------------------------------------------------------------------------
cat("Loading ECM gene lists...\n")
ecm_data_dir <- "data/ECM_related_genes"
gene_list_files <- list.files(ecm_data_dir, pattern = "\\.xlsx$", full.names = TRUE)

if (length(gene_list_files) == 0) {
  stop("No ECM gene list files found in:", ecm_data_dir)
}

# Load only first gene list for testing
test_file <- gene_list_files[1]
list_name <- gsub("\\.xlsx$", "", basename(test_file))
gene_data <- read.xlsx(test_file, sheet = 1)
genes <- gene_data[[1]]
cat("Testing with gene list:", list_name, "(", length(genes), "genes)\n")

# Take only first 20 genes for testing
test_genes <- genes[1:20]
cat("Subsetting to", length(test_genes), "genes for testing\n")

# ---------------------------------------------------------------------------
# 3. AUC ANALYSIS FUNCTION (updated with JoinLayers)
# ---------------------------------------------------------------------------
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
# 4. RUN ANALYSIS FOR TEST GENE LIST
# ---------------------------------------------------------------------------
cat("\n=== Analyzing gene list:", list_name, "(test subset) ===\n")
cat("Testing", length(test_genes), "genes\n")

# Calculate AUC
auc_results <- calculate_auc_by_celltype(sobj, test_genes, min_cells_per_group = 10)

if (nrow(auc_results) == 0) {
  cat("No results for this gene list\n")
} else {
  # Save results
  output_file <- file.path(output_dir, "Documents", 
                          paste0(list_name, "_AUC_per_celltype_TEST.csv"))
  write.csv(auc_results, output_file, row.names = FALSE)
  cat("Saved test results to:", output_file, "\n")
  
  # ---------------------------------------------------------------------------
  # 5. CREATE HEATMAP (if enough data)
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
                             paste0(list_name, "_AUC_heatmap_TEST.pdf"))
    
    pdf(heatmap_file, width = max(6, ncol(auc_matrix) * 0.3), 
        height = max(4, nrow(auc_matrix) * 0.5))
    
    pheatmap(
      auc_matrix,
      main = paste("AUC of", list_name, "Genes Across Cell Types (Test)"),
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
    cat("Created test heatmap:", heatmap_file, "\n")
  } else {
    cat("Insufficient data for heatmap (need at least 2 cell types and 2 genes)\n")
  }
}

# ---------------------------------------------------------------------------
# 6. CLEAN UP
# ---------------------------------------------------------------------------
rm(sobj)
gc()

cat("\n=== Test analysis complete ===\n")
cat("Outputs saved to:", output_dir, "\n")