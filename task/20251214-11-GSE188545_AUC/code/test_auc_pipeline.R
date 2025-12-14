# Test AUC pipeline with small subset of cells
library(Seurat)
library(presto)
library(dplyr)
library(tidyr)
library(openxlsx)
library(here)

# Set working directory to project root
project_root <- here::here()
setwd(project_root)

cat("Testing AUC pipeline with 5000 random cells...\n")

# Load annotated object
sobj_path <- "task/20251214-10-GSE188545_Annotation/GSE188545_sobj_annotated_fast.rds"
if (!file.exists(sobj_path)) stop("Object not found")

cat("Loading object...\n")
sobj <- readRDS(sobj_path)
cat("Original object:", ncol(sobj), "cells,", nrow(sobj), "genes\n")

# Take random sample of 5000 cells (stratified by cell type to preserve diversity)
set.seed(123)
cell_indices <- sample(ncol(sobj), 5000)
sobj_small <- sobj[, cell_indices]
cat("Subset object:", ncol(sobj_small), "cells\n")

# Join layers for Seurat v5 compatibility (on subset)
cat("Joining layers...\n")
sobj_small <- JoinLayers(sobj_small)

# Load first ECM gene list (first 50 genes)
ecm_file <- "data/ECM_related_genes/ECM_Genes_Grouped_Final.xlsx"
gene_data <- read.xlsx(ecm_file, sheet = 1)
genes <- gene_data[[1]][1:50]  # first 50 genes
cat("Testing with", length(genes), "genes\n")

# Define AUC function (with JoinLayers already done)
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

# Run AUC analysis
cat("\nRunning AUC analysis...\n")
auc_results <- calculate_auc_by_celltype(sobj_small, genes, min_cells_per_group = 5)

cat("\nResults summary:\n")
print(head(auc_results))
cat("Total rows:", nrow(auc_results), "\n")
cat("Cell types analyzed:", toString(unique(auc_results$subclass)), "\n")

# Save test results
output_dir <- "task/20251214-11-GSE188545_AUC/Documents"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
test_file <- file.path(output_dir, "auc_pipeline_test.csv")
write.csv(auc_results, test_file, row.names = FALSE)
cat("Test results saved to:", test_file, "\n")

# Clean up
rm(sobj, sobj_small)
gc()

cat("\n=== Pipeline test complete ===\n")