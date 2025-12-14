# Test AUC function with small subset
library(Seurat)
library(presto)
library(dplyr)
library(tidyr)
library(here)

# Set working directory to project root
project_root <- here::here()
setwd(project_root)

cat("Testing AUC calculation function...\n")

# Load annotated object (full)
sobj_path <- "task/20251214-10-GSE188545_Annotation/GSE188545_sobj_annotated_fast.rds"
if (!file.exists(sobj_path)) {
  stop("Annotated object not found")
}

cat("Loading object...\n")
sobj <- readRDS(sobj_path)
cat("Object loaded:", ncol(sobj), "cells,", nrow(sobj), "genes\n")

# Take a small subset for testing: 2 cell types, 1000 cells each
set.seed(123)
cell_types_to_test <- c("Oligodendrocytes", "Astrocytes")
cells_to_keep <- c()

for (ct in cell_types_to_test) {
  ct_cells <- which(sobj$celltype == ct)
  if (length(ct_cells) > 1000) {
    ct_cells <- sample(ct_cells, 1000)
  }
  cells_to_keep <- c(cells_to_keep, ct_cells)
}

sobj_test <- sobj[, cells_to_keep]
cat("Test subset:", ncol(sobj_test), "cells\n")

# Load a small gene list (first 10 genes from ECM list)
library(openxlsx)
ecm_file <- "data/ECM_related_genes/ECM_Genes_Grouped_Final.xlsx"
ecm_genes <- read.xlsx(ecm_file, sheet = 1)[[1]]
test_genes <- ecm_genes[1:10]
cat("Testing with", length(test_genes), "genes:", toString(test_genes), "\n")

# Define the AUC function (copied from template)
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

# Run test
cat("\nRunning AUC calculation...\n")
test_results <- calculate_auc_by_celltype(sobj_test, test_genes, min_cells_per_group = 5)

cat("\nTest results:\n")
print(test_results)
cat("Number of results:", nrow(test_results), "\n")

if (nrow(test_results) > 0) {
  cat("\nTest passed! Function works correctly.\n")
} else {
  cat("\nTest failed: no results returned.\n")
}

# Clean up
rm(sobj, sobj_test)
gc()

cat("\n=== Test complete ===\n")