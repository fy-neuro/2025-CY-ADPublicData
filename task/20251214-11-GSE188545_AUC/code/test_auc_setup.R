# Test script for GSE188545 AUC analysis setup
# Verifies that required files and packages are available

library(Seurat)
library(here)

# Set working directory to project root
project_root <- here::here()
setwd(project_root)
cat("Working directory:", getwd(), "\n")

# Check for annotated Seurat object
annotated_path <- "task/20251214-10-GSE188545_Annotation/GSE188545_sobj_annotated_fast.rds"
cat("Checking for annotated object:", annotated_path, "\n")

if (file.exists(annotated_path)) {
  cat("✓ Annotated object found\n")
  
  # Load object (just metadata to save time)
  cat("Loading object metadata...\n")
  sobj <- readRDS(annotated_path)
  
  cat("Object dimensions:", ncol(sobj), "cells,", nrow(sobj), "genes\n")
  cat("Available metadata columns:", toString(names(sobj@meta.data)), "\n")
  
  if ("celltype" %in% names(sobj@meta.data)) {
    cell_types <- unique(sobj$celltype)
    cat("✓ Cell type annotation found:", length(cell_types), "types\n")
    cat("  Cell types:", toString(cell_types), "\n")
    
    # Count cells per type
    ct_counts <- table(sobj$celltype)
    cat("\nCell type distribution:\n")
    print(ct_counts)
  } else {
    cat("✗ 'celltype' metadata column not found\n")
  }
  
  if ("condition" %in% names(sobj@meta.data)) {
    diag_counts <- table(sobj$condition)
    cat("\n✓ condition metadata found:\n")
    print(diag_counts)
    
    # Check for both AD and HC
    if (all(c("AD", "HC") %in% names(diag_counts))) {
      cat("✓ Both AD and HC groups present\n")
    } else {
      cat("✗ Missing AD or HC groups\n")
    }
  } else {
    cat("✗ 'condition' metadata column not found\n")
  }
  
} else {
  cat("✗ Annotated object not found\n")
  cat("  Expected at:", annotated_path, "\n")
}

# Check for ECM gene lists
ecm_dir <- "data/ECM_related_genes"
cat("\nChecking for ECM gene lists in:", ecm_dir, "\n")

if (dir.exists(ecm_dir)) {
  ecm_files <- list.files(ecm_dir, pattern = "\\.xlsx$")
  cat("Found", length(ecm_files), "Excel files:\n")
  for (f in ecm_files) {
    cat("  -", f, "\n")
  }
} else {
  cat("✗ ECM directory not found\n")
}

# Check required packages
cat("\nChecking required packages...\n")
required_pkgs <- c("Seurat", "presto", "dplyr", "tidyr", "openxlsx", "pheatmap", "viridis")
installed <- sapply(required_pkgs, requireNamespace, quietly = TRUE)

cat("Package status:\n")
for (i in seq_along(installed)) {
  pkg <- names(installed)[i]
  status <- ifelse(installed[i], "✓", "✗")
  cat("  ", status, pkg, "\n")
}

if (!all(installed)) {
  cat("\nMissing packages:", toString(names(installed)[!installed]), "\n")
  cat("Install with: install.packages(c(", 
      paste0('"', names(installed)[!installed], '"', collapse = ", "), 
      "))\n")
}

# Check output directories
output_dir <- "task/20251214-11-GSE188545_AUC"
cat("\nChecking output directory structure:", output_dir, "\n")

subdirs <- c("code", "plot", "Documents")
for (sd in subdirs) {
  path <- file.path(output_dir, sd)
  if (dir.exists(path)) {
    cat("  ✓", sd, "directory exists\n")
  } else {
    cat("  ✗", sd, "directory missing\n")
  }
}

cat("\n=== Test complete ===\n")
cat("If all checks pass, AUC analysis can proceed.\n")
cat("Run the main script: Rscript code/ecm_auc_analysis.R\n")