# Quick test of annotation pipeline
# Loads integrated object and performs basic checks

cat("=== Testing GSE188545 Annotation Pipeline ===\n")

# Load required core packages
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)

# Set working directory
project_root <- here::here()
setwd(project_root)
cat("Working directory:", getwd(), "\n")

# Check if integrated object exists
input_dir <- "task/20251214-9-GSE188545-Load_Data/"
integrated_object <- file.path(input_dir, "GSE188545_sobj_integrated.rds")

if (!file.exists(integrated_object)) {
  stop("Integrated object not found at: ", integrated_object)
}
cat("Integrated object found:", integrated_object, "\n")

# Load object (this may take time and memory)
cat("Loading integrated object...\n")
start_time <- Sys.time()
sobj <- readRDS(integrated_object)
load_time <- Sys.time() - start_time
cat("  Load time:", round(load_time, 2), "seconds\n")

# Basic object checks
cat("\nObject information:\n")
cat("  Project name:", sobj@project.name, "\n")
cat("  Number of cells:", ncol(sobj), "\n")
cat("  Number of features:", nrow(sobj), "\n")
cat("  Assays available:", names(sobj@assays), "\n")
cat("  Reductions available:", names(sobj@reductions), "\n")

# Check metadata
cat("\nMetadata columns:\n")
print(colnames(sobj@meta.data))

# Check clustering
if ("RNA_snn_res.0.5" %in% colnames(sobj@meta.data)) {
  clusters <- unique(sobj$RNA_snn_res.0.5)
  cat("  Clusters (res 0.5):", length(clusters), "\n")
}

# Check condition
if ("condition" %in% colnames(sobj@meta.data)) {
  cat("  Condition distribution:\n")
  print(table(sobj$condition))
} else if ("sample" %in% colnames(sobj@meta.data)) {
  # Infer condition from sample names
  sobj$condition <- ifelse(grepl("_AD", sobj$sample), "AD",
                           ifelse(grepl("_HC", sobj$sample), "HC", "Unknown"))
  cat("  Condition inferred from sample names:\n")
  print(table(sobj$condition))
}

# Check marker genes from src/1stAnnotation.R
cat("\nChecking marker genes...\n")
source(here("src/1stAnnotation.R"))

# Test with a few key markers
key_markers <- c("P2RY12", "TMEM119", "AQP4", "GFAP", "MBP", "MOG", "PDGFRA", "VWF")
available <- key_markers[key_markers %in% rownames(sobj)]
missing <- setdiff(key_markers, rownames(sobj))

cat("  Available markers:", length(available), "/", length(key_markers), "\n")
if (length(missing) > 0) {
  cat("  Missing markers:", paste(missing, collapse = ", "), "\n")
}

# Quick UMAP visualization (optional, comment out if slow)
cat("\nGenerating quick UMAP visualization...\n")
if ("umap" %in% names(sobj@reductions)) {
  # Simple UMAP by sample
  umap_plot <- DimPlot(sobj, 
                      reduction = "umap",
                      group.by = "sample",
                      pt.size = 0.1,
                      raster = TRUE) +
    ggtitle("Samples")
  
  output_dir <- "task/20251214-10-GSE188545_Annotation/"
  dir.create(file.path(output_dir, "plot/test"), showWarnings = FALSE, recursive = TRUE)
  
  ggsave(file.path(output_dir, "plot/test/umap_samples_test.pdf"),
         plot = umap_plot, width = 10, height = 8)
  cat("  Test plot saved to plot/test/umap_samples_test.pdf\n")
}

cat("\n=== Test completed successfully ===\n")
cat("Object size in memory:", format(object.size(sobj), units = "auto"), "\n")
cat("Total memory used:", format(memory.size(), units = "auto"), "\n")