# Load GSE188545 10X Genomics data and create Seurat object
# Direct loading without symlinks

# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(here)
library(ggplot2)
library(patchwork)

# Source project utilities
source(here("src/getdir.R"))

# Set working directory to project root for consistent paths
project_root <- here::here()
setwd(project_root)
cat("Working directory set to:", getwd(), "\n")

# Configuration
data_dir <- "data/GSE188545/GSE188545_GEM/"
output_dir <- "task/20251214-9-GSE188545-Load_Data/"

# Create directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plot"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)

cat("=== GSE188545 Data Loading Pipeline ===\n")
cat("Data directory:", data_dir, "\n")
cat("Output directory:", output_dir, "\n")

# Step 1: Identify samples
cat("\n1. Identifying samples...\n")

# List all barcode files to extract sample prefixes
barcode_files <- list.files(data_dir, pattern = "_barcodes.tsv.gz$", full.names = FALSE)
sample_prefixes <- gsub("_barcodes.tsv.gz$", "", barcode_files)

cat("Found", length(sample_prefixes), "samples:\n")
print(sample_prefixes)

# Step 2: Custom function to read 10X data with file prefixes
read_10x_with_prefix <- function(sample_prefix, data_dir) {
  # Construct file paths
  barcode_file <- file.path(data_dir, paste0(sample_prefix, "_barcodes.tsv.gz"))
  gene_file <- file.path(data_dir, paste0(sample_prefix, "_genes.tsv.gz"))
  matrix_file <- file.path(data_dir, paste0(sample_prefix, "_matrix.mtx.gz"))
  
  cat("  Reading:", sample_prefix, "\n")
  cat("    Barcodes:", barcode_file, "\n")
  cat("    Genes:", gene_file, "\n")
  cat("    Matrix:", matrix_file, "\n")
  
  # Check if files exist
  if (!all(file.exists(c(barcode_file, gene_file, matrix_file)))) {
    cat("    ERROR: One or more files missing\n")
    return(NULL)
  }
  
  # Read the three files
  # Read barcodes
  barcodes <- read.table(gzfile(barcode_file), header = FALSE, stringsAsFactors = FALSE)[,1]
  
  # Read genes/features (two columns: gene_id, gene_name)
  genes <- read.table(gzfile(gene_file), header = FALSE, stringsAsFactors = FALSE)
  # Some 10X files have 3 columns (with gene type), we need column 2 for gene names
  if (ncol(genes) >= 2) {
    gene_names <- genes[,2]
  } else {
    gene_names <- genes[,1]
  }
  
  # Read matrix (MatrixMarket format)
  matrix <- readMM(gzfile(matrix_file))
  
  # Set row and column names
  rownames(matrix) <- gene_names
  colnames(matrix) <- barcodes
  
  # Remove duplicate gene names by keeping first occurrence
  # (Only 24 duplicates out of 33k genes, aggregation would be complex)
  if (any(duplicated(rownames(matrix)))) {
    dup_count <- sum(duplicated(rownames(matrix)))
    cat("    Removing", dup_count, "duplicate gene names (keeping first occurrence)\n")
    matrix <- matrix[!duplicated(rownames(matrix)), ]
  }
  
  cat("    Matrix dimensions:", dim(matrix), "\n")
  return(matrix)
}

# Step 3: Load samples
cat("\n2. Loading samples...\n")

seurat_list <- list()
sample_stats <- data.frame()

for (i in seq_along(sample_prefixes)) {
  sample <- sample_prefixes[i]
  cat("\n[", i, "/", length(sample_prefixes), "]", sample, "\n")
  
  # Read count matrix
  counts <- read_10x_with_prefix(sample, data_dir)
  
  if (is.null(counts)) {
    cat("  Failed to load sample\n")
    next
  }
  
  # Create Seurat object
  sobj <- CreateSeuratObject(
    counts = counts,
    project = sample,
    min.cells = 3,
    min.features = 200
  )
  
  # Add sample metadata
  sobj$sample <- sample
  sobj$condition <- ifelse(grepl("^AD", sample), "AD", "HC")
  sobj$sample_id <- gsub("^(AD|HC)(\\d+).*", "\\2", sample)
  sobj$tissue <- "MTG"
  
  # Calculate QC metrics
  sobj$percent.mt <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj$percent.ribo <- PercentageFeatureSet(sobj, pattern = "^RP[SL]")
  
  cat("  Cells:", ncol(sobj), "Features:", nrow(sobj), "\n")
  cat("  Median nCount_RNA:", round(median(sobj$nCount_RNA)), "\n")
  cat("  Median nFeature_RNA:", round(median(sobj$nFeature_RNA)), "\n")
  cat("  Median % MT:", round(median(sobj$percent.mt), 2), "%\n")
  
  # Store object
  seurat_list[[sample]] <- sobj
  
  # Collect stats
  stats <- data.frame(
    sample = sample,
    condition = sobj$condition[1],
    n_cells = ncol(sobj),
    n_features = nrow(sobj),
    median_nCount = median(sobj$nCount_RNA),
    median_nFeature = median(sobj$nFeature_RNA),
    median_pct_mt = median(sobj$percent.mt)
  )
  sample_stats <- rbind(sample_stats, stats)
  
  # Clear memory
  rm(counts, sobj)
  gc()
}

# Save sample statistics
write.csv(sample_stats, file.path(output_dir, "Documents/sample_statistics.csv"), row.names = FALSE)
cat("\nSample statistics saved.\n")

# Step 4: Merge samples
cat("\n3. Merging samples...\n")

if (length(seurat_list) == 0) {
  stop("No samples loaded successfully!")
}

# Merge all samples
sobj_merged <- merge(seurat_list[[1]], y = seurat_list[-1], 
                     add.cell.ids = names(seurat_list),
                     project = "GSE188545")

cat("  Total cells:", ncol(sobj_merged), "\n")
cat("  Total features:", nrow(sobj_merged), "\n")
cat("  Samples:", length(unique(sobj_merged$sample)), "\n")
cat("  AD cells:", sum(sobj_merged$condition == "AD"), "\n")
cat("  HC cells:", sum(sobj_merged$condition == "HC"), "\n")

# Step 5: Save raw merged object
cat("\n4. Saving raw merged object...\n")
saveRDS(sobj_merged, file = file.path(output_dir, "GSE188545_sobj_raw.rds"))
cat("  Saved to:", file.path(output_dir, "GSE188545_sobj_raw.rds"), "\n")

# Step 6: Generate QC plots
cat("\n5. Generating QC plots...\n")

# Violin plots by sample
qc_violin <- VlnPlot(sobj_merged,
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     group.by = "sample",
                     pt.size = 0,
                     ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "plot/QC_violin_by_sample.pdf"),
       plot = qc_violin, width = 15, height = 6)

# Scatter plots
qc_scatter1 <- FeatureScatter(sobj_merged, 
                             feature1 = "nCount_RNA", 
                             feature2 = "percent.mt",
                             group.by = "condition")

qc_scatter2 <- FeatureScatter(sobj_merged,
                             feature1 = "nCount_RNA",
                             feature2 = "nFeature_RNA",
                             group.by = "condition")

qc_scatter_combined <- qc_scatter1 + qc_scatter2 +
  plot_layout(guides = "collect")

ggsave(file.path(output_dir, "plot/QC_scatter_plots.pdf"),
       plot = qc_scatter_combined, width = 12, height = 6)

# QC summary by condition
qc_summary <- sobj_merged@meta.data %>%
  group_by(condition) %>%
  summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_pct_mt = median(percent.mt),
    .groups = "drop"
  )

write.csv(qc_summary, file.path(output_dir, "Documents/QC_summary.csv"), row.names = FALSE)

# Step 7: Final summary
cat("\n=== GSE188545 Data Loading Complete ===\n")
cat("\nSummary:\n")
cat("  Cells:", ncol(sobj_merged), "\n")
cat("  Features:", nrow(sobj_merged), "\n")
cat("  Samples:", length(unique(sobj_merged$sample)), "\n")
cat("  AD samples:", sum(sobj_merged$condition == "AD"), "cells\n")
cat("  HC samples:", sum(sobj_merged$condition == "HC"), "cells\n")
cat("\nOutput files:\n")
cat("  Raw object:", file.path(output_dir, "GSE188545_sobj_raw.rds"), "\n")
cat("  QC plots:", file.path(output_dir, "plot"), "\n")
cat("  Documents:", file.path(output_dir, "Documents"), "\n")

# Save session info
sink(file.path(output_dir, "Documents/session_info.txt"))
sessionInfo()
sink()

cat("\nNext steps: Run 2.0QC_Filtering_Integration.R for QC filtering and integration.\n")
cat("\nScript completed successfully.\n")