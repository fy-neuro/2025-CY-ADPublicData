# Fixed Pseudobulk Analysis for ECM Genes in GSE188545

library(Seurat)
library(dplyr)
library(tidyr)
library(openxlsx)
library(DESeq2)
library(here)

# Set working directory
setwd(here::here())

# Load data
sobj <- readRDS("data/GSE188545/GSE188545_sobj_annotated_fast.rds")
cat("Loaded GSE188545 object:", ncol(sobj), "cells,", nrow(sobj), "genes\n")

# Join layers for Seurat v5
if ("Seurat" %in% class(sobj) && packageVersion("Seurat") >= "5.0.0") {
  if (length(Layers(sobj)) > 1) {
    sobj <- JoinLayers(sobj)
  }
}

# Load gene list
gene_data <- read.xlsx("data/ECM_related_genes/ECM_communication_genesV3.xlsx", sheet = 1)
gene_list <- unique(gene_data[[1]][!is.na(gene_data[[1]])])
cat("Loaded", length(gene_list), "genes\n")

# Create output directory
task_dir <- "task/20260104-15-ECM_comparison_wilcoxauc_pseudobulk"
dir.create(file.path(task_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)

# Pseudobulk analysis function
run_pseudobulk_fixed <- function(seurat_obj, genes_of_interest, min_cells = 10) {
  
  cell_types <- unique(seurat_obj$celltype)
  results <- data.frame()
  
  for (ct in cell_types) {
    cat("Processing", ct, "...\n")
    
    # Subset cell type
    ct_obj <- subset(seurat_obj, subset = celltype == ct)
    
    # Check conditions and samples
    sample_counts <- table(ct_obj$sample, ct_obj$condition)
    cat("  Sample counts:\n")
    print(sample_counts)
    
    # Need at least 2 samples per condition for DESeq2
    ad_samples <- sum(sample_counts[, "AD"] > 0)
    hc_samples <- sum(sample_counts[, "HC"] > 0)
    
    if (ad_samples < 2 || hc_samples < 2) {
      cat("  Skipping - insufficient samples (need ≥2 per condition)\n")
      cat("  AD samples:", ad_samples, "HC samples:", hc_samples, "\n")
      next
    }
    
    tryCatch({
      # Create pseudobulk counts
      pseudobulk_counts <- AggregateExpression(
        ct_obj,
        assays = "RNA",
        slot = "counts",
        group.by = "sample",
        return.seurat = FALSE
      )
      
      count_matrix <- pseudobulk_counts$RNA
      cat("  Count matrix dimensions:", dim(count_matrix), "\n")
      
      # Filter to genes of interest
      available_genes <- intersect(genes_of_interest, rownames(count_matrix))
      if (length(available_genes) < 3) {
        cat("  Skipping - insufficient genes (", length(available_genes), ")\n")
        next
      }
      
      count_matrix <- count_matrix[available_genes, ]
      cat("  Genes after filtering:", nrow(count_matrix), "\n")
      
      # Create sample metadata - CRITICAL STEP
      # Get samples that actually have cells in this subset
      samples_in_subset <- colnames(count_matrix)
      condition_for_samples <- ct_obj@meta.data$condition[match(samples_in_subset, ct_obj$sample)]
      
      # Remove duplicates and ensure same order
      unique_order <- !duplicated(samples_in_subset)
      samples_in_subset <- samples_in_subset[unique_order]
      condition_for_samples <- condition_for_samples[unique_order]
      
      # Create clean metadata
      sample_meta <- data.frame(
        row.names = samples_in_subset,
        condition = condition_for_samples,
        stringsAsFactors = FALSE
      )
      
      # Ensure count matrix columns match metadata
      count_matrix <- count_matrix[, samples_in_subset, drop = FALSE]
      
      cat("  Final matrix dimensions:", dim(count_matrix), "\n")
      cat("  Sample metadata:\n")
      print(sample_meta)
      
      # Double-check for NA values
      if (any(is.na(sample_meta$condition))) {
        cat("  ERROR: NA values in condition metadata!\n")
        next
      }
      
      # Create DESeq2 dataset
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = sample_meta,
        design = ~ condition
      )
      
      # Pre-filter and run
      keep <- rowSums(counts(dds) >= 10) >= 2
      dds <- dds[keep, ]
      cat("  Genes after filtering:", nrow(dds), "\n")
      
      if (nrow(dds) < 3) {
        cat("  Skipping - too few genes after filtering\n")
        next
      }
      
      dds <- DESeq(dds, quiet = TRUE)
      
      # Get results
      res <- results(dds, name = "condition_AD_vs_HC")
      
      # Format results
      ct_result <- data.frame(
        gene = rownames(res),
        log2_fold_change = res$log2FoldChange,
        p_value = res$pvalue,
        padj = res$padj,
        baseMean = res$baseMean,
        auc_approx = 1 / (1 + exp(-res$log2FoldChange)),
        celltype = ct,
        method = "pseudobulk_DESeq2",
        stringsAsFactors = FALSE
      )
      
      # Remove NA rows
      ct_result <- ct_result[complete.cases(ct_result), ]
      
      if (nrow(ct_result) > 0) {
        results <- rbind(results, ct_result)
        cat("  Analyzed", nrow(ct_result), "genes successfully\n")
        cat("  Significant genes (padj < 0.05):", sum(ct_result$padj < 0.05), "\n")
        cat("  Genes with |log2FC| > 1:", sum(abs(ct_result$log2_fold_change) > 1), "\n")
      } else {
        cat("  No valid results after filtering\n")
      }
      
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
      # Print more debugging info
      cat("  Debug: sample_meta structure\n")
      if (exists("sample_meta")) {
        print(str(sample_meta))
      }
    })
  }
  
  return(results)
}

# Run fixed pseudobulk analysis
cat("\n=== Running Fixed Pseudobulk Analysis ===\n")
pseudobulk_results <- run_pseudobulk_fixed(sobj, gene_list)

if (nrow(pseudobulk_results) > 0) {
  # Save results
  output_file <- file.path(task_dir, "Documents", "pseudobulk_results_fixed.csv")
  write.csv(pseudobulk_results, output_file, row.names = FALSE)
  cat("\nSaved pseudobulk results:", output_file, "\n")
  
  # Basic summary
  cat("\n=== Pseudobulk Summary ===\n")
  summary_stats <- pseudobulk_results %>%
    group_by(celltype) %>%
    summarise(
      genes_tested = n(),
      mean_auc = mean(auc_approx, na.rm = TRUE),
      significant = sum(padj < 0.05, na.rm = TRUE),
      lfc_gt_1 = sum(abs(log2_fold_change) > 1, na.rm = TRUE),
      lfc_gt_2 = sum(abs(log2_fold_change) > 2, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_stats)
  
  # Save summary
  summary_file <- file.path(task_dir, "Documents", "pseudobulk_summary_fixed.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("Saved summary:", summary_file, "\n")
  
  # Top significant genes
  if (any(pseudobulk_results$padj < 0.05)) {
    top_genes <- pseudobulk_results %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      head(10)
    
    cat("\nTop 10 significant genes:\n")
    print(top_genes)
    
    # Save top genes
    top_file <- file.path(task_dir, "Documents", "pseudobulk_top_significant.csv")
    write.csv(top_genes, top_file, row.names = FALSE)
    cat("Saved top significant genes:", top_file, "\n")
  } else {
    cat("\nNo significant genes found (padj < 0.05)\n")
  }
  
} else {
  cat("\nNo pseudobulk results generated\n")
}

cat("\n=== Pseudobulk Analysis Complete ===\n")