# FINAL Fixed Pseudobulk Analysis - Correct Sample ID Issue

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

# FINAL Pseudobulk analysis function with SAMPLE ID FIX
run_pseudobulk_final <- function(seurat_obj, genes_of_interest, min_cells = 10) {
  
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
      cat("  Column names in count matrix:", head(colnames(count_matrix)), "\n")
      
      # Filter to genes of interest
      available_genes <- intersect(genes_of_interest, rownames(count_matrix))
      if (length(available_genes) < 3) {
        cat("  Skipping - insufficient genes (", length(available_genes), ")\n")
        next
      }
      
      count_matrix <- count_matrix[available_genes, ]
      cat("  Genes after filtering:", nrow(count_matrix), "\n")
      
      # FIX: Create sample metadata correctly
      # The issue is that AggregateExpression replaces underscores with dashes
      # We need to map the original sample IDs to the new column names
      
      # Get original sample-condition mapping for this cell type
      orig_metadata <- ct_obj@meta.data %>% 
        select(sample, condition) %>% 
        distinct()
      
      cat("  Original samples in metadata:\n")
      print(head(orig_metadata, 10))
      
      # Create mapping from original to modified sample IDs
      # Replace underscores with dashes to match column names
      sample_mapping <- data.frame(
        original_id = orig_metadata$sample,
        modified_id = gsub("_", "-", orig_metadata$sample),
        condition = orig_metadata$condition,
        stringsAsFactors = FALSE
      )
      
      # Remove duplicates
      sample_mapping <- sample_mapping[!duplicated(sample_mapping$modified_id), ]
      
      cat("  Sample mapping (original -> modified):\n")
      print(head(sample_mapping, 10))
      
      # Filter mapping to only samples present in count matrix
      samples_in_matrix <- colnames(count_matrix)
      sample_mapping <- sample_mapping[sample_mapping$modified_id %in% samples_in_matrix, ]
      
      cat("  Filtered sample mapping:\n")
      print(sample_mapping)
      
      # Create sample metadata for DESeq2
      sample_meta <- data.frame(
        row.names = sample_mapping$modified_id,
        condition = sample_mapping$condition,
        stringsAsFactors = FALSE
      )
      
      # Ensure count matrix columns match metadata
      count_matrix <- count_matrix[, sample_mapping$modified_id, drop = FALSE]
      
      cat("  Final matrix dimensions:", dim(count_matrix), "\n")
      cat("  Sample metadata for DESeq2:\n")
      print(sample_meta)
      
      # Check for NA values
      if (any(is.na(sample_meta$condition))) {
        cat("  ERROR: NA values in condition metadata!\n")
        stop("NA values in metadata")
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
        cat("  ✓ Analyzed", nrow(ct_result), "genes successfully\n")
        cat("  ✓ Significant genes (padj < 0.05):", sum(ct_result$padj < 0.05), "\n")
        cat("  ✓ Genes with |log2FC| > 1:", sum(abs(ct_result$log2_fold_change) > 1), "\n")
        
        # Show top significant genes
        if (any(ct_result$padj < 0.05)) {
          top_sig <- ct_result %>% filter(padj < 0.05) %>% arrange(padj) %>% head(3)
          cat("  Top significant genes:\n")
          print(top_sig[, c("gene", "log2_fold_change", "padj")])
        }
      } else {
        cat("  No valid results after filtering\n")
      }
      
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
      return(NULL)
    })
  }
  
  return(results)
}

# Run FINAL fixed pseudobulk analysis
cat("\n=== RUNNING FINAL FIXED PSEUDOBULK ANALYSIS ===\n")
pseudobulk_results <- run_pseudobulk_final(sobj, gene_list)

if (nrow(pseudobulk_results) > 0) {
  # Save results
  output_file <- file.path(task_dir, "Documents", "pseudobulk_results_final.csv")
  write.csv(pseudobulk_results, output_file, row.names = FALSE)
  cat("\n✅ Saved pseudobulk results:", output_file, "\n")
  
  # Basic summary
  cat("\n=== PSEUDOBULK FINAL SUMMARY ===\n")
  summary_stats <- pseudobulk_results %>%
    group_by(celltype) %>%
    summarise(
      genes_tested = n(),
      mean_auc = mean(auc_approx, na.rm = TRUE),
      sd_auc = sd(auc_approx, na.rm = TRUE),
      significant = sum(padj < 0.05, na.rm = TRUE),
      lfc_gt_1 = sum(abs(log2_fold_change) > 1, na.rm = TRUE),
      lfc_gt_2 = sum(abs(log2_fold_change) > 2, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_stats)
  
  # Save summary
  summary_file <- file.path(task_dir, "Documents", "pseudobulk_summary_final.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("✅ Saved summary:", summary_file, "\n")
  
  # Top significant genes overall
  if (any(pseudobulk_results$padj < 0.05)) {
    top_genes <- pseudobulk_results %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      head(10)
    
    cat("\nTop 10 significant genes overall:\n")
    print(top_genes[, c("gene", "celltype", "log2_fold_change", "padj", "auc_approx")])
    
    # Save top genes
    top_file <- file.path(task_dir, "Documents", "pseudobulk_top_significant_final.csv")
    write.csv(top_genes, top_file, row.names = FALSE)
    cat("✅ Saved top significant genes:", top_file, "\n")
  } else {
    cat("\n❌ No significant genes found (padj < 0.05)\n")
  }
  
  # Compare with wilcoxauc results if available
  wilcoxauc_file <- file.path(task_dir, "Documents", "wilcoxauc_results.csv")
  if (file.exists(wilcoxauc_file)) {
    cat("\n=== COMPARING PSEUDOBULK VS WILCOXAUC ===\n")
    wilcoxauc_results <- read.csv(wilcoxauc_file)
    
    # Prepare for comparison
    wilcoxauc_comp <- wilcoxauc_results %>%
      select(gene = feature, auc, celltype) %>%
      rename(wilcoxauc_auc = auc)
    
    pseudobulk_comp <- pseudobulk_results %>%
      select(gene, auc_approx, celltype) %>%
      rename(pseudobulk_auc = auc_approx)
    
    # Merge for comparison
    comparison <- left_join(wilcoxauc_comp, pseudobulk_comp, 
                         by = c("gene", "celltype")) %>%
      filter(!is.na(pseudobulk_auc))
    
    if (nrow(comparison) > 10) {
      # Calculate correlation
      cor_test <- cor.test(comparison$wilcoxauc_auc, comparison$pseudobulk_auc)
      
      cat("Method comparison:\n")
      cat("  Correlation (r):", round(cor_test$estimate, 4), "\n")
      cat("  P-value:", signif(cor_test$p.value, 6), "\n")
      cat("  N comparisons:", nrow(comparison), "\n")
      
      # Save comparison
      comparison_file <- file.path(task_dir, "Documents", "method_comparison_final.csv")
      write.csv(comparison, comparison_file, row.names = FALSE)
      cat("✅ Saved method comparison:", comparison_file, "\n")
    }
  }
  
} else {
  cat("\n❌ No pseudobulk results generated\n")
}

cat("\n=== 🎉 PSEUDOBULK ANALYSIS COMPLETE ===\n")
cat("Output files saved to:", task_dir, "\n")