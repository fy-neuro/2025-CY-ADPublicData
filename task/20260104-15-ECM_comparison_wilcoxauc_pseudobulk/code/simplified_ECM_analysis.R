# Simplified ECM Genes Analysis: Presto wilcoxauc vs Pseudobulk
# For GSE188545 Dataset

# Load required libraries
library(Seurat)
library(presto)
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)
library(pheatmap)
library(viridis)
library(ggplot2)
library(here)
library(DESeq2)

# Set working directory
project_root <- here::here()
setwd(project_root)

# Create output directories
task_dir <- "task/20260104-15-ECM_comparison_wilcoxauc_pseudobulk"
dir.create(file.path(task_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(task_dir, "plot"), showWarnings = FALSE, recursive = TRUE)

cat("=== Simplified ECM Genes Analysis ===\n")

# ---------------------------------------------------------------------------
# 1. LOAD GSE188545 DATA
# ---------------------------------------------------------------------------
cat("\n1. Loading GSE188545 data...\n")

sobj <- readRDS("data/GSE188545/GSE188545_sobj_annotated_fast.rds")
cat("Loaded GSE188545 object:", ncol(sobj), "cells,", nrow(sobj), "genes\n")

# Join layers for Seurat v5
if ("Seurat" %in% class(sobj) && packageVersion("Seurat") >= "5.0.0") {
  if (length(Layers(sobj)) > 1) {
    sobj <- JoinLayers(sobj)
  }
}

# ---------------------------------------------------------------------------
# 2. LOAD ONE ECM GENE LIST FOR TESTING
# ---------------------------------------------------------------------------
cat("\n2. Loading ECM communication V3 gene list...\n")

gene_data <- read.xlsx("data/ECM_related_genes/ECM_communication_genesV3.xlsx", sheet = 1)
gene_list <- unique(gene_data[[1]][!is.na(gene_data[[1]])])
cat("Loaded", length(gene_list), "genes\n")

# ---------------------------------------------------------------------------
# 3. SIMPLE WILCOXAUC ANALYSIS
# ---------------------------------------------------------------------------
cat("\n3. Running wilcoxauc analysis...\n")

# Function for wilcoxauc
run_wilcoxauc <- function(seurat_obj, genes_of_interest, min_cells = 10) {
  
  cell_types <- unique(seurat_obj$celltype)
  results <- data.frame()
  
  for (ct in cell_types) {
    cat("  Processing", ct, "\n")
    
    # Subset cell type
    ct_obj <- subset(seurat_obj, subset = celltype == ct)
    
    # Check conditions
    cell_counts <- table(ct_obj$condition)
    if (length(cell_counts) < 2 || any(cell_counts < min_cells)) {
      cat("    Skipping - insufficient cells\n")
      next
    }
    
    # Run wilcoxauc
    Idents(ct_obj) <- "condition"
    
    tryCatch({
      auc_result <- presto::wilcoxauc(ct_obj, seurat_assay = "RNA")
      
      # Check result structure
      if (is.data.frame(auc_result) && "feature" %in% colnames(auc_result)) {
        
        # Filter for genes of interest
        ct_result <- auc_result[auc_result$feature %in% genes_of_interest, ]
        
        if (nrow(ct_result) > 0) {
          ct_result$celltype <- ct
          ct_result$method <- "wilcoxauc"
          results <- rbind(results, ct_result)
          cat("    Found", nrow(ct_result), "genes\n")
        } else {
          cat("    No matching genes found\n")
        }
      }
    }, error = function(e) {
      cat("    Error:", e$message, "\n")
    })
  }
  
  return(results)
}

# Run analysis
wilcoxauc_results <- run_wilcoxauc(sobj, gene_list)

if (nrow(wilcoxauc_results) > 0) {
  # Save results
  output_file <- file.path(task_dir, "Documents", "wilcoxauc_results.csv")
  write.csv(wilcoxauc_results, output_file, row.names = FALSE)
  cat("Saved wilcoxauc results:", output_file, "\n")
  
  # Basic summary
  cat("\n=== Wilcoxauc Summary ===\n")
  summary_stats <- wilcoxauc_results %>%
    group_by(celltype) %>%
    summarise(
      genes_tested = n(),
      mean_auc = mean(auc, na.rm = TRUE),
      genes_auc_gt_0.6 = sum(auc > 0.6, na.rm = TRUE),
      genes_auc_gt_0.7 = sum(auc > 0.7, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_stats)
  
  # Save summary
  summary_file <- file.path(task_dir, "Documents", "wilcoxauc_summary.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# 4. SIMPLE PSEUDOBULK ANALYSIS
# ---------------------------------------------------------------------------
cat("\n4. Running pseudobulk analysis...\n")

# Function for pseudobulk analysis
run_pseudobulk <- function(seurat_obj, genes_of_interest, min_cells = 10) {
  
  cell_types <- unique(seurat_obj$celltype)
  results <- data.frame()
  
  for (ct in cell_types) {
    cat("  Processing", ct, "\n")
    
    # Subset cell type
    ct_obj <- subset(seurat_obj, subset = celltype == ct)
    
    # Check conditions and samples
    sample_counts <- table(ct_obj$sample, ct_obj$condition)
    if (nrow(sample_counts) < 2) {
      cat("    Skipping - insufficient samples\n")
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
      
      # Filter to genes of interest
      available_genes <- intersect(genes_of_interest, rownames(count_matrix))
      if (length(available_genes) < 5) {
        cat("    Skipping - insufficient genes\n")
        next
      }
      
      count_matrix <- count_matrix[available_genes, ]
      
      # Create DESeq2 dataset
      sample_meta <- data.frame(
        row.names = colnames(count_matrix),
        condition = ct_obj@meta.data$condition[match(colnames(count_matrix), ct_obj$sample)]
      )
      sample_meta <- sample_meta[!duplicated(rownames(sample_meta)), , drop = FALSE]
      count_matrix <- count_matrix[, rownames(sample_meta), drop = FALSE]
      
      # Run DESeq2
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = sample_meta,
        design = ~ condition
      )
      
      # Pre-filter and run
      keep <- rowSums(counts(dds) >= 10) >= 2
      dds <- dds[keep, ]
      dds <- DESeq(dds, quiet = TRUE)
      
      # Get results
      res <- results(dds, name = "condition_AD_vs_HC")
      
      # Format results
      ct_result <- data.frame(
        feature = rownames(res),
        gene = rownames(res),
        log2_fold_change = res$log2FoldChange,
        p_value = res$pvalue,
        padj = res$padj,
        auc_approx = 1 / (1 + exp(-res$log2FoldChange)),
        celltype = ct,
        method = "pseudobulk_DESeq2",
        stringsAsFactors = FALSE
      )
      
      # Remove NA rows
      ct_result <- ct_result[complete.cases(ct_result), ]
      
      if (nrow(ct_result) > 0) {
        results <- rbind(results, ct_result)
        cat("    Analyzed", nrow(ct_result), "genes\n")
      }
      
    }, error = function(e) {
      cat("    Error:", e$message, "\n")
    })
  }
  
  return(results)
}

# Run pseudobulk analysis (this might take a while)
pseudobulk_results <- run_pseudobulk(sobj, gene_list)

if (nrow(pseudobulk_results) > 0) {
  # Save results
  output_file <- file.path(task_dir, "Documents", "pseudobulk_results.csv")
  write.csv(pseudobulk_results, output_file, row.names = FALSE)
  cat("Saved pseudobulk results:", output_file, "\n")
  
  # Basic summary
  cat("\n=== Pseudobulk Summary ===\n")
  summary_stats <- pseudobulk_results %>%
    group_by(celltype) %>%
    summarise(
      genes_tested = n(),
      mean_auc = mean(auc_approx, na.rm = TRUE),
      significant = sum(padj < 0.05, na.rm = TRUE),
      lfc_gt_1 = sum(abs(log2_fold_change) > 1, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_stats)
  
  # Save summary
  summary_file <- file.path(task_dir, "Documents", "pseudobulk_summary.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# 5. COMPARE RESULTS IF BOTH AVAILABLE
# ---------------------------------------------------------------------------
if (nrow(wilcoxauc_results) > 0 && nrow(pseudobulk_results) > 0) {
  cat("\n5. Comparing methods...\n")
  
  # Merge results for comparison
  comparison <- wilcoxauc_results %>%
    select(feature, auc, celltype) %>%
    rename(wilcoxauc_auc = auc) %>%
    left_join(
      pseudobulk_results %>%
        select(feature, auc_approx, celltype) %>%
        rename(pseudobulk_auc = auc_approx),
      by = c("feature", "celltype")
    )
  
  # Remove rows with missing values
  comparison <- comparison[complete.cases(comparison), ]
  
  if (nrow(comparison) > 0) {
    # Calculate correlation
    cor_test <- cor.test(comparison$wilcoxauc_auc, comparison$pseudobulk_auc)
    
    cat("Correlation between methods:\n")
    cat("  r =", cor_test$estimate, "\n")
    cat("  p =", cor_test$p.value, "\n")
    cat("  n =", nrow(comparison), "\n")
    
    # Create scatter plot
    p <- ggplot(comparison, aes(x = wilcoxauc_auc, y = pseudobulk_auc)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      labs(
        title = "Method Comparison: Wilcoxauc vs Pseudobulk",
        subtitle = paste0("r = ", round(cor_test$estimate, 3), ", p = ", signif(cor_test$p.value, 3)),
        x = "AUC (Wilcoxauc)",
        y = "Approx. AUC (Pseudobulk)"
      ) +
      theme_minimal()
    
    # Save plot
    plot_file <- file.path(task_dir, "plot", "method_comparison.pdf")
    ggsave(plot_file, plot = p, width = 8, height = 6)
    cat("Saved comparison plot:", plot_file, "\n")
    
    # Save comparison data
    comparison_file <- file.path(task_dir, "Documents", "method_comparison.csv")
    write.csv(comparison, comparison_file, row.names = FALSE)
    cat("Saved comparison data:", comparison_file, "\n")
  }
}

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", task_dir, "\n")