# MINIMAL Pseudobulk Analysis - Workaround DESeq2 Issues

library(Seurat)
library(dplyr)
library(tidyr)
library(openxlsx)
library(edgeR)  # Try edgeR instead of DESeq2
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

# Try edgeR instead of DESeq2
run_pseudobulk_edger <- function(seurat_obj, genes_of_interest, min_cells = 10) {
  
  cell_types <- unique(seurat_obj$celltype)
  results <- data.frame()
  
  for (ct in cell_types) {
    cat("Processing", ct, "with edgeR...\n")
    
    # Subset cell type
    ct_obj <- subset(seurat_obj, subset = celltype == ct)
    
    # Check conditions and samples
    sample_counts <- table(ct_obj$sample, ct_obj$condition)
    cat("  Sample counts:\n")
    print(sample_counts)
    
    # Need at least 2 samples per condition
    ad_samples <- sum(sample_counts[, "AD"] > 0)
    hc_samples <- sum(sample_counts[, "HC"] > 0)
    
    if (ad_samples < 2 || hc_samples < 2) {
      cat("  Skipping - insufficient samples\n")
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
        cat("  Skipping - insufficient genes\n")
        next
      }
      
      count_matrix <- count_matrix[available_genes, ]
      cat("  Genes after filtering:", nrow(count_matrix), "\n")
      
      # Create sample metadata
      samples_in_matrix <- colnames(count_matrix)
      
      # Get condition for each sample
      conditions <- sapply(samples_in_matrix, function(sample_id) {
        # Convert back to underscore format for matching
        original_id <- gsub("-", "_", sample_id)
        condition <- unique(ct_obj@meta.data$condition[ct_obj@meta.data$sample == original_id])
        if (length(condition) == 0) return(NA)
        return(condition[1])
      })
      
      # Create metadata
      sample_meta <- data.frame(
        row.names = samples_in_matrix,
        condition = conditions,
        stringsAsFactors = FALSE
      )
      
      cat("  Sample metadata:\n")
      print(sample_meta)
      
      # Check for NA values
      if (any(is.na(sample_meta$condition))) {
        cat("  ERROR: NA values in condition metadata!\n")
        next
      }
      
      # Create DGEList for edgeR
      group <- factor(sample_meta$condition)
      dge <- DGEList(counts = count_matrix, group = group)
      
      # Filter lowly expressed genes
      keep <- filterByExpr(dge, group = group, min.count = 10)
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      cat("  Genes after filtering:", nrow(dge), "\n")
      
      if (nrow(dge) < 3) {
        cat("  Skipping - too few genes after filtering\n")
        next
      }
      
      # Normalize
      dge <- calcNormFactors(dge)
      
      # Design matrix
      design <- model.matrix(~0 + group)
      colnames(design) <- levels(group)
      
      # Estimate dispersion
      dge <- estimateDisp(dge, design)
      
      # Fit model
      fit <- glmFit(dge, design)
      
      # Contrast (AD vs HC)
      contrast <- makeContrasts(AD_vs_HC = AD - HC, levels = design)
      fit2 <- glmFit(dge, contrast)
      lrt <- glmLRT(fit2, contrast = contrast)
      
      # Get results
      top <- topTags(lrt, n = Inf)
      res <- top$table
      
      # Format results
      ct_result <- data.frame(
        gene = rownames(res),
        log2_fold_change = res$logFC,
        p_value = res$PValue,
        padj = res$FDR,
        logCPM = res$logCPM,
        auc_approx = 1 / (1 + exp(-res$logFC)),
        celltype = ct,
        method = "pseudobulk_edgeR",
        stringsAsFactors = FALSE
      )
      
      # Remove NA rows
      ct_result <- ct_result[complete.cases(ct_result), ]
      
      if (nrow(ct_result) > 0) {
        results <- rbind(results, ct_result)
        cat("  ✓ Analyzed", nrow(ct_result), "genes successfully\n")
        cat("  ✓ Significant genes (FDR < 0.05):", sum(ct_result$padj < 0.05), "\n")
        cat("  ✓ Genes with |log2FC| > 1:", sum(abs(ct_result$log2_fold_change) > 1), "\n")
        
        # Show top significant genes
        if (any(ct_result$padj < 0.05)) {
          top_sig <- ct_result %>% filter(padj < 0.05) %>% arrange(padj) %>% head(3)
          cat("  Top significant genes:\n")
          print(top_sig[, c("gene", "log2_fold_change", "padj")])
        }
      }
      
    }, error = function(e) {
      cat("  ERROR with edgeR:", e$message, "\n")
    })
  }
  
  return(results)
}

# Install edgeR if needed
if (!requireNamespace("edgeR", quietly = TRUE)) {
  cat("Installing edgeR package...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("edgeR", update = FALSE, ask = FALSE)
}

# Run edgeR pseudobulk analysis
cat("\n=== RUNNING EDGER PSEUDOBULK ANALYSIS ===\n")
pseudobulk_results <- run_pseudobulk_edger(sobj, gene_list)

if (nrow(pseudobulk_results) > 0) {
  # Save results
  output_file <- file.path(task_dir, "Documents", "pseudobulk_results_edger.csv")
  write.csv(pseudobulk_results, output_file, row.names = FALSE)
  cat("\n✅ Saved edgeR pseudobulk results:", output_file, "\n")
  
  # Basic summary
  cat("\n=== EDGER PSEUDOBULK SUMMARY ===\n")
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
  summary_file <- file.path(task_dir, "Documents", "pseudobulk_summary_edger.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("✅ Saved edgeR summary:", summary_file, "\n")
  
  # Top significant genes overall
  if (any(pseudobulk_results$padj < 0.05)) {
    top_genes <- pseudobulk_results %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      head(10)
    
    cat("\nTop 10 significant genes overall:\n")
    print(top_genes[, c("gene", "celltype", "log2_fold_change", "padj", "auc_approx")])
    
    # Save top genes
    top_file <- file.path(task_dir, "Documents", "pseudobulk_top_significant_edger.csv")
    write.csv(top_genes, top_file, row.names = FALSE)
    cat("✅ Saved top significant genes:", top_file, "\n")
  } else {
    cat("\n❌ No significant genes found (FDR < 0.05)\n")
  }
  
} else {
  cat("\n❌ No edgeR pseudobulk results generated\n")
}

cat("\n=== 🎉 EDGER PSEUDOBULK ANALYSIS COMPLETE ===\n")
cat("Output files saved to:", task_dir, "\n")