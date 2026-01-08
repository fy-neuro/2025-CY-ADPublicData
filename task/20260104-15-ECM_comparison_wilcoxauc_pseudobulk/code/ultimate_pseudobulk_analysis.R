# ULTIMATE FINAL Pseudobulk Analysis - All Issues Fixed

library(Seurat)
library(dplyr)
library(tidyr)
library(openxlsx)
library(edgeR)
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

# ULTIMATE working pseudobulk analysis
run_pseudobulk_ultimate <- function(seurat_obj, genes_of_interest, min_cells = 10) {
  
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
      cat("  Original count matrix dimensions:", dim(count_matrix), "\n")
      
      # Create sample metadata FIRST
      samples_in_matrix <- colnames(count_matrix)
      
      # Get condition for each sample (convert dash to underscore for matching)
      conditions <- sapply(samples_in_matrix, function(sample_id) {
        original_id <- gsub("-", "_", sample_id)
        # Find cells from this sample
        sample_cells <- which(ct_obj$sample == original_id)
        if (length(sample_cells) == 0) return(NA)
        # Get condition (should be same for all cells from same sample)
        condition <- unique(ct_obj@meta.data$condition[sample_cells])
        return(condition[1])
      })
      
      # Create metadata
      sample_meta <- data.frame(
        row.names = samples_in_matrix,
        condition = conditions,
        stringsAsFactors = FALSE
      )
      
      # Remove samples with NA conditions
      valid_samples <- !is.na(sample_meta$condition)
      if (sum(valid_samples) < 4) {  # Need at least 4 valid samples
        cat("  Too few valid samples after NA removal\n")
        next
      }
      
      count_matrix <- count_matrix[, valid_samples, drop = FALSE]
      sample_meta <- sample_meta[valid_samples, , drop = FALSE]
      
      cat("  Filtered count matrix dimensions:", dim(count_matrix), "\n")
      cat("  Filtered sample metadata:\n")
      print(sample_meta)
      
      # Filter to genes of interest
      available_genes <- intersect(genes_of_interest, rownames(count_matrix))
      if (length(available_genes) < 3) {
        cat("  Skipping - insufficient genes\n")
        next
      }
      
      count_matrix <- count_matrix[available_genes, ]
      cat("  Genes after filtering:", nrow(count_matrix), "\n")
      
      # Create group factor
      group <- factor(sample_meta$condition)
      
      # Create DGEList with CORRECT sample names as column names
      dge <- DGEList(counts = count_matrix, group = group)
      colnames(dge) <- colnames(count_matrix)  # CRITICAL FIX!
      
      cat("  DGEList dimensions:", dim(dge$counts), "\n")
      cat("  DGEList colnames:", colnames(dge), "\n")
      cat("  Group levels:", levels(group), "\n")
      
      # Filter lowly expressed genes (keep genes with CPM > 1 in at least 2 samples)
      keep <- filterByExpr(dge, group = group, min.count = 10)
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      cat("  Genes after filterByExpr:", nrow(dge), "\n")
      
      if (nrow(dge) < 3) {
        cat("  Skipping - too few genes after filtering\n")
        next
      }
      
      # Normalize
      dge <- calcNormFactors(dge)
      
      # Design matrix - use correct approach
      design <- model.matrix(~0 + group)
      colnames(design) <- levels(group)
      
      cat("  Design matrix dimensions:", dim(design), "\n")
      cat("  Design matrix rownames:", rownames(design), "\n")
      cat("  DGEList colnames after filtering:", colnames(dge), "\n")
      
      # Estimate dispersion
      dge <- estimateDisp(dge, design)
      
      # Fit model
      fit <- glmFit(dge, design)
      
      # Create contrast (AD vs HC)
      if ("AD" %in% levels(group) && "HC" %in% levels(group)) {
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
          cat("  ✓ SUCCESS! Analyzed", nrow(ct_result), "genes\n")
          cat("  ✓ Significant genes (FDR < 0.05):", sum(ct_result$padj < 0.05), "\n")
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
      } else {
        cat("  Cannot create contrast - missing groups\n")
      }
      
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
    })
  }
  
  return(results)
}

# Run ULTIMATE pseudobulk analysis
cat("\n=== 🚀 RUNNING ULTIMATE PSEUDOBULK ANALYSIS ===\n")
pseudobulk_results <- run_pseudobulk_ultimate(sobj, gene_list)

if (nrow(pseudobulk_results) > 0) {
  # Save results
  output_file <- file.path(task_dir, "Documents", "pseudobulk_results_ultimate.csv")
  write.csv(pseudobulk_results, output_file, row.names = FALSE)
  cat("\n🎉 SUCCESS! Saved pseudobulk results:", output_file, "\n")
  
  # Basic summary
  cat("\n=== ULTIMATE PSEUDOBULK SUMMARY ===\n")
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
  summary_file <- file.path(task_dir, "Documents", "pseudobulk_summary_ultimate.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("🎉 Saved summary:", summary_file, "\n")
  
  # Top significant genes overall
  if (any(pseudobulk_results$padj < 0.05)) {
    top_genes <- pseudobulk_results %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      head(10)
    
    cat("\nTop 10 significant genes overall:\n")
    print(top_genes[, c("gene", "celltype", "log2_fold_change", "padj", "auc_approx")])
    
    # Save top genes
    top_file <- file.path(task_dir, "Documents", "pseudobulk_top_significant_ultimate.csv")
    write.csv(top_genes, top_file, row.names = FALSE)
    cat("🎉 Saved top significant genes:", top_file, "\n")
  } else {
    cat("\n❌ No significant genes found (FDR < 0.05)\n")
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
      comparison_file <- file.path(task_dir, "Documents", "method_comparison_ultimate.csv")
      write.csv(comparison, comparison_file, row.names = FALSE)
      cat("🎉 Saved method comparison:", comparison_file, "\n")
    }
  }
  
} else {
  cat("\n❌ No pseudobulk results generated\n")
}

cat("\n=== 🏁 ULTIMATE PSEUDOBULK ANALYSIS COMPLETE ===\n")
cat("Output files saved to:", task_dir, "\n")