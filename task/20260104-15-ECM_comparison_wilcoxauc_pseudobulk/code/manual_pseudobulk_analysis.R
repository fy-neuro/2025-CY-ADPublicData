# MANUAL PSEUDOBULK ANALYSIS - Working Solution

library(Seurat)
library(dplyr)
library(tidyr)
library(openxlsx)
library(edgeR)
library(here)

# Set working directory
setwd(here::here())

# Load data
cat("Loading GSE188545 data...\n")
sobj <- readRDS("data/GSE188545/GSE188545_sobj_annotated_fast.rds")
cat("Loaded GSE188545 object:", ncol(sobj), "cells,", nrow(sobj), "genes\n")

# Check Seurat structure - we have separate layers for each sample!
cat("Available layers:\n")
print(head(Layers(sobj[['RNA'])), 10))

# Load gene list
cat("Loading ECM gene list...\n")
gene_data <- read.xlsx("data/ECM_related_genes/ECM_communication_genesV3.xlsx", sheet = 1)
gene_list <- unique(gene_data[[1]][!is.na(gene_data[[1]])])
cat("Loaded", length(gene_list), "genes\n")

# Create output directory
task_dir <- "task/20260104-15-ECM_comparison_wilcoxauc_pseudobulk"
dir.create(file.path(task_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)

# MANUAL pseudobulk aggregation
run_manual_pseudobulk <- function(seurat_obj, genes_of_interest, min_cells = 10) {
  
  cell_types <- unique(seurat_obj$celltype)
  results <- data.frame()
  
  for (ct in cell_types) {
    cat("Processing", ct, "...\n")
    
    # Subset cell type
    ct_obj <- subset(seurat_obj, subset = celltype == ct)
    cat("  Subset created with", ncol(ct_obj), "cells\n")
    
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
      # Manual pseudobulk aggregation
      cat("  Creating manual pseudobulk counts...\n")
      
      # Get unique samples in this cell type
      unique_samples <- unique(ct_obj$sample)
      
      # Get all ECM genes that exist in the data
      available_genes <- intersect(genes_of_interest, rownames(ct_obj))
      
      if (length(available_genes) < 3) {
        cat("  Skipping - insufficient available genes\n")
        next
      }
      
      cat("  Available ECM genes:", length(available_genes), "\n")
      cat("  Samples to aggregate:", length(unique_samples), "\n")
      
      # Initialize pseudobulk matrix
      pseudobulk_matrix <- matrix(0, nrow = length(available_genes), 
                               ncol = length(unique_samples))
      rownames(pseudobulk_matrix) <- available_genes
      colnames(pseudobulk_matrix) <- unique_samples
      
      # Aggregate counts manually
      for (i in seq_along(unique_samples)) {
        sample_id <- unique_samples[i]
        cat("    Processing sample:", sample_id, "\n")
        
        # Find cells from this sample
        sample_cells <- which(ct_obj$sample == sample_id)
        
        if (length(sample_cells) == 0) {
          cat("      No cells found\n")
          next
        }
        
        cat("      Found", length(sample_cells), "cells\n")
        
        # Extract counts from the appropriate layer
        # The layer name is "counts.SAMPLEID"
        layer_name <- paste0("counts.", sample_id)
        
        if (layer_name %in% Layers(ct_obj[['RNA']])) {
          # Get counts for this sample's layer
          sample_counts <- GetAssayData(ct_obj, assay = "RNA", layer = layer_name)
          
          # Filter to available ECM genes
          sample_counts_ecm <- sample_counts[available_genes, , drop = FALSE]
          
          # Sum across all cells from this sample
          pseudobulk_matrix[, i] <- rowSums(sample_counts_ecm)
          
          cat("      Aggregated", sum(pseudobulk_matrix[, i] > 0), "expressed genes\n")
        } else {
          cat("      WARNING: Layer", layer_name, "not found!\n")
        }
      }
      
      cat("  Final pseudobulk matrix dimensions:", dim(pseudobulk_matrix), "\n")
      cat("  Total counts per sample:\n")
      print(colSums(pseudobulk_matrix))
      
      # Create sample metadata
      sample_conditions <- sapply(colnames(pseudobulk_matrix), function(sample_id) {
        sample_cells <- which(ct_obj$sample == sample_id)
        if (length(sample_cells) == 0) return(NA)
        condition <- unique(ct_obj@meta.data$condition[sample_cells])
        return(condition[1])
      })
      
      sample_meta <- data.frame(
        row.names = colnames(pseudobulk_matrix),
        condition = sample_conditions,
        stringsAsFactors = FALSE
      )
      
      cat("  Sample metadata:\n")
      print(sample_meta)
      
      # Remove samples with NA conditions or zero counts
      valid_samples <- !is.na(sample_meta$condition) & colSums(pseudobulk_matrix) > 0
      
      if (sum(valid_samples) < 4) {
        cat("  Too few valid samples after filtering\n")
        next
      }
      
      pseudobulk_matrix <- pseudobulk_matrix[, valid_samples, drop = FALSE]
      sample_meta <- sample_meta[valid_samples, , drop = FALSE]
      
      cat("  After filtering:", dim(pseudobulk_matrix), "\n")
      
      # Create group factor
      group <- factor(sample_meta$condition)
      
      # Create DGEList
      dge <- DGEList(counts = pseudobulk_matrix, group = group)
      
      # Filter lowly expressed genes
      keep <- filterByExpr(dge, group = group, min.count = 10)
      cat("  Genes to keep:", sum(keep), "out of", length(keep), "\n")
      
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      
      if (nrow(dge) < 3) {
        cat("  Skipping - too few genes after filtering\n")
        next
      }
      
      cat("  Final DGEList dimensions:", dim(dge), "\n")
      
      # Normalize
      dge <- calcNormFactors(dge)
      cat("  Normalization completed\n")
      
      # Create design matrix manually
      samples <- colnames(dge)
      n_samples <- length(samples)
      design <- matrix(0, nrow = n_samples, ncol = 2)
      rownames(design) <- samples
      colnames(design) <- c("AD", "HC")
      
      # Fill design matrix
      ad_samples_idx <- which(group == "AD")
      hc_samples_idx <- which(group == "HC")
      design[ad_samples_idx, "AD"] <- 1
      design[hc_samples_idx, "HC"] <- 1
      
      cat("  Design matrix:\n")
      print(design)
      
      # Estimate dispersion
      dge <- estimateDisp(dge, design)
      cat("  Dispersion estimation completed\n")
      
      # Fit model
      fit <- glmFit(dge, design)
      cat("  Model fitting completed\n")
      
      # Create contrast (AD vs HC)
      contrast <- c(1, -1)  # AD - HC
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
        method = "manual_pseudobulk",
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
      
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
      cat("  Error call stack:\n")
      print(traceback())
    })
  }
  
  return(results)
}

# Run MANUAL pseudobulk analysis
cat("\n=== 🚀 RUNNING MANUAL PSEUDOBULK ANALYSIS ===\n")
pseudobulk_results <- run_manual_pseudobulk(sobj, gene_list)

if (nrow(pseudobulk_results) > 0) {
  # Save results
  output_file <- file.path(task_dir, "Documents", "pseudobulk_results_manual.csv")
  write.csv(pseudobulk_results, output_file, row.names = FALSE)
  cat("\n🎉 SUCCESS! Saved manual pseudobulk results:", output_file, "\n")
  
  # Basic summary
  cat("\n=== MANUAL PSEUDOBULK SUMMARY ===\n")
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
  summary_file <- file.path(task_dir, "Documents", "pseudobulk_summary_manual.csv")
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
    top_file <- file.path(task_dir, "Documents", "pseudobulk_top_significant_manual.csv")
    write.csv(top_genes, top_file, row.names = FALSE)
    cat("🎉 Saved top significant genes:", top_file, "\n")
  } else {
    cat("\n❌ No significant genes found (FDR < 0.05)\n")
  }
  
  # Compare with wilcoxauc results if available
  wilcoxauc_file <- file.path(task_dir, "Documents", "wilcoxauc_results.csv")
  if (file.exists(wilcoxauc_file)) {
    cat("\n=== COMPARING MANUAL PSEUDOBULK VS WILCOXAUC ===\n")
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
      comparison_file <- file.path(task_dir, "Documents", "method_comparison_manual.csv")
      write.csv(comparison, comparison_file, row.names = FALSE)
      cat("🎉 Saved method comparison:", comparison_file, "\n")
    }
  }
  
} else {
  cat("\n❌ No manual pseudobulk results generated\n")
}

cat("\n=== 🏁 MANUAL PSEUDOBULK ANALYSIS COMPLETE ===\n")
cat("Output files saved to:", task_dir, "\n")