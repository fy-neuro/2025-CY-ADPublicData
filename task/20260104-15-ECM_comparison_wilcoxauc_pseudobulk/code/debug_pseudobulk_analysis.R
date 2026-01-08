# Pseudobulk Analysis with Detailed Logging

library(Seurat)
library(dplyr)
library(tidyr)
library(openxlsx)
library(edgeR)
library(here)

# Set working directory
setwd(here::here())

# Create log file
log_file <- "pseudobulk_debug_log.txt"
sink(log_file, split = TRUE, append = FALSE)
cat("=== PSEUDOBULK ANALYSIS DEBUG LOG ===\n")
cat("Timestamp:", format(Sys.time()), "\n\n")

# Load data
cat("Loading GSE188545 data...\n")
sobj <- readRDS("data/GSE188545/GSE188545_sobj_annotated_fast.rds")
cat("Loaded GSE188545 object:", ncol(sobj), "cells,", nrow(sobj), "genes\n")

# Join layers for Seurat v5
if ("Seurat" %in% class(sobj) && packageVersion("Seurat") >= "5.0.0") {
  if (length(Layers(sobj)) > 1) {
    sobj <- JoinLayers(sobj)
  }
}

# Load gene list
cat("Loading ECM gene list...\n")
gene_data <- read.xlsx("data/ECM_related_genes/ECM_communication_genesV3.xlsx", sheet = 1)
gene_list <- unique(gene_data[[1]][!is.na(gene_data[[1]])])
cat("Loaded", length(gene_list), "genes\n")

# Create output directory
task_dir <- "task/20260104-15-ECM_comparison_wilcoxauc_pseudobulk"
dir.create(file.path(task_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)

# Debug pseudobulk analysis with extensive logging
run_pseudobulk_debug <- function(seurat_obj, genes_of_interest, min_cells = 10) {
  
  cell_types <- unique(seurat_obj$celltype)
  results <- data.frame()
  
  for (ct in cell_types) {
    cat("\n", "="*80, "\n")
    cat("PROCESSING CELL TYPE:", ct, "\n")
    cat("="*80, "\n")
    
    # Subset cell type
    ct_obj <- subset(seurat_obj, subset = celltype == ct)
    cat("Subset created with", ncol(ct_obj), "cells\n")
    
    # Check conditions and samples
    sample_counts <- table(ct_obj$sample, ct_obj$condition)
    cat("Sample counts:\n")
    print(sample_counts)
    
    # Need at least 2 samples per condition
    ad_samples <- sum(sample_counts[, "AD"] > 0)
    hc_samples <- sum(sample_counts[, "HC"] > 0)
    cat("AD samples with cells:", ad_samples, "\n")
    cat("HC samples with cells:", hc_samples, "\n")
    
    if (ad_samples < 2 || hc_samples < 2) {
      cat("SKIPPING - insufficient samples (need ≥2 per condition)\n")
      next
    }
    
    tryCatch({
      # Step 1: Create pseudobulk counts
      cat("\n--- STEP 1: Creating pseudobulk counts ---\n")
      pseudobulk_counts <- AggregateExpression(
        ct_obj,
        assays = "RNA",
        slot = "counts",
        group.by = "sample",
        return.seurat = FALSE
      )
      
      count_matrix <- pseudobulk_counts$RNA
      cat("Original count matrix dimensions:", dim(count_matrix), "\n")
      cat("Original count matrix column names:\n")
      print(colnames(count_matrix))
      cat("Original count matrix row names (first 10):\n")
      print(head(rownames(count_matrix), 10))
      
      # Step 2: Create sample metadata
      cat("\n--- STEP 2: Creating sample metadata ---\n")
      samples_in_matrix <- colnames(count_matrix)
      cat("Samples in matrix:", paste(samples_in_matrix, collapse = ", "), "\n")
      
      # Get condition for each sample
      conditions <- sapply(samples_in_matrix, function(sample_id) {
        original_id <- gsub("-", "_", sample_id)
        cat("  Processing sample:", sample_id, "->", original_id, "\n")
        
        # Find cells from this sample
        sample_cells <- which(ct_obj$sample == original_id)
        cat("    Found", length(sample_cells), "cells from this sample\n")
        
        if (length(sample_cells) == 0) {
          cat("    ERROR: No cells found for sample\n")
          return(NA)
        }
        
        # Get condition (should be same for all cells from same sample)
        conditions_found <- unique(ct_obj@meta.data$condition[sample_cells])
        cat("    Conditions found:", paste(conditions_found, collapse = ", "), "\n")
        
        if (length(conditions_found) == 0) {
          cat("    ERROR: No conditions found\n")
          return(NA)
        }
        
        return(conditions_found[1])
      })
      
      cat("Extracted conditions for all samples:\n")
      print(data.frame(sample = samples_in_matrix, condition = conditions))
      
      # Create metadata
      sample_meta <- data.frame(
        row.names = samples_in_matrix,
        condition = conditions,
        stringsAsFactors = FALSE
      )
      
      cat("Sample metadata created:\n")
      print(sample_meta)
      
      # Remove samples with NA conditions
      valid_samples <- !is.na(sample_meta$condition)
      cat("Valid samples:", sum(valid_samples), "out of", length(valid_samples), "\n")
      
      if (sum(valid_samples) < 4) {
        cat("SKIPPING - too few valid samples after NA removal\n")
        next
      }
      
      count_matrix <- count_matrix[, valid_samples, drop = FALSE]
      sample_meta <- sample_meta[valid_samples, , drop = FALSE]
      
      cat("Filtered count matrix dimensions:", dim(count_matrix), "\n")
      cat("Filtered sample metadata:\n")
      print(sample_meta)
      
      # Step 3: Filter to genes of interest
      cat("\n--- STEP 3: Filtering to ECM genes ---\n")
      available_genes <- intersect(genes_of_interest, rownames(count_matrix))
      cat("Genes in list:", length(genes_of_interest), "\n")
      cat("Genes available in data:", length(available_genes), "\n")
      cat("Available genes:", paste(available_genes, collapse = ", "), "\n")
      
      if (length(available_genes) < 3) {
        cat("SKIPPING - insufficient genes\n")
        next
      }
      
      count_matrix <- count_matrix[available_genes, ]
      cat("Count matrix after gene filtering:", dim(count_matrix), "\n")
      
      # Step 4: Create DGEList
      cat("\n--- STEP 4: Creating DGEList ---\n")
      group <- factor(sample_meta$condition)
      cat("Group factor:\n")
      print(group)
      cat("Group levels:", levels(group), "\n")
      
      dge <- DGEList(counts = count_matrix, group = group)
      cat("Initial DGEList created\n")
      cat("DGEList counts dimensions:", dim(dge$counts), "\n")
      
      # CRITICAL: Set column names properly
      cat("Setting DGEList column names...\n")
      colnames(dge) <- colnames(count_matrix)
      cat("DGEList column names after setting:\n")
      print(colnames(dge))
      
      # Step 5: Filter lowly expressed genes
      cat("\n--- STEP 5: Filtering lowly expressed genes ---\n")
      keep <- filterByExpr(dge, group = group, min.count = 10)
      cat("Genes to keep:", sum(keep), "out of", length(keep), "\n")
      cat("Genes being filtered:", paste(rownames(dge)[!keep], collapse = ", "), "\n")
      
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      cat("DGEList after filtering:", dim(dge), "\n")
      cat("DGEList column names after filtering:\n")
      print(colnames(dge))
      
      if (nrow(dge) < 3) {
        cat("SKIPPING - too few genes after filtering\n")
        next
      }
      
      # Step 6: Normalize
      cat("\n--- STEP 6: Normalization ---\n")
      dge <- calcNormFactors(dge)
      cat("Normalization completed\n")
      
      # Step 7: Design matrix
      cat("\n--- STEP 7: Creating design matrix ---\n")
      design <- model.matrix(~0 + group)
      colnames(design) <- levels(group)
      
      cat("Design matrix created:\n")
      cat("Dimensions:", dim(design), "\n")
      cat("Row names:", rownames(design), "\n")
      cat("Column names:", colnames(design), "\n")
      print(design)
      
      # Verify dimensions match
      cat("Dimension check:\n")
      cat("  Design matrix rows:", nrow(design), "\n")
      cat("  Design matrix cols:", ncol(design), "\n")
      cat("  DGEList counts rows:", nrow(dge$counts), "\n")
      cat("  DGEList counts cols:", ncol(dge$counts), "\n")
      cat("  DGEList column names:", colnames(dge), "\n")
      
      if (nrow(design) != ncol(dge$counts)) {
        cat("ERROR: Design matrix rows don't match DGEList columns!\n")
        stop("Dimension mismatch")
      }
      
      # Step 8: Estimate dispersion
      cat("\n--- STEP 8: Estimating dispersion ---\n")
      dge <- estimateDisp(dge, design)
      cat("Dispersion estimation completed\n")
      
      # Step 9: Fit model
      cat("\n--- STEP 9: Fitting model ---\n")
      fit <- glmFit(dge, design)
      cat("Model fitting completed\n")
      
      # Step 10: Create contrast and test
      cat("\n--- STEP 10: Differential expression testing ---\n")
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
          cat("SUCCESS! Analyzed", nrow(ct_result), "genes\n")
          cat("Significant genes (FDR < 0.05):", sum(ct_result$padj < 0.05), "\n")
          cat("Genes with |log2FC| > 1:", sum(abs(ct_result$log2_fold_change) > 1), "\n")
          
          # Show top significant genes
          if (any(ct_result$padj < 0.05)) {
            top_sig <- ct_result %>% filter(padj < 0.05) %>% arrange(padj) %>% head(3)
            cat("Top significant genes:\n")
            print(top_sig[, c("gene", "log2_fold_change", "padj")])
          }
        } else {
          cat("No valid results after filtering\n")
        }
      } else {
        cat("ERROR: Cannot create contrast - missing groups\n")
        cat("Available groups:", levels(group), "\n")
      }
      
    }, error = function(e) {
      cat("ERROR CATCHED:", e$message, "\n")
      cat("Error call stack:\n")
      print(traceback())
    })
  }
  
  return(results)
}

# Run debug analysis
cat("\n", "="*80, "\n")
cat("STARTING DEBUG PSEUDOBULK ANALYSIS\n")
cat("="*80, "\n")

pseudobulk_results <- run_pseudobulk_debug(sobj, gene_list)

if (nrow(pseudobulk_results) > 0) {
  # Save results
  output_file <- file.path(task_dir, "Documents", "pseudobulk_results_debug.csv")
  write.csv(pseudobulk_results, output_file, row.names = FALSE)
  cat("\nSUCCESS! Saved pseudobulk results:", output_file, "\n")
} else {
  cat("\nFAILED: No pseudobulk results generated\n")
}

cat("\n", "="*80, "\n")
cat("DEBUG ANALYSIS COMPLETE\n")
cat("="*80, "\n")
cat("Final timestamp:", format(Sys.time()), "\n")

# Close log file
sink()

cat("Debug log saved to:", log_file, "\n")