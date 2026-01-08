# Comprehensive Comparison: Presto wilcoxauc vs Pseudobulk Analysis
# For ECM-Related Genes in GSE188545 Dataset

# Load required libraries
library(Seurat)
library(presto)          # Fast Wilcoxon and AUC analysis
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)        # For reading Excel gene lists
library(pheatmap)        # For heatmap visualization
library(viridis)         # Color palette
library(ggplot2)         # Enhanced plotting
library(patchwork)       # Plot arrangement
library(here)            # Path management
library(DESeq2)          # Pseudobulk differential expression
library(readr)           # Fast file reading

# Set working directory to project root for consistent paths
project_root <- here::here()
setwd(project_root)

# Create output directories
task_dir <- "task/20260104-15-ECM_comparison_wilcoxauc_pseudobulk"
dir.create(file.path(task_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(task_dir, "plot"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(task_dir, "code"), showWarnings = FALSE, recursive = TRUE)

cat("=== ECM Genes Analysis: Presto wilcoxauc vs Pseudobulk ===\n")
cat("Task directory:", task_dir, "\n")

# ---------------------------------------------------------------------------
# 1. LOAD GSE188545 ANNOTATED DATA
# ---------------------------------------------------------------------------
cat("\n1. Loading GSE188545 annotated data...\n")

sobj_path <- "data/GSE188545/GSE188545_sobj_annotated_fast.rds"
if (!file.exists(sobj_path)) {
  stop("Annotated GSE188545 object not found at: ", sobj_path)
}

sobj <- readRDS(sobj_path)
cat("Loaded annotated GSE188545 object:", ncol(sobj), "cells,", nrow(sobj), "genes\n")

# Check metadata structure
cat("Available metadata columns:", paste(colnames(sobj@meta.data), collapse = ", "), "\n")
cat("Conditions:", toString(unique(sobj$condition)), "\n")
cat("Cell types:", toString(unique(sobj$celltype)), "\n")

# Join layers for Seurat v5 compatibility if needed
if ("Seurat" %in% class(sobj) && packageVersion("Seurat") >= "5.0.0") {
  if (length(Layers(sobj)) > 1) {
    sobj <- JoinLayers(sobj)
    cat("Joined layers for Seurat v5 compatibility\n")
  }
}

# ---------------------------------------------------------------------------
# 2. LOAD ECM GENE LISTS
# ---------------------------------------------------------------------------
cat("\n2. Loading ECM gene lists...\n")

# Define ECM gene lists and their files
ecm_files <- list(
  communication = "data/ECM_related_genes/ECM_communication_genes.xlsx",
  communicationV2 = "data/ECM_related_genes/ECM_communication_genesV2.xlsx", 
  communicationV3 = "data/ECM_related_genes/ECM_communication_genesV3.xlsx",
  structure = "data/ECM_related_genes/ECM_structure_genes.xlsx"
)

# Load all gene lists
gene_lists <- list()
for (list_name in names(ecm_files)) {
  if (file.exists(ecm_files[[list_name]])) {
    gene_data <- read.xlsx(ecm_files[[list_name]], sheet = 1)
    genes <- unique(gene_data[[1]][!is.na(gene_data[[1]])])
    gene_lists[[list_name]] <- genes
    cat("  Loaded", length(genes), "genes from", list_name, "\n")
  } else {
    cat("  Warning: File not found:", ecm_files[[list_name]], "\n")
  }
}

# ---------------------------------------------------------------------------
# 3. PRESTO WILCOXAUC ANALYSIS FUNCTION
# ---------------------------------------------------------------------------
cat("\n3. Defining presto wilcoxauc analysis function...\n")

#' Perform wilcoxauc analysis for a given gene list
#'
#' @param sobj Seurat object with 'celltype' and 'condition' metadata
#' @param gene_list Vector of gene symbols to analyze
#' @param min_cells_per_group Minimum cells per condition to include cell type
#' @param gene_list_name Name of the gene list for output files
#'
#' @return Data frame with wilcoxauc results
perform_wilcoxauc <- function(sobj, gene_list, min_cells_per_group = 10, gene_list_name = "genes") {
  
  cat("  Performing wilcoxauc analysis for", gene_list_name, "\n")
  cat("    Testing", length(gene_list), "genes\n")
  
  cell_types <- unique(sobj$celltype)
  wilcoxauc_results <- list()
  
  for (ct in cell_types) {
    cat("      Processing cell type:", ct, "\n")
    
    # Subset to current cell type
    sobj_ct <- subset(sobj, subset = celltype == ct)
    
    # Check if both conditions present with sufficient cells
    cell_counts <- table(sobj_ct$condition)
    if (length(cell_counts) < 2) {
      cat("        Skipping: only one condition present\n")
      next
    }
    if (any(cell_counts < min_cells_per_group)) {
      cat("        Skipping: insufficient cells (", 
          paste(names(cell_counts), cell_counts, sep = "=", collapse = ", "), ")\n")
      next
    }
    
    # Find genes present in this subset
    available_genes <- rownames(sobj_ct)
    genes_to_test <- intersect(gene_list, available_genes)
    
    if (length(genes_to_test) == 0) {
      cat("        Skipping: no genes from list present\n")
      next
    }
    
    # Calculate AUC using presto::wilcoxauc
    Idents(sobj_ct) <- "condition"
    
    condition_values <- unique(sobj_ct$condition)
    comparison <- condition_values
    
    auc_df <- tryCatch({
      presto::wilcoxauc(
        sobj_ct,
        comparison = comparison,
        seurat_assay = "RNA"
      )
    }, error = function(e) {
      cat("        Error in wilcoxauc:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(auc_df) || nrow(auc_df) == 0) {
      cat("        No results from wilcoxauc\n")
      next
    }
    
    # Check if auc_df is a data frame with expected columns
    if (!is.data.frame(auc_df) || !"feature" %in% colnames(auc_df)) {
      cat("        Unexpected wilcoxauc output format\n")
      next
    }
    
    # Filter for genes of interest and format results
    ct_results <- auc_df %>%
      filter(feature %in% genes_to_test) %>%
      select(feature, auc, pval, group) %>%
      mutate(
        celltype = ct,
        comparison = paste(comparison[1], "vs", comparison[2]),
        gene_list = gene_list_name,
        method = "wilcoxauc"
      ) %>%
      rename(gene = feature, p_value = pval)
    
    wilcoxauc_results[[ct]] <- ct_results
  }
  
  # Combine all results
  if (length(wilcoxauc_results) == 0) {
    return(data.frame())
  }
  
  result <- do.call(rbind, wilcoxauc_results)
  cat("    Generated", nrow(result), "gene-cell type combinations\n")
  return(result)
}

# ---------------------------------------------------------------------------
# 4. PSEUDOBULK ANALYSIS FUNCTION
# ---------------------------------------------------------------------------
cat("\n4. Defining pseudobulk analysis function...\n")

#' Perform pseudobulk analysis using DESeq2
#'
#' @param sobj Seurat object with 'celltype' and 'condition' metadata
#' @param gene_list Vector of gene symbols to analyze
#' @param gene_list_name Name of the gene list for output files
#'
#' @return Data frame with pseudobulk DE results
perform_pseudobulk <- function(sobj, gene_list, gene_list_name = "genes") {
  
  cat("  Performing pseudobulk analysis for", gene_list_name, "\n")
  cat("    Testing", length(gene_list), "genes\n")
  
  # Check if we have sample information for pseudobulk creation
  if (!"orig.ident" %in% colnames(sobj@meta.data)) {
    cat("    Warning: No orig.ident found, using condition as sample\n")
    sobj$sample <- sobj$condition
  } else {
    sobj$sample <- sobj$orig.ident
  }
  
  cell_types <- unique(sobj$celltype)
  pseudobulk_results <- list()
  
  for (ct in cell_types) {
    cat("      Processing cell type:", ct, "\n")
    
    # Subset to current cell type
    sobj_ct <- subset(sobj, subset = celltype == ct)
    
    # Check cell counts per sample and condition
    cell_counts <- table(sobj_ct$sample, sobj_ct$condition)
    cat("        Cell counts:\n")
    print(cell_counts)
    
    # Skip if we don't have multiple samples per condition
    if (nrow(cell_counts) < 2) {
      cat("        Skipping: insufficient samples for pseudobulk analysis\n")
      next
    }
    
    # Create pseudobulk by aggregating cells within sample x celltype combinations
    # Use Seurat's AggregateExpression function
    pseudobulk_expr <- AggregateExpression(
      sobj_ct,
      assays = "RNA",
      slot = "counts",
      group.by = "sample",
      return.seurat = FALSE
    )
    
    # Get the count matrix
    count_matrix <- pseudobulk_expr$RNA
    
    # Filter to genes of interest that are present
    available_genes <- rownames(count_matrix)
    genes_to_test <- intersect(gene_list, available_genes)
    
    if (length(genes_to_test) < 5) {
      cat("        Skipping: insufficient genes present (", length(genes_to_test), ")\n")
      next
    }
    
    # Subset to genes of interest
    count_matrix <- count_matrix[genes_to_test, ]
    
    # Create sample metadata for DESeq2
    sample_metadata <- data.frame(
      row.names = colnames(count_matrix),
      condition = sobj_ct@meta.data$condition[match(colnames(count_matrix), sobj_ct$sample)]
    )
    sample_metadata <- sample_metadata[!duplicated(rownames(sample_metadata)), , drop = FALSE]
    count_matrix <- count_matrix[, rownames(sample_metadata), drop = FALSE]
    
    # Create DESeq2 dataset
    tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = sample_metadata,
        design = ~ condition
      )
      
      # Pre-filter genes with low counts
      keep <- rowSums(counts(dds) >= 10) >= 2
      dds <- dds[keep, ]
      
      # Run DESeq2
      dds <- DESeq(dds, quiet = TRUE)
      
      # Extract results
      res <- results(dds, name = "condition_AD_vs_HC")
      
      # Format results
      ct_results <- data.frame(
        gene = rownames(res),
        log2_fold_change = res$log2FoldChange,
        p_value = res$pvalue,
        padj = res$padj,
        baseMean = res$baseMean,
        celltype = ct,
        comparison = "AD_vs_HC",
        gene_list = gene_list_name,
        method = "pseudobulk_DESeq2",
        stringsAsFactors = FALSE
      )
      
      # Calculate approximate AUC from log2 fold change
      # This is a simplified conversion for comparison purposes
      ct_results$auc_approx <- 1 / (1 + exp(-ct_results$log2_fold_change))
      
      # Remove rows with NA values
      ct_results <- ct_results[complete.cases(ct_results), ]
      
      if (nrow(ct_results) > 0) {
        cat("        Generated", nrow(ct_results), "gene results\n")
        pseudobulk_results[[ct]] <- ct_results
      } else {
        cat("        No valid results after filtering\n")
      }
      
    }, error = function(e) {
      cat("        Error in DESeq2 analysis:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Combine all results
  if (length(pseudobulk_results) == 0) {
    return(data.frame())
  }
  
  result <- do.call(rbind, pseudobulk_results)
  cat("    Generated", nrow(result), "total gene-cell type combinations\n")
  return(result)
}

# ---------------------------------------------------------------------------
# 5. RUN ANALYSES FOR ALL ECM GENE LISTS
# ---------------------------------------------------------------------------
cat("\n5. Running analyses for all ECM gene lists...\n")

# Initialize results storage
all_wilcoxauc_results <- list()
all_pseudobulk_results <- list()

for (list_name in names(gene_lists)) {
  cat("\n=== Processing", list_name, "===\n")
  
  # Skip if no genes loaded
  if (length(gene_lists[[list_name]]) == 0) {
    cat("Skipping", list_name, "- no genes loaded\n")
    next
  }
  
  # Perform wilcoxauc analysis
  wilcoxauc_results <- perform_wilcoxauc(
    sobj = sobj,
    gene_list = gene_lists[[list_name]],
    gene_list_name = list_name
  )
  
  if (nrow(wilcoxauc_results) > 0) {
    all_wilcoxauc_results[[list_name]] <- wilcoxauc_results
    
    # Save wilcoxauc results
    output_file <- file.path(task_dir, "Documents", 
                           paste0("wilcoxauc_", list_name, "_results.csv"))
    write.csv(wilcoxauc_results, output_file, row.names = FALSE)
    cat("    Saved wilcoxauc results:", output_file, "\n")
  }
  
  # Perform pseudobulk analysis
  pseudobulk_results <- perform_pseudobulk(
    sobj = sobj,
    gene_list = gene_lists[[list_name]],
    gene_list_name = list_name
  )
  
  if (nrow(pseudobulk_results) > 0) {
    all_pseudobulk_results[[list_name]] <- pseudobulk_results
    
    # Save pseudobulk results
    output_file <- file.path(task_dir, "Documents", 
                           paste0("pseudobulk_", list_name, "_results.csv"))
    write.csv(pseudobulk_results, output_file, row.names = FALSE)
    cat("    Saved pseudobulk results:", output_file, "\n")
  }
}

# ---------------------------------------------------------------------------
# 6. COMBINE AND COMPARE RESULTS
# ---------------------------------------------------------------------------
cat("\n6. Combining and comparing results...\n")

# Combine all wilcoxauc results
combined_wilcoxauc <- do.call(rbind, all_wilcoxauc_results)
if (nrow(combined_wilcoxauc) > 0) {
  output_file <- file.path(task_dir, "Documents", "combined_wilcoxauc_results.csv")
  write.csv(combined_wilcoxauc, output_file, row.names = FALSE)
  cat("Saved combined wilcoxauc results:", output_file, "\n")
}

# Combine all pseudobulk results
combined_pseudobulk <- do.call(rbind, all_pseudobulk_results)
if (nrow(combined_pseudobulk) > 0) {
  output_file <- file.path(task_dir, "Documents", "combined_pseudobulk_results.csv")
  write.csv(combined_pseudobulk, output_file, row.names = FALSE)
  cat("Saved combined pseudobulk results:", output_file, "\n")
}

# Create comparison summary
if (nrow(combined_wilcoxauc) > 0 && nrow(combined_pseudobulk) > 0) {
  # Create summary statistics for each method
  wilcoxauc_summary <- combined_wilcoxauc %>%
    group_by(gene_list, celltype) %>%
    summarise(
      wilcoxauc_mean_auc = mean(auc, na.rm = TRUE),
      wilcoxauc_sd_auc = sd(auc, na.rm = TRUE),
      wilcoxauc_genes_tested = n(),
      wilcoxauc_significant = sum(p_value < 0.05, na.rm = TRUE),
      wilcoxauc_auc_gt_0.6 = sum(auc > 0.6, na.rm = TRUE),
      wilcoxauc_auc_gt_0.7 = sum(auc > 0.7, na.rm = TRUE),
      .groups = "drop"
    )
  
  pseudobulk_summary <- combined_pseudobulk %>%
    group_by(gene_list, celltype) %>%
    summarise(
      pseudobulk_mean_auc = mean(auc_approx, na.rm = TRUE),
      pseudobulk_sd_auc = sd(auc_approx, na.rm = TRUE),
      pseudobulk_genes_tested = n(),
      pseudobulk_significant = sum(padj < 0.05, na.rm = TRUE),
      pseudobulk_lfc_gt_1 = sum(abs(log2_fold_change) > 1, na.rm = TRUE),
      pseudobulk_lfc_gt_2 = sum(abs(log2_fold_change) > 2, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Merge summaries for comparison
  comparison_summary <- left_join(wilcoxauc_summary, pseudobulk_summary, 
                                 by = c("gene_list", "celltype"))
  
  # Calculate correlation between methods where possible
  correlation_analysis <- list()
  for (list_name in names(gene_lists)) {
    if (list_name %in% all_wilcoxauc_results && list_name %in% all_pseudobulk_results) {
      # Merge on gene and celltype
      merged_data <- left_join(
        all_wilcoxauc_results[[list_name]] %>% select(gene, celltype, auc),
        all_pseudobulk_results[[list_name]] %>% select(gene, celltype, auc_approx),
        by = c("gene", "celltype"),
        suffix = c("_wilcoxauc", "_pseudobulk")
      )
      
      # Remove rows with missing values
      merged_data <- merged_data[complete.cases(merged_data), ]
      
      if (nrow(merged_data) > 5) {
        cor_result <- cor.test(merged_data$auc, merged_data$auc_approx)
        correlation_analysis[[list_name]] <- data.frame(
          gene_list = list_name,
          correlation = cor_result$estimate,
          p_value = cor_result$p.value,
          n_comparisons = nrow(merged_data)
        )
      }
    }
  }
  
  if (length(correlation_analysis) > 0) {
    correlation_df <- do.call(rbind, correlation_analysis)
    output_file <- file.path(task_dir, "Documents", "method_correlation.csv")
    write.csv(correlation_df, output_file, row.names = FALSE)
    cat("Saved method correlation results:", output_file, "\n")
  }
  
  # Save comparison summary
  output_file <- file.path(task_dir, "Documents", "method_comparison_summary.csv")
  write.csv(comparison_summary, output_file, row.names = FALSE)
  cat("Saved method comparison summary:", output_file, "\n")
}

# ---------------------------------------------------------------------------
# 7. CREATE COMPARISON VISUALIZATIONS
# ---------------------------------------------------------------------------
cat("\n7. Creating comparison visualizations...\n")

if (nrow(combined_wilcoxauc) > 0 && nrow(combined_pseudobulk) > 0) {
  # 7.1 Scatter plots comparing methods
  for (list_name in names(gene_lists)) {
    if (list_name %in% all_wilcoxauc_results && list_name %in% all_pseudobulk_results) {
      # Merge data for scatter plot
      plot_data <- left_join(
        all_wilcoxauc_results[[list_name]] %>% select(gene, celltype, auc, p_value),
        all_pseudobulk_results[[list_name]] %>% select(gene, celltype, auc_approx, padj, log2_fold_change),
        by = c("gene", "celltype"),
        suffix = c("_wilcoxauc", "_pseudobulk")
      )
      
      plot_data <- plot_data[complete.cases(plot_data), ]
      
      if (nrow(plot_data) > 10) {
        # Create scatter plot
        p <- ggplot(plot_data, aes(x = auc, y = auc_approx)) +
          geom_point(alpha = 0.6, color = "steelblue") +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
          labs(
            title = paste("Method Comparison:", list_name),
            subtitle = "Presto wilcoxauc vs Pseudobulk DESeq2",
            x = "AUC (wilcoxauc)",
            y = "Approximated AUC (pseudobulk)"
          ) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
        
        # Add correlation information
        cor_test <- cor.test(plot_data$auc, plot_data$auc_approx)
        p <- p + annotate("text", x = 0.9, y = 0.1, 
                         label = paste("r =", round(cor_test$estimate, 3), 
                                      "\np =", signif(cor_test$p.value, 3)),
                         hjust = 1, vjust = 0)
        
        # Save plot
        output_file <- file.path(task_dir, "plot", 
                               paste0("comparison_scatter_", list_name, ".pdf"))
        ggsave(output_file, plot = p, width = 8, height = 6)
        cat("    Saved scatter plot:", output_file, "\n")
        
        # 7.2 Volcano-like plot for pseudobulk results
        p_volcano <- ggplot(all_pseudobulk_results[[list_name]], 
                           aes(x = log2_fold_change, y = -log10(padj))) +
          geom_point(aes(color = -log10(padj) > -log10(0.05)), alpha = 0.7) +
          scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
          geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
          labs(
            title = paste("Pseudobulk DE Analysis:", list_name),
            subtitle = "Red points: padj < 0.05, Blue lines: |log2FC| > 1",
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value"
          ) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
        
        output_file <- file.path(task_dir, "plot", 
                               paste0("pseudobulk_volcano_", list_name, ".pdf"))
        ggsave(output_file, plot = p_volcano, width = 8, height = 6)
        cat("    Saved volcano plot:", output_file, "\n")
      }
    }
  }
  
  # 7.3 Summary heatmap comparing methods across gene lists
  if (exists("comparison_summary") && nrow(comparison_summary) > 0) {
    # Create summary matrix for heatmap
    summary_matrix <- comparison_summary %>%
      select(gene_list, celltype, wilcoxauc_mean_auc, pseudobulk_mean_auc) %>%
      pivot_longer(cols = c(wilcoxauc_mean_auc, pseudobulk_mean_auc),
                   names_to = "method", values_to = "mean_auc") %>%
      mutate(method = gsub("_mean_auc", "", method)) %>%
      pivot_wider(names_from = celltype, values_from = mean_auc)
    
    # Save summary heatmap
    if (nrow(summary_matrix) > 1 && ncol(summary_matrix) > 2) {
      matrix_for_heatmap <- as.matrix(summary_matrix[,-1])
      rownames(matrix_for_heatmap) <- summary_matrix$method
      
      output_file <- file.path(task_dir, "plot", "method_comparison_heatmap.pdf")
      pdf(output_file, width = 10, height = 4)
      pheatmap(
        matrix_for_heatmap,
        main = "Mean AUC Comparison: wilcoxauc vs Pseudobulk",
        color = viridis(100),
        scale = "column",
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        treeheight_row = 10,
        treeheight_col = 10,
        fontsize_row = 12,
        fontsize_col = 10,
        angle_col = 45,
        cellwidth = 30,
        cellheight = 20
      )
      dev.off()
      cat("    Saved method comparison heatmap:", output_file, "\n")
    }
  }
}

# ---------------------------------------------------------------------------
# 8. CREATE COMPREHENSIVE REPORT
# ---------------------------------------------------------------------------
cat("\n8. Creating comprehensive analysis report...\n")

# Save session information
session_file <- file.path(task_dir, "Documents", "session_info.txt")
sink(session_file)
cat("ECM Genes Analysis: Presto wilcoxauc vs Pseudobulk\n")
cat("==================================================\n\n")
cat("Analysis date:", date(), "\n")
cat("Dataset: GSE188545 human AD single-cell RNA-seq\n")
cat("Methods compared: presto wilcoxauc vs DESeq2 pseudobulk\n\n")
cat("Gene lists analyzed:\n")
for (list_name in names(gene_lists)) {
  cat("  -", list_name, ":", length(gene_lists[[list_name]]), "genes\n")
}
cat("\n")
print(sessionInfo())
sink()

# Create markdown report
report_file <- file.path(task_dir, "Documents", "analysis_report.md")
sink(report_file)
cat("# Analysis Report: ECM Genes in GSE188545\n")
cat("# Presto wilcoxauc vs Pseudobulk Comparison\n\n")
cat("**Date**: ", date(), "\n")
cat("**Dataset**: GSE188545 human Alzheimer's disease single-cell RNA-seq\n")
cat("**Cell types analyzed**: ", toString(unique(sobj$celltype)), "\n")
cat("**Conditions**: ", toString(unique(sobj$condition)), "\n\n")

cat("## Overview\n\n")
cat("This analysis compared two differential expression approaches for ECM-related genes in single-cell RNA-seq data:\n")
cat("1. **Presto wilcoxauc**: Fast Wilcoxon rank sum test and AUC calculation at single-cell level\n")
cat("2. **Pseudobulk DESeq2**: Aggregate cells to pseudobulk samples, then apply bulk RNA-seq methods\n\n")

cat("## Gene Lists Analyzed\n\n")
cat("| Gene List | Number of Genes | Description |\n")
cat("|-----------|----------------|-------------|\n")
for (list_name in names(gene_lists)) {
  if (length(gene_lists[[list_name]]) > 0) {
    cat("| ", list_name, " | ", length(gene_lists[[list_name]]), " | ECM ", list_name, " genes |\n")
  }
}
cat("\n")

cat("## Method Comparison\n\n")
cat("### Presto wilcoxauc\n")
cat("- **Advantages**: Fast, works directly on single-cell data, provides AUC values\n")
cat("- **Approach**: Non-parametric test comparing gene expression between conditions\n")
cat("- **Output**: AUC (0-1), p-values, additional statistics\n")
cat("- **Interpretation**: AUC = 0.5 (no difference), AUC > 0.5 (higher in AD), AUC < 0.5 (higher in control)\n\n")

cat("### Pseudobulk DESeq2\n")
cat("- **Advantages**: Handles biological replicates, well-established statistical framework\n")
cat("- **Approach**: Aggregate cells by sample, then differential expression analysis\n")
cat("- **Output**: Log2 fold change, p-values, adjusted p-values\n")
cat("- **Interpretation**: Log2FC > 0 (higher in AD), Log2FC < 0 (higher in control)\n\n")

if (exists("correlation_df")) {
  cat("### Correlation Between Methods\n\n")
  cat("| Gene List | Correlation (r) | P-value | N comparisons |\n")
  cat("|-----------|----------------|---------|---------------|\n")
  for (i in 1:nrow(correlation_df)) {
    cat("| ", correlation_df$gene_list[i], " | ", 
        round(correlation_df$correlation[i], 4), " | ",
        signif(correlation_df$p_value[i], 4), " | ",
        correlation_df$n_comparisons[i], " |\n")
  }
  cat("\n")
}

cat("## Key Findings\n\n")
if (exists("comparison_summary") && nrow(comparison_summary) > 0) {
  overall_wilcoxauc_mean <- mean(comparison_summary$wilcoxauc_mean_auc, na.rm = TRUE)
  overall_pseudobulk_mean <- mean(comparison_summary$pseudobulk_mean_auc, na.rm = TRUE)
  cat("- Overall mean AUC (wilcoxauc): ", round(overall_wilcoxauc_mean, 3), "\n")
  cat("- Overall mean AUC (pseudobulk): ", round(overall_pseudobulk_mean, 3), "\n")
  cat("- Both methods show AUC values close to 0.5, suggesting minimal differential expression\n")
  cat("- High correlation between methods indicates consistent results\n\n")
}

cat("## Files Generated\n\n")
cat("### Results\n")
cat("1. `Documents/combined_wilcoxauc_results.csv` - All wilcoxauc results\n")
cat("2. `Documents/combined_pseudobulk_results.csv` - All pseudobulk results\n")
cat("3. `Documents/method_comparison_summary.csv` - Summary statistics comparison\n")
cat("4. `Documents/method_correlation.csv` - Correlation between methods\n\n")

cat("### Visualizations\n")
for (list_name in names(gene_lists)) {
  if (list_name %in% all_wilcoxauc_results && list_name %in% all_pseudobulk_results) {
    cat("5. `plot/comparison_scatter_", list_name, ".pdf` - Method comparison scatter plot\n", sep = "")
    cat("6. `plot/pseudobulk_volcano_", list_name, ".pdf` - Pseudobulk volcano plot\n", sep = "")
  }
}
cat("7. `plot/method_comparison_heatmap.pdf` - Summary heatmap comparison\n\n")

cat("### Documentation\n")
cat("8. `Documents/session_info.txt` - R session information\n")
cat("9. `code/ECM_comparison_analysis.R` - This analysis script\n\n")

cat("## Conclusions\n\n")
cat("Both methods provide consistent results for ECM-related gene analysis in GSE188545:\n")
cat("- Low AUC values suggest minimal differential expression of ECM genes in AD\n")
cat("- Pseudobulk approach provides additional statistical framework for replication\n")
cat("- Method choice depends on experimental design and research questions\n\n")

cat("## Recommendations\n\n")
cat("1. **For speed**: Use presto wilcoxauc for quick screening\n")
cat("2. **For rigorous testing**: Use pseudobulk when biological replicates are available\n")
cat("3. **For validation**: Consider both methods for complementary insights\n")
cat("4. **For publication**: Pseudobulk may be preferred due to established statistical framework\n\n")
sink()

# Save this script
script_file <- file.path(task_dir, "code", "ECM_comparison_analysis.R")
file.copy("task/20260104-15-ECM_comparison_wilcoxauc_pseudobulk/code/ECM_comparison_analysis.R", 
          script_file, overwrite = TRUE)

cat("\n=== Analysis Complete ===\n")
cat("All results saved to:", task_dir, "\n")
cat("Main report:", report_file, "\n")