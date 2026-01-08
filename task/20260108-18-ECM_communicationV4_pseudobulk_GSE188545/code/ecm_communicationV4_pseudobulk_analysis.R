# Pseudobulk DEG Analysis for ECM Communication V4 Genes
# Dataset: GSE188545
# Method: DESeq2 with AggregateExpression using 'sample' variable

# Load required libraries
library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)
library(pheatmap)
library(viridis)
library(ggplot2)
library(here)

# Set working directory to project root
project_root <- here::here()
setwd(project_root)

# Create output directories
task_dir <- "task/20260108-18-ECM_communicationV4_pseudobulk_GSE188545"
dir.create(file.path(task_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(task_dir, "plot"), showWarnings = FALSE, recursive = TRUE)

cat("=== ECM Communication V4 Pseudobulk DEG Analysis ===\n")
cat("Dataset: GSE188545\n\n")

# ---------------------------------------------------------------------------
# 1. LOAD GSE188545 ANNOTATED DATA
# ---------------------------------------------------------------------------
cat("1. Loading GSE188545 annotated data...\n")

sobj_path <- "data/GSE188545/GSE188545_sobj_annotated_fast.rds"
if (!file.exists(sobj_path)) {
  stop("Annotated GSE188545 object not found at: ", sobj_path)
}

sobj <- readRDS(sobj_path)
cat("   Loaded:", ncol(sobj), "cells,", nrow(sobj), "genes\n")

# Check metadata
cat("   Available metadata columns:", paste(colnames(sobj@meta.data), collapse = ", "), "\n")
cat("   Conditions:", toString(unique(sobj$condition)), "\n")
cat("   Cell types:", toString(unique(sobj$celltype)), "\n")
cat("   Samples:", toString(unique(sobj$sample)), "\n")
cat("   Number of samples:", length(unique(sobj$sample)), "\n")

# Check sample distribution
sample_table <- table(sobj$sample, sobj$condition)
cat("\n   Sample distribution:\n")
print(sample_table)

# Join layers for Seurat v5 compatibility if needed
if ("Seurat" %in% class(sobj) && packageVersion("Seurat") >= "5.0.0") {
  if (length(Layers(sobj)) > 1) {
    cat("\n   Joining layers for Seurat v5 compatibility...\n")
    sobj <- JoinLayers(sobj)
  }
}

# ---------------------------------------------------------------------------
# 2. LOAD ECM COMMUNICATION V4 GENE LIST
# ---------------------------------------------------------------------------
cat("\n2. Loading ECM communication V4 gene list...\n")

ecm_file <- "data/ECM_related_genes/ECM_communication_genesV4.xlsx"
if (!file.exists(ecm_file)) {
  stop("ECM communication V4 gene list not found at: ", ecm_file)
}

gene_data <- read.xlsx(ecm_file, sheet = 1)
gene_list <- unique(gene_data[[1]][!is.na(gene_data[[1]])])
cat("   Loaded", length(gene_list), "genes\n\n")

# ---------------------------------------------------------------------------
# 3. PSEUDOBULK ANALYSIS FUNCTION
# ---------------------------------------------------------------------------
cat("3. Defining pseudobulk analysis function...\n")

perform_pseudobulk_deseq2 <- function(seurat_obj, genes_of_interest,
                                      min_samples_per_condition = 2,
                                      min_cells_per_sample = 5) {

  cat("\n   === Pseudobulk DESeq2 Analysis ===\n")

  cell_types <- unique(seurat_obj$celltype)
  all_results <- list()

  for (ct in cell_types) {
    cat("\n   --- Cell type:", ct, "---\n")

    # Skip Unknown cell type
    if (ct == "Unknown") {
      cat("      Skipping: Unknown cell type\n")
      next
    }

    # Subset to current cell type
    ct_obj <- subset(seurat_obj, subset = celltype == ct)
    cat("      Cells in this type:", ncol(ct_obj), "\n")

    # Check sample and condition distribution
    sample_counts <- table(ct_obj$sample, ct_obj$condition)
    cat("      Sample distribution:\n")
    print(sample_counts)

    # Count samples per condition
    samples_per_condition <- table(ct_obj$condition)
    cat("      Samples per condition:", paste(names(samples_per_condition), samples_per_condition,
                                               sep = "=", collapse = ", "), "\n")

    # Check if we have enough samples per condition
    if (any(samples_per_condition < min_samples_per_condition)) {
      cat("      Skipping: insufficient samples per condition (need ≥", min_samples_per_condition, ")\n")
      next
    }

    # Check total cells per sample
    cells_per_sample <- table(ct_obj$sample)
    if (any(cells_per_sample < min_cells_per_sample)) {
      cat("      Skipping: some samples have <", min_cells_per_sample, "cells\n")
      next
    }

    tryCatch({
      # Create pseudobulk counts by aggregating cells from the same sample
      cat("      Creating pseudobulk counts...\n")

      pseudobulk_counts <- AggregateExpression(
        ct_obj,
        assays = "RNA",
        slot = "counts",
        group.by = "sample",
        return.seurat = FALSE
      )

      count_matrix <- pseudobulk_counts$RNA
      cat("      Count matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
      cat("      Sample column names:", paste(colnames(count_matrix), collapse = ", "), "\n")

      # Note: AggregateExpression replaces underscores with dashes in sample IDs
      # We need to create the correct sample metadata mapping

      # Get unique sample-condition pairs for this cell type
      orig_metadata <- ct_obj@meta.data %>%
        select(sample, condition) %>%
        distinct() %>%
        arrange(sample)

      cat("      Original sample metadata:\n")
      print(orig_metadata)

      # Create mapping: replace underscores with dashes to match column names
      sample_mapping <- data.frame(
        original_id = orig_metadata$sample,
        modified_id = gsub("_", "-", orig_metadata$sample),
        condition = orig_metadata$condition,
        stringsAsFactors = FALSE
      )

      # Ensure order matches count matrix columns
      sample_mapping <- sample_mapping %>%
        filter(modified_id %in% colnames(count_matrix)) %>%
        arrange(match(modified_id, colnames(count_matrix)))

      cat("      Sample ID mapping:\n")
      print(sample_mapping)

      # Filter to genes of interest
      available_genes <- intersect(genes_of_interest, rownames(count_matrix))
      if (length(available_genes) < 5) {
        cat("      Skipping: insufficient genes (", length(available_genes), ")\n")
        next
      }

      count_matrix <- count_matrix[available_genes, ]
      cat("      Genes after filtering:", nrow(count_matrix), "\n")

      # Create colData for DESeq2
      coldata <- data.frame(
        row.names = sample_mapping$modified_id,
        condition = sample_mapping$condition
      )

      # Ensure count matrix columns match coldata rows
      count_matrix <- count_matrix[, rownames(coldata), drop = FALSE]

      cat("      DESeq2 colData:\n")
      print(coldata)

      # Create DESeq2 dataset
      cat("      Creating DESeq2 dataset...\n")
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = coldata,
        design = ~ condition
      )

      # Pre-filter low counts
      keep <- rowSums(counts(dds) >= 10) >= 2
      dds <- dds[keep, ]
      cat("      Genes after pre-filtering:", nrow(dds), "\n")

      # Run DESeq2
      cat("      Running DESeq2...\n")
      dds <- DESeq(dds, quiet = TRUE)

      # Get results
      # Check which comparison exists
      res_names <- resultsNames(dds)
      cat("      Available results:", paste(res_names, collapse = ", "), "\n")

      # Try to get AD vs HC comparison
      if ("condition_AD_vs_HC" %in% res_names) {
        res <- results(dds, name = "condition_AD_vs_HC")
      } else {
        # Use the first available result
        res <- results(dds)
      }

      # Format results
      ct_result <- data.frame(
        gene = rownames(res),
        baseMean = res$baseMean,
        log2_fold_change = res$log2FoldChange,
        lfcSE = res$lfcSE,
        stat = res$stat,
        p_value = res$pvalue,
        padj = res$padj,
        celltype = ct,
        comparison = "AD_vs_HC",
        method = "pseudobulk_DESeq2",
        stringsAsFactors = FALSE
      )

      # Remove rows with NA in critical columns
      ct_result <- ct_result[!is.na(ct_result$p_value), ]

      if (nrow(ct_result) > 0) {
        all_results[[ct]] <- ct_result
        cat("      Results:", nrow(ct_result), "genes\n")
        cat("      Significant (padj < 0.05):", sum(ct_result$padj < 0.05, na.rm = TRUE), "genes\n")
        cat("      |log2FC| > 1:", sum(abs(ct_result$log2_fold_change) > 1, na.rm = TRUE), "genes\n")
      }

    }, error = function(e) {
      cat("      ERROR:", e$message, "\n")
    })
  }

  # Combine all results
  if (length(all_results) == 0) {
    return(data.frame())
  }

  combined <- do.call(rbind, all_results)
  cat("\n   Total results:", nrow(combined), "gene-cell type combinations\n")

  return(combined)
}

# ---------------------------------------------------------------------------
# 4. RUN PSEUDOBULK ANALYSIS
# ---------------------------------------------------------------------------
cat("\n4. Running pseudobulk DESeq2 analysis...\n")

pseudobulk_results <- perform_pseudobulk_deseq2(
  sobj,
  gene_list,
  min_samples_per_condition = 2,
  min_cells_per_sample = 5
)

if (nrow(pseudobulk_results) == 0) {
  stop("No pseudobulk results obtained. Please check the data and parameters.")
}

# Save results
output_file <- file.path(task_dir, "Documents", "ECM_communicationV4_pseudobulk_DEG_results.csv")
write.csv(pseudobulk_results, output_file, row.names = FALSE)
cat("\n   Saved results:", output_file, "\n")

# ---------------------------------------------------------------------------
# 5. CREATE SUMMARY STATISTICS
# ---------------------------------------------------------------------------
cat("\n5. Creating summary statistics...\n")

# Summary per cell type
summary_by_celltype <- pseudobulk_results %>%
  group_by(celltype) %>%
  summarise(
    genes_tested = n(),
    mean_log2FC = mean(log2_fold_change, na.rm = TRUE),
    sd_log2FC = sd(log2_fold_change, na.rm = TRUE),
    significant_padj_0.05 = sum(padj < 0.05, na.rm = TRUE),
    significant_padj_0.01 = sum(padj < 0.01, na.rm = TRUE),
    abs_log2FC_gt_1 = sum(abs(log2_fold_change) > 1, na.rm = TRUE),
    abs_log2FC_gt_0.5 = sum(abs(log2_fold_change) > 0.5, na.rm = TRUE),
    sig_and_lfc_gt_1 = sum(padj < 0.05 & abs(log2_fold_change) > 1, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_by_celltype)

summary_file <- file.path(task_dir, "Documents", "summary_statistics_per_celltype.csv")
write.csv(summary_by_celltype, summary_file, row.names = FALSE)
cat("   Saved:", summary_file, "\n")

# Overall summary
overall_summary <- pseudobulk_results %>%
  summarise(
    dataset = "GSE188545",
    cell_types_analyzed = n_distinct(celltype),
    total_tests = n(),
    genes_tested = n_distinct(gene),
    mean_log2FC = mean(log2_fold_change, na.rm = TRUE),
    sd_log2FC = sd(log2_fold_change, na.rm = TRUE),
    significant_padj_0.05 = sum(padj < 0.05, na.rm = TRUE),
    significant_padj_0.01 = sum(padj < 0.01, na.rm = TRUE),
    abs_log2FC_gt_1 = sum(abs(log2_fold_change) > 1, na.rm = TRUE),
    sig_and_lfc_gt_1 = sum(padj < 0.05 & abs(log2_fold_change) > 1, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n   Overall Summary:\n")
print(overall_summary)

overall_summary_file <- file.path(task_dir, "Documents", "overall_summary.csv")
write.csv(overall_summary, overall_summary_file, row.names = FALSE)
cat("   Saved:", overall_summary_file, "\n")

# ---------------------------------------------------------------------------
# 6. IDENTIFY SIGNIFICANT GENES
# ---------------------------------------------------------------------------
cat("\n6. Identifying significant genes...\n")

# Genes with padj < 0.05
sig_genes_0.05 <- pseudobulk_results %>%
  filter(padj < 0.05) %>%
  arrange(padj)

if (nrow(sig_genes_0.05) > 0) {
  cat("   Found", nrow(sig_genes_0.05), "significant genes (padj < 0.05)\n")

  sig_file <- file.path(task_dir, "Documents", "significant_genes_padj_0.05.csv")
  write.csv(sig_genes_0.05, sig_file, row.names = FALSE)
  cat("   Saved:", sig_file, "\n")

  # Print top 20
  cat("\n   Top 20 significant genes:\n")
  print(head(sig_genes_0.05, 20))
} else {
  cat("   No significant genes found at padj < 0.05\n")
}

# Genes with padj < 0.05 and |log2FC| > 1
sig_lfc_genes <- pseudobulk_results %>%
  filter(padj < 0.05 & abs(log2_fold_change) > 1) %>%
  arrange(padj)

if (nrow(sig_lfc_genes) > 0) {
  cat("\n   Found", nrow(sig_lfc_genes), "genes with padj < 0.05 AND |log2FC| > 1\n")

  sig_lfc_file <- file.path(task_dir, "Documents", "significant_genes_padj_0.05_abs_log2FC_gt_1.csv")
  write.csv(sig_lfc_genes, sig_lfc_file, row.names = FALSE)
  cat("   Saved:", sig_lfc_file, "\n")

  # Print top 20
  cat("\n   Top 20 significant genes with large fold changes:\n")
  print(head(sig_lfc_genes, 20))
} else {
  cat("\n   No genes found with padj < 0.05 AND |log2FC| > 1\n")
}

# ---------------------------------------------------------------------------
# 7. CREATE VISUALIZATIONS
# ---------------------------------------------------------------------------
cat("\n7. Creating visualizations...\n")

# 7a. Volcano plots for each cell type
for (ct in unique(pseudobulk_results$celltype)) {
  cat("   Creating volcano plot for", ct, "\n")

  ct_data <- pseudobulk_results %>%
    filter(celltype == ct) %>%
    mutate(
      significant = case_when(
        padj < 0.05 & abs(log2_fold_change) > 1 ~ "Both",
        padj < 0.05 ~ "Significant",
        abs(log2_fold_change) > 1 ~ "Large FC",
        TRUE ~ "NS"
      )
    )

  p <- ggplot(ct_data, aes(x = log2_fold_change, y = -log10(p_value))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(
      "Both" = "red",
      "Significant" = "orange",
      "Large FC" = "blue",
      "NS" = "gray"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    labs(
      title = paste("Volcano Plot -", ct),
      subtitle = "ECM Communication V4 Genes",
      x = "Log2 Fold Change (AD vs HC)",
      y = "-Log10 P-value",
      color = "Significance"
    ) +
    theme_minimal()

  volcano_file <- file.path(task_dir, "plot", paste0("volcano_", gsub(" ", "_", ct), ".pdf"))
  ggsave(volcano_file, plot = p, width = 8, height = 6)
  cat("      Saved:", volcano_file, "\n")
}

# 7b. Heatmap of log2 fold changes
cat("\n   Creating heatmap of log2 fold changes...\n")

heatmap_data <- pseudobulk_results %>%
  select(gene, celltype, log2_fold_change) %>%
  filter(!is.na(log2_fold_change)) %>%
  pivot_wider(
    names_from = celltype,
    values_from = log2_fold_change,
    values_fill = 0
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

if (nrow(heatmap_data) > 1 && ncol(heatmap_data) > 1) {
  heatmap_file <- file.path(task_dir, "plot", "log2_fold_change_heatmap.pdf")

  pdf(heatmap_file, width = max(10, ncol(heatmap_data) * 0.5),
      height = max(8, nrow(heatmap_data) * 0.1))

  pheatmap(
    heatmap_data,
    main = "Log2 Fold Change (AD vs HC) - ECM Communication V4 Genes",
    color = viridis(100),
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 10,
    treeheight_col = 10,
    fontsize_row = 8,
    fontsize_col = 10,
    angle_col = 45,
    cellwidth = 15,
    cellheight = 3,
    border_color = NA,
    breaks = seq(-3, 3, length.out = 100)
  )

  dev.off()
  cat("      Saved:", heatmap_file, "\n")
}

# 7c. Summary bar plot
p_summary <- ggplot(summary_by_celltype, aes(x = celltype, y = genes_tested)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = genes_tested), vjust = -0.5) +
  labs(
    title = "Genes Tested per Cell Type",
    x = "Cell Type",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary_bar_file <- file.path(task_dir, "plot", "summary_genes_tested.pdf")
ggsave(summary_bar_file, plot = p_summary, width = 10, height = 6)
cat("      Saved:", summary_bar_file, "\n")

# 7d. Significant genes bar plot
p_sig <- ggplot(summary_by_celltype, aes(x = celltype, y = sig_and_lfc_gt_1)) +
  geom_bar(stat = "identity", fill = "coral") +
  geom_text(aes(label = sig_and_lfc_gt_1), vjust = -0.5) +
  labs(
    title = "Significant Genes (padj < 0.05, |log2FC| > 1) per Cell Type",
    x = "Cell Type",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sig_bar_file <- file.path(task_dir, "plot", "summary_significant_genes.pdf")
ggsave(sig_bar_file, plot = p_sig, width = 10, height = 6)
cat("      Saved:", sig_bar_file, "\n")

# ---------------------------------------------------------------------------
# 8. CREATE ANALYSIS SUMMARY DOCUMENT
# ---------------------------------------------------------------------------
cat("\n8. Creating analysis summary document...\n")

summary_doc <- file.path(task_dir, "Documents", "analysis_summary.md")
sink(summary_doc)

cat("# Pseudobulk DEG Analysis Summary: ECM Communication V4 Genes\n\n")
cat("**Task Directory**: `20260108-18-ECM_communicationV4_pseudobulk_GSE188545`\n")
cat("**Date**: ", date(), "\n")
cat("**Dataset**: GSE188545 human Alzheimer's disease single-cell RNA-seq from middle temporal gyrus (MTG)\n\n")

cat("## Overview\n\n")
cat("This task performed pseudobulk differential expression analysis using DESeq2 to identify ECM communication V4 genes that show differential expression between Alzheimer's Disease (AD) and Healthy Control (HC) samples across different brain cell types. The analysis aggregates single-cell counts by sample before performing DESeq2, which accounts for sample-level variability and provides more robust results.\n\n")

cat("## Methodology\n\n")
cat("### Gene List Analyzed\n")
cat("- **ECM communication genes V4** (`ECM_communication_genesV4.xlsx`)\n")
cat("- ", length(gene_list), " genes in list\n\n")

cat("### Pseudobulk Method\n")
cat("- **Aggregation**: Used `AggregateExpression()` with `group.by = \"sample\"`\n")
cat("- **DE analysis**: DESeq2 with design `~ condition`\n")
cat("- **Filtering**: \n")
cat("  - Minimum 2 samples per condition per cell type\n")
cat("  - Minimum 5 cells per sample\n")
cat("  - Pre-filtering: genes with ≥10 counts in ≥2 samples\n\n")

cat("### Sample Metadata Handling\n")
cat("- `AggregateExpression()` replaces underscores with dashes in sample IDs\n")
cat("- Created mapping between original sample IDs and modified column names\n")
cat("- Ensured proper matching of count matrix columns with condition metadata\n\n")

cat("## Key Results\n\n")
cat("- **Cell types analyzed**: ", overall_summary$cell_types_analyzed, "\n")
cat("- **Total tests performed**: ", overall_summary$total_tests, "\n")
cat("- **Unique genes tested**: ", overall_summary$genes_tested, "\n")
cat("- **Significant genes (padj < 0.05)**: ", overall_summary$significant_padj_0.05, "\n")
cat("- **Significant genes (padj < 0.01)**: ", overall_summary$significant_padj_0.01, "\n")
cat("- **Genes with |log2FC| > 1**: ", overall_summary$abs_log2FC_gt_1, "\n")
cat("- **Significant AND |log2FC| > 1**: ", overall_summary$sig_and_lfc_gt_1, "\n\n")

cat("### Per Cell Type Statistics\n\n")
for (i in 1:nrow(summary_by_celltype)) {
  cat("#### ", summary_by_celltype$celltype[i], "\n", sep = "")
  cat("- Genes tested: ", summary_by_celltype$genes_tested[i], "\n")
  cat("- Mean log2FC: ", round(summary_by_celltype$mean_log2FC[i], 3), "\n")
  cat("- Significant (padj < 0.05): ", summary_by_celltype$significant_padj_0.05[i], "\n")
  cat("- Significant AND |log2FC| > 1: ", summary_by_celltype$sig_and_lfc_gt_1[i], "\n\n")
}

cat("## Files Generated\n\n")
cat("### Results Tables\n")
cat("1. `Documents/ECM_communicationV4_pseudobulk_DEG_results.csv` - Complete DEG results\n")
cat("2. `Documents/summary_statistics_per_celltype.csv` - Per-celltype summary\n")
cat("3. `Documents/overall_summary.csv` - Overall summary statistics\n")
cat("4. `Documents/significant_genes_padj_0.05.csv` - Significant genes (padj < 0.05)\n")
cat("5. `Documents/significant_genes_padj_0.05_abs_log2FC_gt_1.csv` - Significant genes with large fold change\n\n")

cat("### Visualizations\n")
cat("6. `plot/volcano_<celltype>.pdf` - Volcano plots for each cell type\n")
cat("7. `plot/log2_fold_change_heatmap.pdf` - Heatmap of log2 fold changes\n")
cat("8. `plot/summary_genes_tested.pdf` - Bar plot of genes tested per cell type\n")
cat("9. `plot/summary_significant_genes.pdf` - Bar plot of significant genes per cell type\n\n")

cat("## Interpretation\n\n")
cat("### Pseudobulk Advantages\n")
cat("- Accounts for sample-to-sample variability\n")
cat("- Reduces false positives from single-cell level analysis\n")
cat("- More robust for detecting consistent differences across samples\n\n")

cat("### Log2 Fold Change Interpretation\n")
cat("- **log2FC > 0**: Higher expression in AD\n")
cat("- **log2FC < 0**: Higher expression in HC\n")
cat("- **|log2FC| > 1**: 2-fold difference or more\n")
cat("- **|log2FC| > 0.5**: ~1.4-fold difference or more\n\n")

cat("### P-value Interpretation\n")
cat("- **padj < 0.05**: Statistically significant after multiple testing correction\n")
(cat("- **padj < 0.01**: Highly statistically significant\n"))

sink()

cat("   Saved:", summary_doc, "\n")

# ---------------------------------------------------------------------------
# 9. SESSION INFO
# ---------------------------------------------------------------------------
session_file <- file.path(task_dir, "Documents", "session_info.txt")
sink(session_file)
cat("ECM Communication V4 Pseudobulk DEG Analysis\n")
cat("============================================\n\n")
cat("Dataset: GSE188545\n")
cat("Analysis date:", date(), "\n\n")
print(sessionInfo())
sink()

cat("\n=== Analysis Complete ===\n")
cat("All results saved to:", task_dir, "\n")
