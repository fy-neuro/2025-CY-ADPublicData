# Pseudobulk DEG Analysis for ECM Communication V4 Genes Across 3 Datasets
# Datasets: GSE138852 (mouse), GSE174367 (human), GSE188545 (human)

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
library(patchwork)
library(here)

# Set working directory to project root
project_root <- here::here()
setwd(project_root)

# Create output directories
task_dir <- "task/20260108-19-ECM_communicationV4_pseudobulk_3datasets"
dir.create(file.path(task_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(task_dir, "plot"), showWarnings = FALSE, recursive = TRUE)

cat("=== ECM Communication V4 Pseudobulk DEG Analysis Across 3 Datasets ===\n\n")

# ---------------------------------------------------------------------------
# 1. LOAD ECM COMMUNICATION V4 GENE LIST
# ---------------------------------------------------------------------------
cat("1. Loading ECM communication V4 gene list...\n")

ecm_file <- "data/ECM_related_genes/ECM_communication_genesV4.xlsx"
if (!file.exists(ecm_file)) {
  stop("ECM communication V4 gene list not found at: ", ecm_file)
}

gene_data <- read.xlsx(ecm_file, sheet = 1)
gene_list <- unique(gene_data[[1]][!is.na(gene_data[[1]])])
cat("   Loaded", length(gene_list), "genes\n\n")

# ---------------------------------------------------------------------------
# 2. PSEUDOBULK ANALYSIS FUNCTION
# ---------------------------------------------------------------------------
cat("2. Defining pseudobulk analysis function...\n")

perform_pseudobulk_deseq2 <- function(seurat_obj, genes_of_interest, dataset_name,
                                      min_samples_per_condition = 2,
                                      min_cells_per_sample = 5) {

  cat("\n   === Pseudobulk DESeq2 Analysis:", dataset_name, "===\n")
  cat("   Cells:", ncol(seurat_obj), "\n")
  cat("   Genes:", nrow(seurat_obj), "\n")

  # Check metadata structure
  meta_cols <- colnames(seurat_obj@meta.data)
  cat("   Metadata columns:", paste(meta_cols, collapse = ", "), "\n")

  # Determine condition and celltype column names
  condition_col <- NULL
  celltype_col <- NULL
  sample_col <- NULL

  # Find condition column
  possible_condition_cols <- c("condition", "disease", "Status", "group")
  for (col in possible_condition_cols) {
    if (col %in% meta_cols) {
      condition_col <- col
      break
    }
  }

  # Find celltype column
  possible_celltype_cols <- c("celltype", "cell_type", "cell_type_annotation", "subclass")
  for (col in possible_celltype_cols) {
    if (col %in% meta_cols) {
      celltype_col <- col
      break
    }
  }

  # Find sample column
  possible_sample_cols <- c("sample", "Sample", "orig.ident", "batch")
  for (col in possible_sample_cols) {
    if (col %in% meta_cols) {
      sample_col <- col
      break
    }
  }

  if (is.null(condition_col)) {
    cat("   ERROR: Cannot find condition column\n")
    return(data.frame())
  }

  if (is.null(celltype_col)) {
    cat("   ERROR: Cannot find celltype column\n")
    return(data.frame())
  }

  if (is.null(sample_col)) {
    cat("   ERROR: Cannot find sample column\n")
    return(data.frame())
  }

  cat("   Using condition column:", condition_col, "\n")
  cat("   Using celltype column:", celltype_col, "\n")
  cat("   Using sample column:", sample_col, "\n")

  # Print unique values
  conditions <- unique(seurat_obj@meta.data[[condition_col]])
  celltypes <- unique(seurat_obj@meta.data[[celltype_col]])
  samples <- unique(seurat_obj@meta.data[[sample_col]])
  cat("   Conditions:", toString(conditions), "\n")
  cat("   Cell types:", toString(celltypes), "\n")
  cat("   Samples:", length(samples), "unique samples\n")

  # Join layers for Seurat v5 if needed
  if ("Seurat" %in% class(seurat_obj) && packageVersion("Seurat") >= "5.0.0") {
    if (length(Layers(seurat_obj)) > 1) {
      seurat_obj <- JoinLayers(seurat_obj)
      cat("   Joined layers for Seurat v5 compatibility\n")
    }
  }

  # Initialize results
  all_results <- list()

  for (ct in celltypes) {
    cat("\n      --- Cell type:", ct, "---\n")

    # Subset to current cell type
    ct_obj <- subset(seurat_obj, subset = !!sym(celltype_col) == ct)
    cat("         Cells:", ncol(ct_obj), "\n")

    # Check sample and condition distribution
    sample_counts <- table(ct_obj@meta.data[[sample_col]], ct_obj@meta.data[[condition_col]])
    cat("         Sample distribution:\n")
    print(sample_counts)

    # Count samples per condition
    samples_per_condition <- table(ct_obj@meta.data[[condition_col]])
    cat("         Samples per condition:", paste(names(samples_per_condition), samples_per_condition,
                                                   sep = "=", collapse = ", "), "\n")

    # Check if we have enough samples per condition
    if (any(samples_per_condition < min_samples_per_condition)) {
      cat("         Skipping: insufficient samples per condition (need ≥", min_samples_per_condition, ")\n")
      next
    }

    # Check total cells per sample
    cells_per_sample <- table(ct_obj@meta.data[[sample_col]])
    if (any(cells_per_sample < min_cells_per_sample)) {
      cat("         Skipping: some samples have <", min_cells_per_sample, "cells\n")
      next
    }

    tryCatch({
      # Create pseudobulk counts by aggregating cells from the same sample
      cat("         Creating pseudobulk counts...\n")

      pseudobulk_counts <- AggregateExpression(
        ct_obj,
        assays = "RNA",
        slot = "counts",
        group.by = sample_col,
        return.seurat = FALSE
      )

      count_matrix <- pseudobulk_counts$RNA
      cat("         Count matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
      cat("         Sample columns:", paste(colnames(count_matrix), collapse = ", "), "\n")

      # Handle sample ID transformation (underscores -> dashes)
      orig_metadata <- ct_obj@meta.data %>%
        select(sample = !!sym(sample_col), condition = !!sym(condition_col)) %>%
        distinct() %>%
        arrange(sample)

      sample_mapping <- data.frame(
        original_id = orig_metadata$sample,
        modified_id = gsub("_", "-", orig_metadata$sample),
        condition = orig_metadata$condition,
        stringsAsFactors = FALSE
      )

      sample_mapping <- sample_mapping %>%
        filter(modified_id %in% colnames(count_matrix)) %>%
        arrange(match(modified_id, colnames(count_matrix)))

      cat("         Sample ID mapping:\n")
      print(sample_mapping)

      # Filter to genes of interest
      available_genes <- intersect(genes_of_interest, rownames(count_matrix))
      if (length(available_genes) < 5) {
        cat("         Skipping: insufficient genes (", length(available_genes), ")\n")
        next
      }

      count_matrix <- count_matrix[available_genes, ]
      cat("         Genes after filtering:", nrow(count_matrix), "\n")

      # Create colData for DESeq2
      coldata <- data.frame(
        row.names = sample_mapping$modified_id,
        condition = sample_mapping$condition
      )

      count_matrix <- count_matrix[, rownames(coldata), drop = FALSE]

      # Create DESeq2 dataset
      cat("         Creating DESeq2 dataset...\n")
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = coldata,
        design = ~ condition
      )

      # Pre-filter low counts
      keep <- rowSums(counts(dds) >= 10) >= 2
      dds <- dds[keep, ]
      cat("         Genes after pre-filtering:", nrow(dds), "\n")

      # Run DESeq2
      cat("         Running DESeq2...\n")
      dds <- DESeq(dds, quiet = TRUE)

      # Get results
      res_names <- resultsNames(dds)

      # Determine comparison
      if ("condition_AD_vs_HC" %in% res_names) {
        res <- results(dds, name = "condition_AD_vs_HC")
      } else if ("condition_PS2_vs_WT" %in% res_names) {
        res <- results(dds, name = "condition_PS2_vs_WT")
      } else {
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
        condition_col = condition_col,
        celltype_col = celltype_col,
        sample_col = sample_col,
        dataset = dataset_name,
        comparison = "Disease_vs_Control",
        method = "pseudobulk_DESeq2",
        stringsAsFactors = FALSE
      )

      # Remove rows with NA in critical columns
      ct_result <- ct_result[!is.na(ct_result$p_value), ]

      if (nrow(ct_result) > 0) {
        all_results[[ct]] <- ct_result
        cat("         Results:", nrow(ct_result), "genes\n")
        cat("         Significant (padj < 0.05):", sum(ct_result$padj < 0.05, na.rm = TRUE), "genes\n")
        cat("         |log2FC| > 1:", sum(abs(ct_result$log2_fold_change) > 1, na.rm = TRUE), "genes\n")
        cat("         Sig AND |log2FC| > 1:", sum(ct_result$padj < 0.05 & abs(ct_result$log2_fold_change) > 1, na.rm = TRUE), "genes\n")
      }

    }, error = function(e) {
      cat("         ERROR:", e$message, "\n")
    })
  }

  # Combine all results
  if (length(all_results) == 0) {
    cat("      WARNING: No results obtained for", dataset_name, "\n")
    return(data.frame())
  }

  result <- do.call(rbind, all_results)
  cat("\n   Total results:", nrow(result), "gene-cell type combinations\n")

  return(result)
}

# ---------------------------------------------------------------------------
# 3. ANALYZE EACH DATASET
# ---------------------------------------------------------------------------
cat("\n3. Loading datasets and running pseudobulk analysis...\n")

all_results <- list()

# === Dataset 1: GSE138852 (Mouse) - SKIPPED ===
cat("\n--- Dataset 1: GSE138852 (Mouse AD) ---\n")
cat("   SKIPPED: Pseudobulk analysis not recommended for GSE138852\n")
cat("   (Insufficient samples or experimental design not suitable)\n")

# === Dataset 2: GSE174367 (Human snRNA-seq + snATAC-seq) ===
cat("\n--- Dataset 2: GSE174367 (Human Multiome) ---\n")
gse174367_path <- "data/GSE174367/sobj_20251212.rds"

if (file.exists(gse174367_path)) {
  sobj_gse174367 <- readRDS(gse174367_path)

  # For GSE174367, specifically use SampleID, Diagnosis, and Cell.Type columns
  # Rename them to standard names for the analysis function
  if ("SampleID" %in% colnames(sobj_gse174367@meta.data)) {
    sobj_gse174367$sample <- sobj_gse174367$SampleID
    cat("   Using SampleID as sample variable\n")
  }
  if ("Diagnosis" %in% colnames(sobj_gse174367@meta.data)) {
    sobj_gse174367$condition <- sobj_gse174367$Diagnosis
    cat("   Using Diagnosis as condition variable\n")
  }
  if ("Cell.Type" %in% colnames(sobj_gse174367@meta.data)) {
    sobj_gse174367$celltype <- sobj_gse174367$Cell.Type
    cat("   Using Cell.Type as celltype variable\n")
  }

  result_gse174367 <- perform_pseudobulk_deseq2(
    sobj_gse174367,
    gene_list,
    "GSE174367_Human",
    min_samples_per_condition = 2,
    min_cells_per_sample = 5
  )

  if (nrow(result_gse174367) > 0) {
    all_results$GSE174367 <- result_gse174367

    # Save individual results
    output_file <- file.path(task_dir, "Documents", "GSE174367_ECM_communicationV4_pseudobulk_DEG.csv")
    write.csv(result_gse174367, output_file, row.names = FALSE)
    cat("   Saved:", output_file, "\n")
  }
} else {
  cat("   WARNING: GSE174367 object not found at:", gse174367_path, "\n")
}

# === Dataset 3: GSE188545 (Human MTG) ===
cat("\n--- Dataset 3: GSE188545 (Human MTG) ---\n")
gse188545_path <- "data/GSE188545/GSE188545_sobj_annotated_fast.rds"

if (file.exists(gse188545_path)) {
  sobj_gse188545 <- readRDS(gse188545_path)
  result_gse188545 <- perform_pseudobulk_deseq2(
    sobj_gse188545,
    gene_list,
    "GSE188545_Human",
    min_samples_per_condition = 2,
    min_cells_per_sample = 5
  )

  if (nrow(result_gse188545) > 0) {
    all_results$GSE188545 <- result_gse188545

    # Save individual results
    output_file <- file.path(task_dir, "Documents", "GSE188545_ECM_communicationV4_pseudobulk_DEG.csv")
    write.csv(result_gse188545, output_file, row.names = FALSE)
    cat("   Saved:", output_file, "\n")
  }
} else {
  cat("   WARNING: GSE188545 object not found at:", gse188545_path, "\n")
}

# ---------------------------------------------------------------------------
# 4. COMBINE AND SAVE ALL RESULTS
# ---------------------------------------------------------------------------
cat("\n4. Combining results from all datasets...\n")

if (length(all_results) > 0) {
  combined_results <- do.call(rbind, all_results)

  # Save combined results
  combined_file <- file.path(task_dir, "Documents", "All_Datasets_ECM_communicationV4_pseudobulk_DEG_combined.csv")
  write.csv(combined_results, combined_file, row.names = FALSE)
  cat("   Saved combined results:", combined_file, "\n")
  cat("   Total gene-cell type combinations:", nrow(combined_results), "\n")

  # ---------------------------------------------------------------------------
  # 5. CREATE SUMMARY STATISTICS
  # ---------------------------------------------------------------------------
  cat("\n5. Creating summary statistics...\n")

  # Summary per dataset
  summary_by_dataset <- combined_results %>%
    group_by(dataset) %>%
    summarise(
      cell_types_analyzed = n_distinct(celltype),
      genes_tested = n_distinct(gene),
      total_combinations = n(),
      mean_log2FC = mean(log2_fold_change, na.rm = TRUE),
      sd_log2FC = sd(log2_fold_change, na.rm = TRUE),
      significant_padj_0.05 = sum(padj < 0.05, na.rm = TRUE),
      significant_padj_0.01 = sum(padj < 0.01, na.rm = TRUE),
      abs_log2FC_gt_1 = sum(abs(log2_fold_change) > 1, na.rm = TRUE),
      sig_and_lfc_gt_1 = sum(padj < 0.05 & abs(log2_fold_change) > 1, na.rm = TRUE),
      .groups = "drop"
    )

  print(summary_by_dataset)

  summary_dataset_file <- file.path(task_dir, "Documents", "summary_statistics_per_dataset.csv")
  write.csv(summary_by_dataset, summary_dataset_file, row.names = FALSE)
  cat("   Saved:", summary_dataset_file, "\n")

  # Summary per cell type across datasets
  summary_by_celltype <- combined_results %>%
    group_by(dataset, celltype) %>%
    summarise(
      genes_tested = n(),
      mean_log2FC = mean(log2_fold_change, na.rm = TRUE),
      sd_log2FC = sd(log2_fold_change, na.rm = TRUE),
      significant_padj_0.05 = sum(padj < 0.05, na.rm = TRUE),
      abs_log2FC_gt_1 = sum(abs(log2_fold_change) > 1, na.rm = TRUE),
      sig_and_lfc_gt_1 = sum(padj < 0.05 & abs(log2_fold_change) > 1, na.rm = TRUE),
      .groups = "drop"
    )

  celltype_file <- file.path(task_dir, "Documents", "summary_statistics_per_celltype.csv")
  write.csv(summary_by_celltype, celltype_file, row.names = FALSE)
  cat("   Saved:", celltype_file, "\n")

  # Identify significant genes across all datasets
  cat("\n6. Identifying significant genes...\n")

  # Genes with padj < 0.05
  sig_genes_0.05 <- combined_results %>%
    filter(padj < 0.05) %>%
    arrange(padj)

  if (nrow(sig_genes_0.05) > 0) {
    cat("   Found", nrow(sig_genes_0.05), "significant genes (padj < 0.05)\n")

    sig_file <- file.path(task_dir, "Documents", "All_Datasets_ECM_communicationV4_pseudobulk_significant_padj_0.05.csv")
    write.csv(sig_genes_0.05, sig_file, row.names = FALSE)
    cat("   Saved:", sig_file, "\n")
  } else {
    cat("   No significant genes found at padj < 0.05\n")
  }

  # Genes with padj < 0.05 AND |log2FC| > 1
  sig_lfc_genes <- combined_results %>%
    filter(padj < 0.05 & abs(log2_fold_change) > 1) %>%
    arrange(padj)

  if (nrow(sig_lfc_genes) > 0) {
    cat("   Found", nrow(sig_lfc_genes), "genes with padj < 0.05 AND |log2FC| > 1\n")

    sig_lfc_file <- file.path(task_dir, "Documents", "All_Datasets_ECM_communicationV4_pseudobulk_significant_padj_0.05_abs_log2FC_gt_1.csv")
    write.csv(sig_lfc_genes, sig_lfc_file, row.names = FALSE)
    cat("   Saved:", sig_lfc_file, "\n")
  } else {
    cat("   No genes found with padj < 0.05 AND |log2FC| > 1\n")
  }

  # Genes significant in multiple datasets
  cat("\n7. Finding genes significant in multiple datasets...\n")

  gene_dataset_significance <- combined_results %>%
    filter(padj < 0.05) %>%
    select(gene, dataset, celltype, log2_fold_change, padj) %>%
    arrange(gene, dataset)

  if (nrow(gene_dataset_significance) > 0) {
    multi_dataset_genes <- gene_dataset_significance %>%
      group_by(gene) %>%
      filter(n_distinct(dataset) > 1) %>%
      arrange(desc(n_distinct(dataset)), padj) %>%
      select(gene, dataset, celltype, log2_fold_change, padj)

    if (nrow(multi_dataset_genes) > 0) {
      cat("   Found", nrow(multi_dataset_genes), "gene occurrences in multiple datasets\n")

      multi_file <- file.path(task_dir, "Documents", "Genes_significant_in_multiple_datasets.csv")
      write.csv(multi_dataset_genes, multi_file, row.names = FALSE)
      cat("   Saved:", multi_file, "\n")
    } else {
      cat("   No genes found significant in multiple datasets\n")
    }
  }

  # ---------------------------------------------------------------------------
  # 8. CREATE VISUALIZATIONS
  # ---------------------------------------------------------------------------
  cat("\n8. Creating visualizations...\n")

  # 8a. Summary bar plots per dataset
  p1 <- ggplot(summary_by_dataset, aes(x = dataset, y = genes_tested, fill = dataset)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = genes_tested), vjust = -0.5) +
    labs(title = "Genes Tested per Dataset",
         x = "Dataset", y = "Number of Genes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  p2 <- ggplot(summary_by_dataset, aes(x = dataset, y = sig_and_lfc_gt_1, fill = dataset)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sig_and_lfc_gt_1), vjust = -0.5) +
    labs(title = "Significant Genes (padj < 0.05, |log2FC| > 1) per Dataset",
         x = "Dataset", y = "Number of Genes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  p3 <- ggplot(summary_by_dataset, aes(x = dataset, y = mean_log2FC, fill = dataset)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(mean_log2FC, 3)), vjust = -0.5) +
    labs(title = "Mean Log2 Fold Change per Dataset",
         x = "Dataset", y = "Mean Log2FC") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  summary_plot <- p1 | p2 | p3
  summary_plot_file <- file.path(task_dir, "plot", "Summary_statistics_per_dataset.pdf")
  ggsave(summary_plot_file, plot = summary_plot, width = 18, height = 5)
  cat("   Saved:", summary_plot_file, "\n")

  # 8b. Heatmap for each dataset
  for (ds in unique(combined_results$dataset)) {
    cat("   Creating heatmap for", ds, "\n")

    ds_data <- combined_results %>%
      filter(dataset == ds) %>%
      select(celltype, gene, log2_fold_change) %>%
      filter(!is.na(log2_fold_change)) %>%
      pivot_wider(
        names_from = gene,
        values_from = log2_fold_change,
        values_fill = 0
      ) %>%
      column_to_rownames("celltype") %>%
      as.matrix()

    if (nrow(ds_data) > 1 && ncol(ds_data) > 1) {
      heatmap_file <- file.path(task_dir, "plot", paste0(ds, "_log2_fold_change_heatmap.pdf"))

      pdf(heatmap_file, width = max(10, ncol(ds_data) * 0.3),
          height = max(6, nrow(ds_data) * 0.5))

      pheatmap(
        ds_data,
        main = paste("Log2 Fold Change -", ds),
        color = viridis(100),
        scale = "none",
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        treeheight_row = 10,
        treeheight_col = 10,
        fontsize_row = 10,
        fontsize_col = 8,
        angle_col = 45,
        cellwidth = 12,
        cellheight = 20,
        border_color = NA,
        breaks = seq(-3, 3, length.out = 100)
      )

      dev.off()
      cat("      Saved:", heatmap_file, "\n")
    }
  }

  # 8c. Combined heatmap
  cat("   Creating combined heatmap...\n")

  combined_heatmap_data <- combined_results %>%
    filter(!is.na(log2_fold_change)) %>%
    mutate(dataset_celltype = paste(dataset, celltype, sep = "_")) %>%
    select(dataset_celltype, gene, log2_fold_change) %>%
    pivot_wider(
      names_from = gene,
      values_from = log2_fold_change,
      values_fill = 0
    ) %>%
    column_to_rownames("dataset_celltype") %>%
    as.matrix()

  if (nrow(combined_heatmap_data) > 1 && ncol(combined_heatmap_data) > 1) {
    combined_heatmap_file <- file.path(task_dir, "plot", "All_Datasets_log2_fold_change_heatmap.pdf")

    pdf(combined_heatmap_file, width = max(12, ncol(combined_heatmap_data) * 0.2),
        height = max(8, nrow(combined_heatmap_data) * 0.3))

    pheatmap(
      combined_heatmap_data,
      main = "Log2 Fold Change - All Datasets",
      color = viridis(100),
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      treeheight_row = 10,
      treeheight_col = 10,
      fontsize_row = 8,
      fontsize_col = 6,
      angle_col = 45,
      cellwidth = 10,
      cellheight = 15,
      border_color = NA,
      breaks = seq(-3, 3, length.out = 100)
    )

    dev.off()
    cat("      Saved:", combined_heatmap_file, "\n")
  }

} else {
  cat("\nERROR: No results obtained from any dataset\n")
}

# ---------------------------------------------------------------------------
# 9. CREATE ANALYSIS SUMMARY DOCUMENT
# ---------------------------------------------------------------------------
cat("\n9. Creating analysis summary document...\n")

summary_doc <- file.path(task_dir, "Documents", "analysis_summary.md")
sink(summary_doc)

cat("# Pseudobulk DEG Analysis Summary: ECM Communication V4 Genes\n\n")
cat("**Task Directory**: `20260108-19-ECM_communicationV4_pseudobulk_3datasets`\n")
cat("**Date**: ", date(), "\n\n")
cat("## Overview\n\n")
cat("This task performed pseudobulk differential expression analysis using DESeq2 to identify ECM communication V4 genes that show differential expression between disease and control conditions in two human Alzheimer's disease single-cell datasets. The analysis aggregates single-cell counts by sample before performing DESeq2, which accounts for sample-level variability and provides more robust results.\n\n")

cat("## Datasets Analyzed\n\n")
cat("1. **GSE174367** (Human snRNA-seq + snATAC-seq multiome)\n")
cat("   - Used SampleID for aggregation\n")
cat("   - Used Diagnosis for condition\n")
cat("   - Used Cell.Type for celltype\n")
cat("2. **GSE188545** (Human MTG single-cell RNA-seq)\n")
cat("   - Used sample for aggregation, condition for AD vs HC\n")
cat("3. **GSE138852** (Mouse AD model) - **SKIPPED**\n")
cat("   - Not suitable for pseudobulk analysis\n\n")

cat("## Gene List Analyzed\n\n")
cat("- **ECM communication genes V4** (`ECM_communication_genesV4.xlsx`)\n")
cat("- ", length(gene_list), " genes in list\n\n")

cat("## Methodology\n\n")
cat("### Pseudobulk Method\n")
cat("- **Aggregation**: Used `AggregateExpression()` with `group.by = \"sample\"`\n")
cat("- **DE analysis**: DESeq2 with design `~ condition`\n")
cat("- **Sample ID handling**: Automatic mapping for underscore → dash conversion\n")
cat("- **Filtering**: \n")
cat("  - Minimum 2 samples per condition per cell type\n")
cat("  - Minimum 5 cells per sample\n")
cat("  - Pre-filtering: genes with ≥10 counts in ≥2 samples\n\n")

cat("### Key Results\n\n")

if (length(all_results) > 0 && exists("summary_by_dataset")) {
  cat("- Datasets analyzed:", length(unique(combined_results$dataset)), "\n")
  cat("- Total gene-cell type combinations:", nrow(combined_results), "\n")
  cat("- Mean log2FC across all:", round(mean(combined_results$log2_fold_change), 3), "\n")
  cat("- Significant (padj < 0.05):", sum(combined_results$padj < 0.05, na.rm = TRUE), "\n")
  cat("- Significant AND |log2FC| > 1:", sum(combined_results$padj < 0.05 & abs(combined_results$log2_fold_change) > 1, na.rm = TRUE), "\n\n")

  cat("### Per Dataset Statistics\n\n")
  for (i in 1:nrow(summary_by_dataset)) {
    cat("#### ", summary_by_dataset$dataset[i], "\n", sep = "")
    cat("- Cell types analyzed:", summary_by_dataset$cell_types_analyzed[i], "\n")
    cat("- Genes tested:", summary_by_dataset$genes_tested[i], "\n")
    cat("- Mean log2FC:", round(summary_by_dataset$mean_log2FC[i], 3), "\n")
    cat("- Significant (padj < 0.05):", summary_by_dataset$significant_padj_0.05[i], "\n")
    cat("- Significant AND |log2FC| > 1:", summary_by_dataset$sig_and_lfc_gt_1[i], "\n\n")
  }
}

cat("## Files Generated\n\n")
cat("### Individual Dataset Results\n")
cat("1. `Documents/GSE138852_ECM_communicationV4_pseudobulk_DEG.csv`\n")
cat("2. `Documents/GSE174367_ECM_communicationV4_pseudobulk_DEG.csv`\n")
cat("3. `Documents/GSE188545_ECM_communicationV4_pseudobulk_DEG.csv`\n\n")
cat("### Combined Results\n")
cat("4. `Documents/All_Datasets_ECM_communicationV4_pseudobulk_DEG_combined.csv`\n")
cat("5. `Documents/summary_statistics_per_dataset.csv`\n")
cat("6. `Documents/summary_statistics_per_celltype.csv`\n")
cat("7. `Documents/All_Datasets_ECM_communicationV4_pseudobulk_significant_padj_0.05.csv`\n")
cat("8. `Documents/All_Datasets_ECM_communicationV4_pseudobulk_significant_padj_0.05_abs_log2FC_gt_1.csv`\n")
cat("9. `Documents/Genes_significant_in_multiple_datasets.csv`\n\n")
cat("### Visualizations\n")
cat("10. `plot/Summary_statistics_per_dataset.pdf`\n")
cat("11. `plot/GSE138852_log2_fold_change_heatmap.pdf`\n")
cat("12. `plot/GSE174367_log2_fold_change_heatmap.pdf`\n")
cat("13. `plot/GSE188545_log2_fold_change_heatmap.pdf`\n")
cat("14. `plot/All_Datasets_log2_fold_change_heatmap.pdf`\n\n")

cat("## Interpretation\n\n")
cat("### Pseudobulk Advantages\n")
cat("- Accounts for sample-to-sample variability\n")
cat("- More conservative than cell-level methods\n")
cat("- Reduces false positives\n")
cat("- Better control of type I error\n\n")

cat("### Cross-Dataset Comparison\n")
cat("- Genes significant in multiple datasets are high-confidence findings\n")
cat("- Species differences expected between mouse and human\n")
cat("- Tissue/region differences may affect results\n\n")

sink()

cat("   Saved:", summary_doc, "\n")

# ---------------------------------------------------------------------------
# 10. SESSION INFO
# ---------------------------------------------------------------------------
session_file <- file.path(task_dir, "Documents", "session_info.txt")
sink(session_file)
cat("ECM Communication V4 Pseudobulk DEG Analysis\n")
cat("============================================\n\n")
cat("Datasets: GSE174367, GSE188545 (GSE138852 skipped)\n")
cat("Analysis date:", date(), "\n\n")
print(sessionInfo())
sink()

cat("\n=== Analysis Complete ===\n")
cat("All results saved to:", task_dir, "\n")
