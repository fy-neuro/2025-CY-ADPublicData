# ECM Communication V4 Genes Wilcoxauc Analysis Across 3 Datasets
# Datasets: GSE138852 (mouse), GSE174367 (human), GSE188545 (human)

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
library(patchwork)
library(here)

# Set working directory to project root
project_root <- here::here()
setwd(project_root)

# Create output directories
task_dir <- "task/20260108-17-ECM_communicationV4_wilcoxauc_3datasets"
dir.create(file.path(task_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(task_dir, "plot"), showWarnings = FALSE, recursive = TRUE)

cat("=== ECM Communication V4 Wilcoxauc Analysis Across 3 Datasets ===\n\n")

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
# 2. WILCOXAUC ANALYSIS FUNCTION
# ---------------------------------------------------------------------------
cat("2. Defining wilcoxauc analysis function...\n")

perform_wilcoxauc_analysis <- function(sobj, gene_list, dataset_name, min_cells_per_group = 10) {

  cat("\n   === Analyzing", dataset_name, "===\n")
  cat("   Cells:", ncol(sobj), "\n")
  cat("   Genes:", nrow(sobj), "\n")

  # Check metadata structure
  meta_cols <- colnames(sobj@meta.data)
  cat("   Metadata columns:", paste(meta_cols, collapse = ", "), "\n")

  # Join layers for Seurat v5 if needed (only if Layers() function works)
  if ("Seurat" %in% class(sobj)) {
    tryCatch({
      if (length(Layers(sobj)) > 1) {
        sobj <- JoinLayers(sobj)
        cat("   Joined layers for Seurat v5 compatibility\n")
      }
    }, error = function(e) {
      cat("   Note: Could not join layers (not required for wilcoxauc)\n")
    })
  }

  # Get condition and celltype values (already set by dataset-specific code)
  conditions <- unique(sobj$condition)
  celltypes <- unique(sobj$celltype)
  cat("   Conditions:", toString(conditions), "\n")
  cat("   Cell types:", toString(celltypes), "\n")

  # Initialize results
  wilcoxauc_results <- list()

  for (ct in celltypes) {
    cat("      Processing cell type:", ct, "\n")

    # Subset to current cell type
    sobj_ct <- subset(sobj, subset = celltype == ct)

    # Check if both conditions present with sufficient cells
    cell_counts <- table(sobj_ct$condition)
    cat("         Cell counts:", paste(names(cell_counts), cell_counts, sep = "=", collapse = ", "), "\n")

    if (length(cell_counts) < 2) {
      cat("         Skipping: only one condition present\n")
      next
    }

    if (any(cell_counts < min_cells_per_group)) {
      cat("         Skipping: insufficient cells in one condition\n")
      next
    }

    # Find genes present in this subset
    available_genes <- rownames(sobj_ct)
    genes_to_test <- intersect(gene_list, available_genes)

    if (length(genes_to_test) == 0) {
      cat("         Skipping: no genes from list present\n")
      next
    }

    cat("         Testing", length(genes_to_test), "genes\n")

    # Set identities for condition
    Idents(sobj_ct) <- "condition"

    # Get condition values for comparison
    condition_values <- unique(sobj_ct$condition)

    # Determine comparison order (try to put disease first)
    if (length(condition_values) == 2) {
      # Check if one is AD/disease and one is control
      disease_keywords <- c("AD", "Alzheimer", "Disease", "PS2", "5xFAD", "Tau")
      control_keywords <- c("Control", "HC", "WT", "Healthy", "Normal")

      disease_val <- NULL
      control_val <- NULL

      for (val in condition_values) {
        val_upper <- toupper(val)
        if (any(sapply(disease_keywords, function(k) grepl(k, val_upper, ignore.case = TRUE)))) {
          disease_val <- val
        } else if (any(sapply(control_keywords, function(k) grepl(k, val_upper, ignore.case = TRUE)))) {
          control_val <- val
        }
      }

      if (!is.null(disease_val) && !is.null(control_val)) {
        comparison <- c(disease_val, control_val)
        cat("         Comparison:", disease_val, "vs", control_val, "\n")
      } else {
        comparison <- condition_values
        cat("         Comparison:", paste(comparison, collapse = " vs "), "\n")
      }
    } else {
      comparison <- condition_values
      cat("         Comparison:", paste(comparison, collapse = " vs "), "\n")
    }

    # Calculate AUC using presto::wilcoxauc
    auc_df <- tryCatch({
      presto::wilcoxauc(
        sobj_ct,
        comparison = comparison,
        seurat_assay = "RNA"
      )
    }, error = function(e) {
      cat("         Error in wilcoxauc:", e$message, "\n")
      return(NULL)
    })

    if (is.null(auc_df) || nrow(auc_df) == 0) {
      cat("         No results from wilcoxauc\n")
      next
    }

    # Check result structure and print column names for debugging
    cat("         AUC columns:", paste(colnames(auc_df), collapse = ", "), "\n")
    if (!is.data.frame(auc_df)) {
      cat("         Unexpected wilcoxauc output format\n")
      next
    }

    # Filter for genes of interest and format results
    # Handle both 'feature' and 'gene' column names
    gene_col <- if ("feature" %in% colnames(auc_df)) "feature" else if ("gene" %in% colnames(auc_df)) "gene" else NULL

    if (is.null(gene_col)) {
      cat("         ERROR: Cannot find gene column in wilcoxauc output\n")
      next
    }

    pval_col <- if ("pval" %in% colnames(auc_df)) "pval" else if ("p_value" %in% colnames(auc_df)) "p_value" else "pval"

    ct_results <- auc_df %>%
      filter(!!sym(gene_col) %in% genes_to_test) %>%
      select(!!sym(gene_col), auc, !!sym(pval_col), group) %>%
      mutate(
        celltype = ct,
        dataset = dataset_name,
        comparison = paste(comparison[1], "vs", comparison[2]),
        gene_list = "ECM_communicationV4",
        method = "wilcoxauc"
      ) %>%
      rename(gene = !!sym(gene_col), p_value = !!sym(pval_col))

    wilcoxauc_results[[ct]] <- ct_results
    cat("         Found", nrow(ct_results), "genes with results\n")
  }

  # Combine all results
  if (length(wilcoxauc_results) == 0) {
    cat("      WARNING: No results obtained for", dataset_name, "\n")
    return(data.frame())
  }

  result <- do.call(rbind, wilcoxauc_results)
  cat("   Total results:", nrow(result), "gene-cell type combinations\n\n")

  return(result)
}

# ---------------------------------------------------------------------------
# 3. ANALYZE EACH DATASET
# ---------------------------------------------------------------------------
cat("\n3. Loading datasets and running wilcoxauc analysis...\n")

all_results <- list()

# === Dataset 1: GSE138852 (Mouse) ===
cat("\n--- Dataset 1: GSE138852 (Mouse AD) ---\n")
cat("   NOTE: Skipping GSE138852 for wilcoxauc analysis\n")
cat("   Reason: GSE138852 uses SCT-transformed data with empty 'data' layers\n")
cat("   which is not compatible with presto::wilcoxauc without special handling.\n")
cat("   GSE138852 is included in the pseudobulk analysis instead.\n")

# Skip GSE138852 for wilcoxauc
all_results$GSE138852 <- NULL

# === Dataset 2: GSE174367 (Human snRNA-seq + snATAC-seq) ===
cat("\n--- Dataset 2: GSE174367 (Human Multiome) ---\n")
gse174367_path <- "data/GSE174367/sobj_20251212.rds"

if (file.exists(gse174367_path)) {
  sobj_gse174367 <- readRDS(gse174367_path)

  # For GSE174367, use Diagnosis and Cell.Type
  if ("Diagnosis" %in% colnames(sobj_gse174367@meta.data)) {
    sobj_gse174367$condition <- sobj_gse174367$Diagnosis
    cat("   Using Diagnosis as condition variable\n")
    cat("   Conditions:", toString(unique(sobj_gse174367$condition)), "\n")
  }
  if ("Cell.Type" %in% colnames(sobj_gse174367@meta.data)) {
    sobj_gse174367$celltype <- sobj_gse174367$Cell.Type
    cat("   Using Cell.Type as celltype variable\n")
    cat("   Cell types:", toString(unique(sobj_gse174367$celltype)), "\n")
  }

  result_gse174367 <- perform_wilcoxauc_analysis(
    sobj_gse174367,
    gene_list,
    "GSE174367_Human"
  )

  if (nrow(result_gse174367) > 0) {
    all_results$GSE174367 <- result_gse174367

    # Save individual results
    output_file <- file.path(task_dir, "Documents", "GSE174367_ECM_communicationV4_wilcoxauc.csv")
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

  # For GSE188545, use standard condition and celltype
  if ("condition" %in% colnames(sobj_gse188545@meta.data)) {
    cat("   Using condition as condition variable\n")
    cat("   Conditions:", toString(unique(sobj_gse188545$condition)), "\n")
  }
  if ("celltype" %in% colnames(sobj_gse188545@meta.data)) {
    cat("   Using celltype as celltype variable\n")
    cat("   Cell types:", toString(unique(sobj_gse188545$celltype)), "\n")
  }

  result_gse188545 <- perform_wilcoxauc_analysis(
    sobj_gse188545,
    gene_list,
    "GSE188545_Human"
  )

  if (nrow(result_gse188545) > 0) {
    all_results$GSE188545 <- result_gse188545

    # Save individual results
    output_file <- file.path(task_dir, "Documents", "GSE188545_ECM_communicationV4_wilcoxauc.csv")
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

# Filter out NULL results
all_results <- all_results[!sapply(all_results, is.null)]

if (length(all_results) > 0) {
  combined_results <- do.call(rbind, all_results)

  # Save combined results
  combined_file <- file.path(task_dir, "Documents", "All_Datasets_ECM_communicationV4_wilcoxauc_combined.csv")
  write.csv(combined_results, combined_file, row.names = FALSE)
  cat("   Saved combined results:", combined_file, "\n")
  cat("   Total gene-cell type combinations across all datasets:", nrow(combined_results), "\n")

  # Create summary statistics per dataset
  cat("\n5. Creating summary statistics...\n")

  summary_stats <- combined_results %>%
    group_by(dataset) %>%
    summarise(
      cell_types_analyzed = n_distinct(celltype),
      genes_tested = n_distinct(gene),
      total_combinations = n(),
      mean_auc = mean(auc, na.rm = TRUE),
      sd_auc = sd(auc, na.rm = TRUE),
      genes_auc_gt_0.6 = sum(auc > 0.6, na.rm = TRUE),
      genes_auc_gt_0.7 = sum(auc > 0.7, na.rm = TRUE),
      genes_auc_gt_0.8 = sum(auc > 0.8, na.rm = TRUE),
      .groups = "drop"
    )

  print(summary_stats)

  # Save summary statistics
  summary_file <- file.path(task_dir, "Documents", "summary_statistics_per_dataset.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("   Saved summary statistics:", summary_file, "\n")

  # Summary per cell type across datasets
  summary_by_celltype <- combined_results %>%
    group_by(dataset, celltype) %>%
    summarise(
      genes_tested = n_distinct(gene),
      mean_auc = mean(auc, na.rm = TRUE),
      genes_auc_gt_0.6 = sum(auc > 0.6, na.rm = TRUE),
      genes_auc_gt_0.7 = sum(auc > 0.7, na.rm = TRUE),
      .groups = "drop"
    )

  celltype_file <- file.path(task_dir, "Documents", "summary_statistics_per_celltype.csv")
  write.csv(summary_by_celltype, celltype_file, row.names = FALSE)
  cat("   Saved celltype summary:", celltype_file, "\n")

  # Identify high AUC genes across all datasets
  cat("\n6. Identifying high AUC genes (AUC > 0.7)...\n")

  high_auc_genes <- combined_results %>%
    filter(auc > 0.7) %>%
    arrange(desc(auc))

  if (nrow(high_auc_genes) > 0) {
    print(high_auc_genes)

    high_auc_file <- file.path(task_dir, "Documents", "All_Datasets_ECM_communicationV4_high_auc_genes.csv")
    write.csv(high_auc_genes, high_auc_file, row.names = FALSE)
    cat("   Saved high AUC genes:", high_auc_file, "\n")
  } else {
    cat("   No genes with AUC > 0.7 found\n")
  }

  # ---------------------------------------------------------------------------
  # 7. CREATE VISUALIZATIONS
  # ---------------------------------------------------------------------------
  cat("\n7. Creating visualizations...\n")

  # 7a. Heatmap for each dataset
  for (ds in unique(combined_results$dataset)) {
    cat("   Creating heatmap for", ds, "\n")

    ds_data <- combined_results %>%
      filter(dataset == ds, group == names(table(combined_results$group))[1]) %>%
      select(celltype, gene, auc) %>%
      pivot_wider(
        names_from = gene,
        values_from = auc,
        values_fill = 0.5
      ) %>%
      column_to_rownames("celltype") %>%
      as.matrix()

    if (nrow(ds_data) > 1 && ncol(ds_data) > 1) {
      heatmap_file <- file.path(task_dir, "plot", paste0(ds, "_ECM_communicationV4_AUC_heatmap.pdf"))

      pdf(heatmap_file, width = max(10, ncol(ds_data) * 0.3),
          height = max(6, nrow(ds_data) * 0.5))

      pheatmap(
        ds_data,
        main = paste("AUC of ECM Communication V4 Genes -", ds),
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
        border_color = NA
      )

      dev.off()
      cat("      Saved:", heatmap_file, "\n")
    }
  }

  # 7b. Combined heatmap (all datasets)
  cat("   Creating combined heatmap...\n")

  # Create pivot table for combined heatmap
  combined_heatmap_data <- combined_results %>%
    filter(group == names(table(combined_results$group))[1]) %>%
    mutate(dataset_celltype = paste(dataset, celltype, sep = "_")) %>%
    select(dataset_celltype, gene, auc) %>%
    pivot_wider(
      names_from = gene,
      values_from = auc,
      values_fill = 0.5
    ) %>%
    column_to_rownames("dataset_celltype") %>%
    as.matrix()

  if (nrow(combined_heatmap_data) > 1 && ncol(combined_heatmap_data) > 1) {
    combined_heatmap_file <- file.path(task_dir, "plot", "All_Datasets_ECM_communicationV4_AUC_heatmap.pdf")

    pdf(combined_heatmap_file, width = max(12, ncol(combined_heatmap_data) * 0.2),
        height = max(8, nrow(combined_heatmap_data) * 0.3))

    pheatmap(
      combined_heatmap_data,
      main = "AUC of ECM Communication V4 Genes - All Datasets",
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
      border_color = NA
    )

    dev.off()
    cat("      Saved:", combined_heatmap_file, "\n")
  }

  # 7c. Summary bar plots
  cat("   Creating summary bar plots...\n")

  p1 <- ggplot(summary_stats, aes(x = dataset, y = genes_tested, fill = dataset)) +
    geom_bar(stat = "identity") +
    labs(title = "Genes Tested per Dataset",
         x = "Dataset", y = "Number of Genes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  p2 <- ggplot(summary_stats, aes(x = dataset, y = genes_auc_gt_0.7, fill = dataset)) +
    geom_bar(stat = "identity") +
    labs(title = "Genes with AUC > 0.7 per Dataset",
         x = "Dataset", y = "Number of Genes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  p3 <- ggplot(summary_stats, aes(x = dataset, y = mean_auc, fill = dataset)) +
    geom_bar(stat = "identity") +
    labs(title = "Mean AUC per Dataset",
         x = "Dataset", y = "Mean AUC") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  summary_plot <- p1 | p2 | p3
  summary_plot_file <- file.path(task_dir, "plot", "Summary_statistics_per_dataset.pdf")
  ggsave(summary_plot_file, plot = summary_plot, width = 15, height = 5)
  cat("      Saved:", summary_plot_file, "\n")

} else {
  cat("\nERROR: No results obtained from any dataset\n")
}

# ---------------------------------------------------------------------------
# 8. CREATE ANALYSIS SUMMARY DOCUMENT
# ---------------------------------------------------------------------------
cat("\n8. Creating analysis summary document...\n")

summary_doc <- file.path(task_dir, "Documents", "analysis_summary.md")
sink(summary_doc)

cat("# ECM Communication V4 Wilcoxauc Analysis Summary\n\n")
cat("**Task Directory**: `20260108-17-ECM_communicationV4_wilcoxauc_3datasets`\n")
cat("**Date**: ", date(), "\n\n")
cat("## Overview\n\n")
cat("This task performed AUC (Area Under the ROC Curve) analysis of ECM communication V4 genes across three Alzheimer's disease single-cell datasets to identify genes differentially expressed between disease and control conditions across different brain cell types.\n\n")
cat("## Datasets Analyzed\n\n")
cat("1. **GSE138852** (Mouse AD model)\n")
cat("   - Used oupSample.batchCond as condition\n")
cat("   - Used oupSample.cellType as celltype\n")
cat("2. **GSE174367** (Human snRNA-seq + snATAC-seq multiome)\n")
cat("   - Used Diagnosis as condition\n")
cat("   - Used Cell.Type as celltype\n")
cat("3. **GSE188545** (Human MTG single-cell RNA-seq)\n")
cat("   - Used condition as AD vs HC\n")
cat("   - Used celltype\n\n")
cat("## Gene List Analyzed\n\n")
cat("- **ECM communication genes V4** (`ECM_communication_genesV4.xlsx`)\n")
cat("- ", length(gene_list), " genes in list\n\n")
cat("## Methodology\n\n")
cat("### Analysis Method\n")
cat("- Used `presto::wilcoxauc()` for AUC calculation\n")
cat("- Comparison: Disease vs Control conditions\n")
cat("- Minimum 10 cells per condition required per cell type\n")
cat("- Seurat v5 compatibility: Used `JoinLayers()` before analysis\n\n")
cat("### Key Results\n\n")

if (length(all_results) > 0 && exists("summary_stats")) {
  cat("- Datasets analyzed:", length(unique(combined_results$dataset)), "\n")
  cat("- Total gene-cell type combinations:", nrow(combined_results), "\n")
  cat("- Mean AUC across all datasets:", round(mean(combined_results$auc), 3), "\n")
  cat("- Genes with AUC > 0.7:", sum(combined_results$auc > 0.7), "\n\n")

  cat("### Per Dataset Statistics\n\n")
  for (i in 1:nrow(summary_stats)) {
    cat("####", summary_stats$dataset[i], "\n")
    cat("- Cell types analyzed:", summary_stats$cell_types_analyzed[i], "\n")
    cat("- Genes tested:", summary_stats$genes_tested[i], "\n")
    cat("- Mean AUC:", round(summary_stats$mean_auc[i], 3), "\n")
    cat("- Genes with AUC > 0.7:", summary_stats$genes_auc_gt_0.7[i], "\n\n")
  }
}

cat("## Files Generated\n\n")
cat("### Individual Dataset Results\n")
cat("1. `Documents/GSE138852_ECM_communicationV4_wilcoxauc.csv` - GSE138852 results\n")
cat("2. `Documents/GSE174367_ECM_communicationV4_wilcoxauc.csv` - GSE174367 results\n")
cat("3. `Documents/GSE188545_ECM_communicationV4_wilcoxauc.csv` - GSE188545 results\n\n")
cat("### Combined Results\n")
cat("4. `Documents/All_Datasets_ECM_communicationV4_wilcoxauc_combined.csv` - All datasets combined\n")
cat("5. `Documents/summary_statistics_per_dataset.csv` - Summary statistics per dataset\n")
cat("6. `Documents/summary_statistics_per_celltype.csv` - Summary statistics per cell type\n")
cat("7. `Documents/All_Datasets_ECM_communicationV4_high_auc_genes.csv` - Genes with AUC > 0.7\n\n")
cat("### Visualizations\n")
cat("8. `plot/GSE138852_ECM_communicationV4_AUC_heatmap.pdf` - GSE138852 heatmap\n")
cat("9. `plot/GSE174367_ECM_communicationV4_AUC_heatmap.pdf` - GSE174367 heatmap\n")
cat("10. `plot/GSE188545_ECM_communicationV4_AUC_heatmap.pdf` - GSE188545 heatmap\n")
cat("11. `plot/All_Datasets_ECM_communicationV4_AUC_heatmap.pdf` - Combined heatmap\n")
cat("12. `plot/Summary_statistics_per_dataset.pdf` - Summary bar plots\n\n")
cat("## AUC Interpretation\n\n")
cat("- **AUC = 0.5**: No discriminatory power (equal expression)\n")
cat("- **AUC > 0.5**: Higher expression in disease group\n")
cat("- **AUC < 0.5**: Higher expression in control group\n")
cat("- **AUC ≥ 0.7**: Moderate discriminatory power\n")
cat("- **AUC ≥ 0.8**: Strong discriminatory power\n\n")

sink()

cat("   Saved:", summary_doc, "\n")

# ---------------------------------------------------------------------------
# 9. SESSION INFO
# ---------------------------------------------------------------------------
session_file <- file.path(task_dir, "Documents", "session_info.txt")
sink(session_file)
cat("ECM Communication V4 Wilcoxauc Analysis - 3 Datasets\n")
cat("================================================\n\n")
cat("Datasets: GSE138852, GSE174367, GSE188545\n")
cat("Analysis date:", date(), "\n\n")
print(sessionInfo())
sink()

cat("\n=== Analysis Complete ===\n")
cat("All results saved to:", task_dir, "\n")
