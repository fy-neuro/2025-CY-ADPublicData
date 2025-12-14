# Differential Expression Analysis for GSE188545
# AD vs HC comparison per annotated cell type
# Uses annotated Seurat object from fast_annotation.R

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(here)

# Set working directory to project root
project_root <- here::here()
setwd(project_root)

cat("=== GSE188545 Differential Expression Analysis (AD vs HC) ===\n")

# Input: annotated object from fast annotation
sobj_path <- "task/20251214-10-GSE188545_Annotation/GSE188545_sobj_annotated_fast.rds"
if (!file.exists(sobj_path)) {
  stop("Annotated object not found at: ", sobj_path)
}

cat("Loading annotated object...\n")
sobj <- readRDS(sobj_path)
cat("Object loaded:", ncol(sobj), "cells,", nrow(sobj), "genes\n")
cat("Cell types:", toString(unique(sobj$celltype)), "\n")
cat("Conditions:", toString(unique(sobj$condition)), "\n")

# Create output directories
output_dir <- "task/20251214-10-GSE188545_Annotation"
dir.create(file.path(output_dir, "plot/differential_expression"), 
           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "Documents/differential_expression"), 
           showWarnings = FALSE, recursive = TRUE)

# Subset to cell types with at least 50 cells in both conditions
celltype_counts <- sobj@meta.data %>%
  group_by(celltype, condition) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = n_cells, values_fill = 0)

celltypes_for_de <- celltype_counts %>%
  filter(AD >= 50 & HC >= 50) %>%
  pull(celltype)

cat("\nCell types with sufficient cells for DE analysis:", length(celltypes_for_de), "\n")
cat(paste(celltypes_for_de, collapse = ", "), "\n")

if (length(celltypes_for_de) == 0) {
  stop("No cell types with sufficient cells in both conditions for DE analysis.")
}

de_results <- list()

for (ct in celltypes_for_de) {
  cat("  Analyzing", ct, "...\n")
  
  # Subset to cell type
  sobj_subset <- subset(sobj, subset = celltype == ct)
  
  # Join layers for Seurat v5 compatibility (required for FindMarkers)
  sobj_subset <- JoinLayers(sobj_subset, assay = "RNA", layer = "data")
  
  # Set identity to condition
  Idents(sobj_subset) <- "condition"
  
  # Find markers between AD and HC
  markers_ct <- FindMarkers(sobj_subset,
                            ident.1 = "AD",
                            ident.2 = "HC",
                            assay = "RNA",
                            logfc.threshold = 0.1,
                            min.pct = 0.1,
                            verbose = FALSE)
  
  # Add gene symbols and cell type
  markers_ct$gene <- rownames(markers_ct)
  markers_ct$celltype <- ct
  markers_ct$direction <- ifelse(markers_ct$avg_log2FC > 0, "UP_in_AD", "DOWN_in_AD")
  
  # Save individual results
  write.csv(markers_ct,
            file.path(output_dir, "Documents/differential_expression", 
                     paste0("DE_", gsub("/", "_", ct), "_AD_vs_HC.csv")),
            row.names = FALSE)
  
  de_results[[ct]] <- markers_ct
  
  # Volcano plot
  if (nrow(markers_ct) > 0) {
    # Select top up/down regulated genes for labeling
    top_up <- markers_ct %>%
      filter(avg_log2FC > 0) %>%
      arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
      head(10)
    
    top_down <- markers_ct %>%
      filter(avg_log2FC < 0) %>%
      arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
      head(10)
    
    volcano <- ggplot(markers_ct, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(aes(color = direction), alpha = 0.6, size = 1) +
      geom_point(data = rbind(top_up, top_down), 
                 aes(x = avg_log2FC, y = -log10(p_val_adj)), 
                 color = "red", size = 2) +
      geom_text_repel(data = rbind(top_up, top_down),
                     aes(label = gene), size = 3, max.overlaps = 20) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
      scale_color_manual(values = c("UP_in_AD" = "red", "DOWN_in_AD" = "blue")) +
      labs(x = "Log2 Fold Change (AD vs HC)", 
           y = "-log10(Adjusted p-value)",
           title = paste("Differential Expression:", ct, "AD vs HC")) +
      theme_minimal()
    
    ggsave(file.path(output_dir, "plot/differential_expression", 
                    paste0("volcano_", gsub("/", "_", ct), "_AD_vs_HC.pdf")),
           plot = volcano, width = 10, height = 8)
  }
}

# Combine all DE results
de_combined <- bind_rows(de_results)
write.csv(de_combined,
          file.path(output_dir, "Documents/differential_expression/DE_all_celltypes_AD_vs_HC.csv"),
          row.names = FALSE)

cat("\n=== Differential expression analysis complete ===\n")
cat("Results saved to:", file.path(output_dir, "Documents/differential_expression"), "\n")
cat("Volcano plots saved to:", file.path(output_dir, "plot/differential_expression"), "\n")

# Session info
session_file <- file.path(output_dir, "Documents/differential_expression/session_info_de.txt")
sink(session_file)
cat("Differential Expression Analysis Session Information\n")
cat("====================================================\n\n")
cat("Analysis date:", date(), "\n\n")
print(sessionInfo())
sink()
cat("Session info saved to:", session_file, "\n")