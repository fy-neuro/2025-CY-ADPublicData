# Fast annotation pipeline for GSE188545
# Performs basic cell type annotation without DE or GO analysis

cat("=== GSE188545 Fast Annotation Pipeline ===\n")

# Load required libraries
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

# Set working directory
project_root <- here::here()
setwd(project_root)
cat("Working directory:", getwd(), "\n")

# Configuration
input_dir <- "task/20251214-9-GSE188545-Load_Data/"
output_dir <- "task/20251214-10-GSE188545_Annotation/"
integrated_object <- file.path(input_dir, "GSE188545_sobj_integrated.rds")
annotated_object <- file.path(output_dir, "GSE188545_sobj_annotated_fast.rds")

# Create output directories
dir.create(file.path(output_dir, "plot/annotation_fast"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "Documents/annotation_fast"), showWarnings = FALSE, recursive = TRUE)

# Step 1: Load integrated object
cat("\n1. Loading integrated object...\n")
sobj <- readRDS(integrated_object)
cat("  Cells:", ncol(sobj), "\n")
cat("  Clusters (res 0.5):", length(unique(sobj$RNA_snn_res.0.5)), "\n")
cat("  AD cells:", sum(sobj$condition == "AD"), "\n")
cat("  HC cells:", sum(sobj$condition == "HC"), "\n")

# Set cluster identity
Idents(sobj) <- "RNA_snn_res.0.5"

# Step 2: Load marker genes
cat("\n2. Loading cell type markers...\n")
# Note: src/1stAnnotation.R contains mouse markers; using human markers directly
# source(here("src/1stAnnotation.R"), encoding = "UTF-8")

# Simplified cell type markers
celltype_markers <- list(
  Microglia = c("P2RY12", "TMEM119", "CX3CR1", "CSF1R", "HEXB", "ITGAM", "CD74", "C3"),
  Astrocytes = c("AQP4", "GFAP", "FGFR3", "NHSL1", "SLC25A18"),
  Oligodendrocytes = c("MBP", "MOBP", "MOG", "OPALIN"),
  OPC = c("PDGFRA", "VCAN", "OLIG1"),
  Endothelial = c("DCN", "FLT1", "LEF1", "VWF"),
  Glutamatergic_Neurons = c("RALYL", "LDB2", "NELL2"),
  GABAergic_Neurons = c("GAD1", "GAD2", "SLC6A1", "GRIP1"),
  T_cells = c("CD2", "THEMIS", "CD3D"),
  Macrophages = c("CD68", "ITGAM", "CD14"),
  Proliferating = c("MKI67", "PCNA", "TOP2A")
)

# Step 3: Score-based annotation
cat("\n3. Performing score-based annotation...\n")

score_celltype <- function(markers, obj) {
  # Get variable features (scale.data contains only variable genes)
  var_features <- rownames(GetAssayData(obj, assay = "RNA", layer = "scale.data"))
  markers_avail <- intersect(markers, var_features)
  if (length(markers_avail) < 2) return(rep(0, ncol(obj)))
  expr <- GetAssayData(obj, assay = "RNA", layer = "scale.data")
  expr_subset <- expr[markers_avail, , drop = FALSE]
  colMeans(expr_subset)
}

# Calculate scores
for (celltype in names(celltype_markers)) {
  sobj[[paste0("score_", celltype)]] <- score_celltype(celltype_markers[[celltype]], sobj)
}

# Assign preliminary annotation
celltype_scores <- sobj@meta.data %>%
  select(starts_with("score_")) %>%
  as.matrix()

celltype_assignments <- apply(celltype_scores, 1, function(x) {
  celltypes <- names(celltype_markers)
  if (max(x) < 0.5) return("Unknown")
  celltypes[which.max(x)]
})

sobj$celltype_prelim <- celltype_assignments

# Refine by cluster consensus
cluster_annotations <- sobj@meta.data %>%
  group_by(RNA_snn_res.0.5, celltype_prelim) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(RNA_snn_res.0.5) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  filter(prop >= 0.3) %>%
  arrange(RNA_snn_res.0.5, desc(prop))

annotation_map <- cluster_annotations %>%
  group_by(RNA_snn_res.0.5) %>%
  slice(1) %>%
  select(cluster = RNA_snn_res.0.5, celltype = celltype_prelim)

# Apply annotation
sobj$celltype <- plyr::mapvalues(sobj$RNA_snn_res.0.5,
                                 from = annotation_map$cluster,
                                 to = annotation_map$celltype)

unannotated <- sobj$celltype == "" | is.na(sobj$celltype)
sobj$celltype[unannotated] <- sobj$celltype_prelim[unannotated]

# Save mapping
write.csv(annotation_map,
          file.path(output_dir, "Documents/annotation_fast/cluster_annotation_mapping.csv"),
          row.names = FALSE)

cat("  Annotation complete.\n")
cat("  Cell type distribution:\n")
print(table(sobj$celltype))

# Step 4: Basic visualizations
cat("\n4. Generating visualizations...\n")

# UMAP by cell type
umap_celltype <- DimPlot(sobj,
                         reduction = "umap",
                         group.by = "celltype",
                         label = TRUE,
                         repel = TRUE,
                         pt.size = 0.1,
                         raster = TRUE) +
  ggtitle("Cell Type Annotation (Fast)")

ggsave(file.path(output_dir, "plot/annotation_fast/umap_celltype.pdf"),
       plot = umap_celltype, width = 10, height = 8)

# Dot plot of key markers
key_markers <- c("P2RY12", "TMEM119", "AQP4", "GFAP", "MBP", "MOG", 
                 "PDGFRA", "VWF", "GAD1", "NELL2", "CD3D")
key_markers_avail <- key_markers[key_markers %in% rownames(sobj)]

if (length(key_markers_avail) > 0) {
  dotplot <- DotPlot(sobj,
                     features = key_markers_avail,
                     group.by = "celltype",
                     dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Marker Expression by Cell Type")
  
  ggsave(file.path(output_dir, "plot/annotation_fast/dotplot_key_markers.pdf"),
         plot = dotplot, width = max(10, length(key_markers_avail) * 0.4),
         height = max(6, length(unique(sobj$celltype)) * 0.5))
}

# Step 5: Save annotated object
cat("\n5. Saving annotated object...\n")
saveRDS(sobj, file = annotated_object)
cat("  Saved to:", annotated_object, "\n")

# Step 6: Summary
cat("\n=== Fast Annotation Complete ===\n")
cat("Total cells:", ncol(sobj), "\n")
cat("Cell types:", length(unique(sobj$celltype)), "\n")
cat("Annotated cells:", sum(sobj$celltype != "Unknown"), "\n")

celltype_counts <- table(sobj$celltype)
for (ct in names(celltype_counts)) {
  cat("  ", ct, ":", celltype_counts[ct], "cells\n")
}

cat("\nOutput directory:", output_dir, "\n")
cat("Next: Run full annotation pipeline for DE and GO analysis.\n")