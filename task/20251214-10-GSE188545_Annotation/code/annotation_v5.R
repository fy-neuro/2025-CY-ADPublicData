# Cell Type Annotation for GSE188545
# Annotates clusters using canonical marker genes

# Load required libraries
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)
library(viridis)

# Set working directory to project root
project_root <- here::here()
setwd(project_root)
cat("Working directory set to:", getwd(), "\n")

# Configuration
input_dir <- "task/20251214-9-GSE188545-Load_Data/"
output_dir <- "task/20251214-10-GSE188545_Annotation/"
integrated_object <- file.path(input_dir, "GSE188545_sobj_integrated.rds")
annotated_object <- file.path(output_dir, "GSE188545_sobj_annotated.rds")

cat("=== GSE188545 Cell Type Annotation ===\n")
cat("Input directory:", input_dir, "\n")

# Create output directories
dir.create(file.path(output_dir, "plot/annotation"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "Documents/annotation"), showWarnings = FALSE, recursive = TRUE)

# Step 1: Load integrated object
cat("\n1. Loading integrated object...\n")
if (!file.exists(integrated_object)) {
  stop("Integrated object not found. Run 2.0QC_Filtering_Integration.R first.")
}

sobj <- readRDS(integrated_object)
cat("  Loaded object:", sobj@project.name, "\n")
cat("  Cells:", ncol(sobj), "\n")
cat("  Clusters (res 0.5):", length(unique(sobj$RNA_snn_res.0.5)), "\n")

# Ensure condition metadata is correct
if ("sample" %in% colnames(sobj@meta.data)) {
  sobj$condition <- ifelse(grepl("_AD", sobj$sample), "AD",
                           ifelse(grepl("_HC", sobj$sample), "HC", "Unknown"))
}
cat("  AD cells:", sum(sobj$condition == "AD"), "\n")
cat("  HC cells:", sum(sobj$condition == "HC"), "\n")

# Set default cluster identity for annotation
Idents(sobj) <- "RNA_snn_res.0.5"
clusters <- levels(Idents(sobj))
cat("  Cluster IDs:", clusters, "\n")

# Step 2: Load cell type marker genes
cat("\n2. Loading cell type marker genes...\n")

# Source the annotation markers from project
source(here("src/1stAnnotation.R"), encoding = "UTF-8")

# Define cell type markers based on CNS cell types
# Using markers from src/1stAnnotation.R
celltype_markers <- list(
  # Microglia markers
  Microglia = c("P2RY12", "TMEM119", "CX3CR1", "CSF1R", "HEXB", "ITGAM", "CD74", "C3"),
  
  # Astrocyte markers
  Astrocytes = c("AQP4", "GFAP", "FGFR3", "NHSL1", "SLC25A18"),
  
  # Oligodendrocyte markers
  Oligodendrocytes = c("MBP", "MOBP", "MOG", "OPALIN"),
  
  # OPC (Oligodendrocyte Precursor Cells)
  OPC = c("PDGFRA", "VCAN", "OLIG1"),
  
  # Endothelial cells (VC = Vascular Cells)
  Endothelial = c("DCN", "FLT1", "LEF1", "VWF"),
  
  # Neurons
  Glutamatergic_Neurons = c("RALYL", "LDB2", "NELL2"),
  GABAergic_Neurons = c("GAD1", "GAD2", "SLC6A1", "GRIP1"),
  
  # T cells
  T_cells = c("CD2", "THEMIS", "CD3D"),
  
  # Additional immune cells
  Macrophages = c("CD68", "ITGAM", "CD14"),
  
  # Proliferating cells
  Proliferating = c("MKI67", "PCNA", "TOP2A")
)

# Create a combined vector of all markers for heatmaps
all_markers <- unique(unlist(celltype_markers))

# Check which markers are present in the dataset
available_markers <- all_markers[all_markers %in% rownames(sobj)]
missing_markers <- setdiff(all_markers, rownames(sobj))

cat("  Available markers:", length(available_markers), "/", length(all_markers), "\n")
if (length(missing_markers) > 0) {
  cat("  Missing markers (", length(missing_markers), "): ", 
      paste(missing_markers[1:min(10, length(missing_markers))], collapse = ", "), 
      ifelse(length(missing_markers) > 10, "...", ""), "\n")
}

# Step 3: Find cluster markers
cat("\n3. Finding cluster markers...\n")

# Find all markers for each cluster
cat("  Finding differentially expressed genes...\n")
# For Seurat v5, need to join layers for differential expression
cat("    Joining data layers...\n")
sobj_joined <- JoinLayers(sobj, assay = "RNA", layer = "data")
cluster_markers <- FindAllMarkers(sobj_joined,
                                  assay = "RNA",
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25,
                                  verbose = FALSE)

# Save cluster markers
write.csv(cluster_markers,
          file.path(output_dir, "Documents/annotation/cluster_markers_all.csv"),
          row.names = FALSE)

# Extract top markers per cluster
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup()

write.csv(top_markers,
          file.path(output_dir, "Documents/annotation/cluster_markers_top10.csv"),
          row.names = FALSE)

cat("  Found", nrow(cluster_markers), "marker genes across", 
    length(unique(cluster_markers$cluster)), "clusters\n")

# Step 4: Annotate clusters based on marker expression
cat("\n4. Annotating clusters...\n")

# Function to score clusters for each cell type (Seurat v5 compatible)
score_celltype <- function(markers, obj) {
  # Subset to available markers
  markers_avail <- markers[markers %in% rownames(obj)]
  if (length(markers_avail) < 2) return(rep(0, ncol(obj)))
  
  # Calculate average expression of markers using joined data layer
  # First join layers to get combined expression matrix
  obj_joined <- JoinLayers(obj, assay = "RNA", layer = "data")
  expr <- GetAssayData(obj_joined, assay = "RNA", layer = "data")
  expr_subset <- expr[markers_avail, , drop = FALSE]
  scores <- colMeans(expr_subset)
  return(scores)
}

# Calculate scores for each cell type
cat("  Calculating cell type scores...\n")
for (celltype in names(celltype_markers)) {
  sobj[[paste0("score_", celltype)]] <- score_celltype(celltype_markers[[celltype]], sobj)
}

# Assign preliminary annotation based on highest score
celltype_scores <- sobj@meta.data %>%
  select(starts_with("score_")) %>%
  as.matrix()

# Find maximum score for each cell
max_scores <- apply(celltype_scores, 1, max)
celltype_assignments <- apply(celltype_scores, 1, function(x) {
  celltypes <- names(celltype_markers)
  if (max(x) < 0.5) return("Unknown")
  celltypes[which.max(x)]
})

sobj$celltype_prelim <- celltype_assignments

# Refine annotation by cluster consensus
cluster_annotations <- sobj@meta.data %>%
  group_by(RNA_snn_res.0.5, celltype_prelim) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(RNA_snn_res.0.5) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  filter(prop >= 0.3) %>%  # At least 30% of cluster agrees
  arrange(RNA_snn_res.0.5, desc(prop))

# Assign final annotation
annotation_map <- cluster_annotations %>%
  group_by(RNA_snn_res.0.5) %>%
  slice(1) %>%  # Take top cell type for each cluster
  select(cluster = RNA_snn_res.0.5, celltype = celltype_prelim)

# Apply annotation
sobj$celltype <- plyr::mapvalues(sobj$RNA_snn_res.0.5,
                                 from = annotation_map$cluster,
                                 to = annotation_map$celltype)

# For clusters without clear annotation, keep preliminary assignment
unannotated <- sobj$celltype == "" | is.na(sobj$celltype)
sobj$celltype[unannotated] <- sobj$celltype_prelim[unannotated]

# Save annotation mapping
write.csv(annotation_map,
          file.path(output_dir, "Documents/annotation/cluster_annotation_mapping.csv"),
          row.names = FALSE)

cat("  Annotation complete.\n")
cat("  Cell type distribution:\n")
print(table(sobj$celltype))

# Step 5: Visualize annotations
cat("\n5. Generating annotation visualizations...\n")

# UMAP by cell type
umap_celltype <- DimPlot(sobj,
                         reduction = "umap",
                         group.by = "celltype",
                         label = TRUE,
                         repel = TRUE,
                         pt.size = 0.1) +
  ggtitle("Cell Type Annotation") +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "plot/annotation/umap_celltype.pdf"),
       plot = umap_celltype, width = 10, height = 8)

# UMAP by cluster with cell type labels
umap_cluster_annotated <- DimPlot(sobj,
                                  reduction = "umap",
                                  group.by = "RNA_snn_res.0.5",
                                  label = TRUE,
                                  pt.size = 0.1) +
  ggtitle("Clusters with Cell Type Annotation")

# Add cell type labels to cluster plot
cluster_centers <- sobj@meta.data %>%
  group_by(RNA_snn_res.0.5, celltype) %>%
  summarise(UMAP_1 = median(sobj@reductions$umap@cell.embeddings[,1]),
            UMAP_2 = median(sobj@reductions$umap@cell.embeddings[,2]),
            .groups = "drop")

umap_cluster_annotated <- umap_cluster_annotated +
  geom_text(data = cluster_centers,
            aes(x = UMAP_1, y = UMAP_2, label = celltype),
            size = 3, fontface = "bold", color = "darkred")

ggsave(file.path(output_dir, "plot/annotation/umap_clusters_with_celltypes.pdf"),
       plot = umap_cluster_annotated, width = 10, height = 8)

# Dot plot of marker expression
cat("  Generating marker expression dot plot...\n")
dotplot_markers <- DotPlot(sobj,
                           assay = "RNA",
                           features = available_markers,
                           group.by = "celltype",
                           cols = c("lightgrey", "blue"),
                           dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marker Expression by Cell Type")

ggsave(file.path(output_dir, "plot/annotation/dotplot_markers.pdf"),
       plot = dotplot_markers, width = max(10, length(available_markers) * 0.3),
       height = max(6, length(unique(sobj$celltype)) * 0.5))

# Feature plots of key markers
key_markers <- c("P2RY12", "TMEM119", "AQP4", "GFAP", "MBP", "MOG",
                 "PDGFRA", "VWF", "GAD1", "NELL2", "CD3D")

key_markers_avail <- key_markers[key_markers %in% rownames(sobj)]
if (length(key_markers_avail) > 0) {
  feature_plots <- FeaturePlot(sobj,
                              features = key_markers_avail,
                              reduction = "umap",
                              pt.size = 0.1,
                              ncol = 3,
                              order = TRUE) &
    scale_color_viridis_c() &
    theme(plot.title = element_text(size = 10))
  
  ggsave(file.path(output_dir, "plot/annotation/feature_plots_key_markers.pdf"),
         plot = feature_plots,
         width = 12,
         height = 4 * ceiling(length(key_markers_avail) / 3))
}

# Violin plots of top markers per cell type
cat("  Generating violin plots...\n")
top_markers_per_celltype <- cluster_markers %>%
  filter(gene %in% available_markers) %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC) %>%
  ungroup() %>%
  select(gene, cluster) %>%
  distinct()

if (nrow(top_markers_per_celltype) > 0) {
  # Limit to 12 markers for readability
  markers_to_plot <- head(unique(top_markers_per_celltype$gene), 12)
  
  violin_plots <- VlnPlot(sobj,
                          features = markers_to_plot,
                          group.by = "celltype",
                          pt.size = 0,
                          ncol = 4,
                          slot = "data") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "plot/annotation/violin_plots_top_markers.pdf"),
         plot = violin_plots,
         width = 16,
         height = 12)
}

# Step 6: Cell type composition analysis
cat("\n6. Analyzing cell type composition...\n")

# Composition by condition
composition_by_condition <- sobj@meta.data %>%
  group_by(condition, celltype) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(percentage = n_cells / sum(n_cells) * 100)

write.csv(composition_by_condition,
          file.path(output_dir, "Documents/annotation/celltype_composition_by_condition.csv"),
          row.names = FALSE)

# Composition by sample
composition_by_sample <- sobj@meta.data %>%
  group_by(sample, condition, celltype) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(percentage = n_cells / sum(n_cells) * 100)

write.csv(composition_by_sample,
          file.path(output_dir, "Documents/annotation/celltype_composition_by_sample.csv"),
          row.names = FALSE)

# Plot composition
comp_plot <- ggplot(composition_by_condition,
                    aes(x = condition, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Condition", y = "Percentage", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "plot/annotation/celltype_composition.pdf"),
       plot = comp_plot, width = 8, height = 6)

# Step 7: Save annotated object
cat("\n7. Saving annotated object...\n")
saveRDS(sobj, file = annotated_object)
cat("  Saved to:", annotated_object, "\n")

# Step 8: Final summary
cat("\n=== GSE188545 Annotation Complete ===\n")
cat("\nSummary:\n")
cat("  Total cells:", ncol(sobj), "\n")
cat("  Cell types identified:", length(unique(sobj$celltype)), "\n")
cat("\nCell type counts:\n")
celltype_counts <- table(sobj$celltype)
for (ct in names(celltype_counts)) {
  cat("  ", ct, ":", celltype_counts[ct], "cells (",
      round(celltype_counts[ct]/ncol(sobj)*100, 1), "%)\n")
}

cat("\nOutput files:\n")
cat("  Annotated object:", annotated_object, "\n")
cat("  Plots:", file.path(output_dir, "plot/annotation"), "\n")
cat("  Annotation documents:", file.path(output_dir, "Documents/annotation"), "\n")

# Save session info
sink(file.path(output_dir, "Documents/annotation/session_info_annotation.txt"))
sessionInfo()
sink()

cat("\nNext steps: Analyze differential expression between conditions.\n")
cat("\nScript completed successfully.\n")