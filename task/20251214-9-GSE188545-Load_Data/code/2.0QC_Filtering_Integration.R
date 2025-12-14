# QC Filtering and Integration for GSE188545
# Performs QC filtering, normalization, Harmony integration, and clustering

# Load required libraries
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(here)

# Set working directory to project root
project_root <- here::here()
setwd(project_root)
cat("Working directory set to:", getwd(), "\n")

# Configuration
input_dir <- "task/20251214-9-GSE188545-Load_Data/"
output_dir <- input_dir  # Same directory for outputs
raw_object <- file.path(input_dir, "GSE188545_sobj_raw.rds")
filtered_object <- file.path(output_dir, "GSE188545_sobj_filtered.rds")
integrated_object <- file.path(output_dir, "GSE188545_sobj_integrated.rds")

cat("=== GSE188545 QC Filtering and Integration ===\n")
cat("Input directory:", input_dir, "\n")
cat("Raw object:", raw_object, "\n")

# Create output directories
dir.create(file.path(output_dir, "plot"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "Documents"), showWarnings = FALSE, recursive = TRUE)

# Step 1: Load raw Seurat object
cat("\n1. Loading raw Seurat object...\n")
if (!file.exists(raw_object)) {
  stop("Raw object not found. Run 1.0Load_Data.R first.")
}

sobj <- readRDS(raw_object)
cat("  Loaded object:", sobj@project.name, "\n")
cat("  Cells:", ncol(sobj), "\n")
cat("  Features:", nrow(sobj), "\n")
cat("  Samples:", length(unique(sobj$sample)), "\n")

# Fix condition metadata (original loading script may have misclassified)
cat("  Correcting condition metadata...\n")
sobj$condition <- ifelse(grepl("_AD", sobj$sample), "AD", 
                         ifelse(grepl("_HC", sobj$sample), "HC", "Unknown"))
cat("    AD cells:", sum(sobj$condition == "AD"), "\n")
cat("    HC cells:", sum(sobj$condition == "HC"), "\n")

# Step 2: QC filtering
cat("\n2. Applying QC filters...\n")

# Calculate additional QC metrics if not present
if (!"percent.mt" %in% colnames(sobj@meta.data)) {
  sobj$percent.mt <- PercentageFeatureSet(sobj, pattern = "^MT-")
}
if (!"percent.ribo" %in% colnames(sobj@meta.data)) {
  sobj$percent.ribo <- PercentageFeatureSet(sobj, pattern = "^RP[SL]")
}

# Plot pre-filtering QC
cat("  Generating pre-filtering QC plots...\n")
pre_filter_dir <- file.path(output_dir, "plot/pre_filter")
dir.create(pre_filter_dir, showWarnings = FALSE, recursive = TRUE)

# Violin plots pre-filtering
pre_violin <- VlnPlot(sobj,
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      group.by = "condition",
                      pt.size = 0,
                      ncol = 3) +
  ggtitle("Pre-filtering QC metrics")

ggsave(file.path(pre_filter_dir, "pre_filter_violin.pdf"),
       plot = pre_violin, width = 12, height = 5)

# Scatter plots pre-filtering
pre_scatter <- FeatureScatter(sobj,
                             feature1 = "nCount_RNA",
                             feature2 = "percent.mt",
                             group.by = "condition") +
  ggtitle("Pre-filtering: nCount_RNA vs %MT")

ggsave(file.path(pre_filter_dir, "pre_filter_scatter.pdf"),
       plot = pre_scatter, width = 8, height = 6)

# Set QC thresholds based on data distribution
cat("  Calculating QC thresholds...\n")

# Analyze distributions to set thresholds
nfeature_stats <- quantile(sobj$nFeature_RNA, probs = c(0.01, 0.99))
ncount_stats <- quantile(sobj$nCount_RNA, probs = c(0.01, 0.99))
mt_stats <- quantile(sobj$percent.mt, probs = c(0.95, 0.99))

cat("  nFeature_RNA distribution (1-99%):", round(nfeature_stats), "\n")
cat("  nCount_RNA distribution (1-99%):", round(ncount_stats), "\n")
cat("  %MT distribution (95-99%):", round(mt_stats, 2), "\n")

# Define QC thresholds (adjust based on data)
# Typical thresholds for human brain single-nucleus RNA-seq
nfeature_min <- 500    # Minimum genes per cell
nfeature_max <- 5000   # Maximum genes per cell (avoid doublets)
ncount_min <- 1000     # Minimum UMIs per cell  
ncount_max <- 30000    # Maximum UMIs per cell
mt_threshold <- 10     # Maximum mitochondrial percentage

cat("  Applying thresholds:\n")
cat("    nFeature_RNA:", nfeature_min, "-", nfeature_max, "\n")
cat("    nCount_RNA:", ncount_min, "-", ncount_max, "\n")
cat("    %MT <", mt_threshold, "\n")

# Apply filtering
cells_pre <- ncol(sobj)
sobj_filtered <- subset(sobj,
                        subset = nFeature_RNA > nfeature_min &
                                 nFeature_RNA < nfeature_max &
                                 nCount_RNA > ncount_min &
                                 nCount_RNA < ncount_max &
                                 percent.mt < mt_threshold)

cells_post <- ncol(sobj_filtered)
cells_removed <- cells_pre - cells_post
cat("  Filtering results:\n")
cat("    Cells before:", cells_pre, "\n")
cat("    Cells after:", cells_post, "\n")
cat("    Cells removed:", cells_removed, "(", 
    round(cells_removed/cells_pre * 100, 1), "%)\n")

# Save filtered object
cat("  Saving filtered object...\n")
saveRDS(sobj_filtered, file = filtered_object)
cat("    Saved to:", filtered_object, "\n")

# Plot post-filtering QC
cat("  Generating post-filtering QC plots...\n")
post_filter_dir <- file.path(output_dir, "plot/post_filter")
dir.create(post_filter_dir, showWarnings = FALSE, recursive = TRUE)

post_violin <- VlnPlot(sobj_filtered,
                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       group.by = "condition",
                       pt.size = 0,
                       ncol = 3) +
  ggtitle("Post-filtering QC metrics")

ggsave(file.path(post_filter_dir, "post_filter_violin.pdf"),
       plot = post_violin, width = 12, height = 5)

# Step 3: Normalization and variable feature selection
cat("\n3. Normalization and variable feature selection...\n")

# Use traditional normalization to save memory
cat("  Normalizing data...\n")
sobj_filtered <- NormalizeData(sobj_filtered, 
                             normalization.method = "LogNormalize", 
                             scale.factor = 10000)

cat("  Finding variable features...\n")
sobj_filtered <- FindVariableFeatures(sobj_filtered, 
                                      selection.method = "vst", 
                                      nfeatures = 2000)

cat("  Scaling data...\n")
sobj_filtered <- ScaleData(sobj_filtered, 
                          vars.to.regress = "percent.mt")

cat("    Variable features:", length(VariableFeatures(sobj_filtered)), "\n")

# Step 4: Dimensionality reduction and Harmony integration
cat("\n4. Dimensionality reduction and Harmony integration...\n")

# Run PCA
cat("  Running PCA...\n")
sobj_filtered <- RunPCA(sobj_filtered,
                        features = VariableFeatures(sobj_filtered),
                        npcs = 50,
                        verbose = FALSE)

# Determine significant PCs (Elbow plot)
cat("  Generating Elbow plot...\n")
elbow_plot <- ElbowPlot(sobj_filtered, ndims = 50) +
  ggtitle("Elbow Plot for PCA")

ggsave(file.path(output_dir, "plot/PCA_elbow_plot.pdf"),
       plot = elbow_plot, width = 8, height = 6)

# Run Harmony for batch correction
cat("  Running Harmony integration...\n")
sobj_integrated <- RunHarmony(sobj_filtered,
                              group.by.vars = "sample",
                              reduction.use = "pca",
                              project.dim = FALSE,
                              verbose = FALSE)

cat("  Harmony completed.\n")

# Step 5: Clustering
cat("\n5. Clustering cells...\n")

# Find neighbors using Harmony dimensions
cat("  Finding neighbors...\n")
sobj_integrated <- FindNeighbors(sobj_integrated,
                                 reduction = "harmony",
                                 dims = 1:30)

# Find clusters at multiple resolutions
cat("  Finding clusters at multiple resolutions...\n")
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
sobj_integrated <- FindClusters(sobj_integrated,
                                resolution = resolutions)

cat("  Cluster resolutions applied:", resolutions, "\n")

# Step 6: UMAP visualization
cat("\n6. Running UMAP...\n")
sobj_integrated <- RunUMAP(sobj_integrated,
                           reduction = "harmony",
                           dims = 1:30,
                           min.dist = 0.3,
                           n.neighbors = 30,
                           verbose = FALSE)

cat("  UMAP completed.\n")

# Step 7: Save integrated object
cat("\n7. Saving integrated object...\n")
saveRDS(sobj_integrated, file = integrated_object)
cat("  Saved to:", integrated_object, "\n")

# Step 8: Generate diagnostic plots
cat("\n8. Generating diagnostic plots...\n")

# UMAP by condition
umap_condition <- DimPlot(sobj_integrated,
                          reduction = "umap",
                          group.by = "condition",
                          pt.size = 0.1) +
  ggtitle("UMAP by Condition (AD vs HC)")

# UMAP by sample
umap_sample <- DimPlot(sobj_integrated,
                       reduction = "umap",
                       group.by = "sample",
                       pt.size = 0.1) +
  theme(legend.text = element_text(size = 8)) +
  ggtitle("UMAP by Sample")

# UMAP by cluster (default resolution 0.5)
umap_cluster <- DimPlot(sobj_integrated,
                        reduction = "umap",
                        group.by = "RNA_snn_res.0.5",
                        label = TRUE,
                        pt.size = 0.1) +
  ggtitle("UMAP by Cluster (resolution 0.5)")

# Combine plots
diagnostic_plots <- (umap_condition + umap_sample) / umap_cluster +
  plot_layout(heights = c(1, 1.5))

ggsave(file.path(output_dir, "plot/integration_diagnostic_plots.pdf"),
       plot = diagnostic_plots, width = 14, height = 12)

# QC metrics on UMAP
feature_plots <- FeaturePlot(sobj_integrated,
                            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                            reduction = "umap",
                            pt.size = 0.1,
                            ncol = 3) &
  scale_color_viridis_c()

ggsave(file.path(output_dir, "plot/QC_feature_plots.pdf"),
       plot = feature_plots, width = 15, height = 5)

# Step 9: Generate summary statistics
cat("\n9. Generating summary statistics...\n")

# Cell counts by sample and condition
cell_counts <- sobj_integrated@meta.data %>%
  group_by(sample, condition) %>%
  summarise(n_cells = n(),
            .groups = "drop")

# Cluster composition
cluster_composition <- sobj_integrated@meta.data %>%
  group_by(RNA_snn_res.0.5, condition) %>%
  summarise(n_cells = n(),
            .groups = "drop") %>%
  group_by(RNA_snn_res.0.5) %>%
  mutate(percentage = n_cells / sum(n_cells) * 100)

# Save statistics
write.csv(cell_counts, 
          file.path(output_dir, "Documents/cell_counts_by_sample.csv"),
          row.names = FALSE)

write.csv(cluster_composition,
          file.path(output_dir, "Documents/cluster_composition.csv"),
          row.names = FALSE)

# Step 10: Final summary
cat("\n=== GSE188545 Integration Complete ===\n")
cat("\nSummary:\n")
cat("  Cells after filtering:", ncol(sobj_integrated), "\n")
cat("  Samples:", length(unique(sobj_integrated$sample)), "\n")
cat("  AD cells:", sum(sobj_integrated$condition == "AD"), "\n")
cat("  HC cells:", sum(sobj_integrated$condition == "HC"), "\n")
cat("  Clusters (res 0.5):", length(unique(sobj_integrated$RNA_snn_res.0.5)), "\n")
cat("\nOutput files:\n")
cat("  Filtered object:", filtered_object, "\n")
cat("  Integrated object:", integrated_object, "\n")
cat("  Plots directory:", file.path(output_dir, "plot"), "\n")
cat("  Documents directory:", file.path(output_dir, "Documents"), "\n")

# Save session info
sink(file.path(output_dir, "Documents/session_info_filtering_integration.txt"))
sessionInfo()
sink()

cat("\nNext steps: Run cell type annotation (3.0Annotation.R)\n")
cat("\nScript completed successfully.\n")