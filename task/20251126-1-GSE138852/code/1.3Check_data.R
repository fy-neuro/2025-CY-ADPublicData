# =============================================================================
# Script: 1.3Check_data.R
# Purpose: Check data structure, update metadata, and create QC plots
# Author: [Your Name]
# Date: 2025-11-26
# =============================================================================

# Check metadata structure before and after merging
# Compare the column names between original Seurat metadata and merged metadata
names(sobj@meta.data)  # Original Seurat metadata columns
names(Meta.X)          # Enhanced metadata with covariates

# Update Seurat object with merged metadata
# This replaces the basic Seurat metadata with our enriched version
sobj@meta.data <- Meta.X

# Verify the number of metadata variables after merging
# This helps confirm all expected variables are present
Meta.X %>% length() %>% cat("Total metadata columns:", ., "\n")

# Explore sample-level metadata structure
# Select key identifying columns to understand experimental design
Meta.X %>% 
    select(1,4:9) %>%    # Select first column and columns 4-9
    unique()             # Show unique combinations of these variables

# Check unique values in key experimental variables
# These help understand the experimental conditions and cell types
unique(Meta.X$oupSample.batchCond)        # Batch and condition combinations
unique(Meta.X$oupSample.cellType)        # Cell type annotations
unique(Meta.X$oupSample.cellType_batchCond)  # Cell type by batch/condition
unique(Meta.X$oupSample.subclustCond)    # Subcluster by condition

# Set cell identity for downstream analysis
# This determines how cells are grouped in plots and analysis
Idents(sobj) <- "oupSample.cellType"

# =============================================================================
# File path management for saving outputs
# =============================================================================

# Check current working directory
getwd()

# Get the directory where this script is located
# This ensures plots are saved relative to the script location
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getSourceEditorContext()$path  # Full script path
  script_dir <- dirname(script_path)                        # Directory only
  cat("Script path (RStudio):", script_path, "\n")
  cat("Script directory:", script_dir, "\n")
}

# Set working directory to script location
# This ensures relative paths work correctly for output files
setwd(script_dir)
getwd()  # Confirm new working directory

# =============================================================================
# Data visualization and quality control
# =============================================================================

# Load plotting library
library(ggplot2)

# Create violin plot for SDC4 gene expression
# Split by batch/condition to compare expression across experimental groups
plot_SDC4 <- VlnPlot(sobj, 
                     features = "SDC4",                    # Gene of interest
                     split.by = "oupSample.batchCond") +   # Group by batch/condition
  theme_minimal()  # Clean plot appearance

# Display the plot
print(plot_SDC4)

# Save the plot to the plots directory
# Using relative path from script location
ggsave("../plot/SDC4_violin_plot.png", 
       plot = plot_SDC4,
       width = 10, 
       height = 6, 
       dpi = 300)

cat("Plot saved to: ../plot/SDC4_violin_plot.png\n")

