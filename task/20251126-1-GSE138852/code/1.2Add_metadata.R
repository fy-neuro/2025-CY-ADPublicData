# Load required libraries for single-cell RNA-seq analysis
library(Seurat)      # Main package for single-cell analysis
library(dplyr)       # Data manipulation and pipe operations
library(Matrix)      # Sparse matrix operations
library(here)        # Project-relative file paths
library(tibble)      # Enhanced data frames

# Extract metadata from Seurat object
# This contains basic information about each cell (nFeature_RNA, nCount_RNA, etc.)
Meta.X <- sobj@meta.data

# Convert row names (cell barcodes) to a column for easier joining
# This is necessary because we need a common key to merge with covariates
Meta.X <- Meta.X %>%
  rownames_to_column(var = "barcode")

# Preview the barcode format to ensure correct structure
Meta.X$barcode %>% head()

# Check the structure of covariates data
# This helps us understand what additional metadata we're adding
str(covariates)

# Convert covariates row names to a column for joining
# Both datasets now have a 'barcode' column as the key
covariates <- covariates %>%
  rownames_to_column(var = "barcode")

# Merge Seurat metadata with external covariates
# left_join ensures we keep all cells from Meta.X and add matching covariates
Meta.X <- Meta.X %>%
    left_join(covariates,
              by = "barcode")

# Check the column names after merging
# This shows all available metadata variables for downstream analysis
names(Meta.X)

# Convert barcode back to row names for Seurat compatibility
# Seurat expects metadata to have cell barcodes as row names
Meta.X <- Meta.X %>%
  column_to_rownames(var = "barcode")

# TODO: Add the merged metadata back to the Seurat object
# sobj@meta.data <- Meta.X
