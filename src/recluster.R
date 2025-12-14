# Load required libraries
library(Seurat)  # For single-cell RNA-seq analysis
library(harmony) # For batch effect correction



# Function to perform reclustering on a Seurat object
# Parameters:
#   sobj: A Seurat object containing single-cell RNA-seq data
# Returns:
#   A Seurat object with updated clustering and dimensionality reduction
recluster <- function(sobj ,
                      reduction_use = "harmony",
                      integrate_by = "orig.ident"
                      ) {
  
  # Specify the dimensionality reduction method to use for clustering
  reduction_use <- "pca"  # Options: "PCA" or "harmony"
  
  
  # Step 1: Normalize the gene expression data
  # Uses log-normalization by default
  sobj <- SCTransform(sobj, verbose = FALSE)
  
  
  # Step 4: Perform Principal Component Analysis (PCA)
  # Uses the variable features identified in step 2
  sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))
  
  # Step 5: Run Harmony for batch effect correction
  # Integrates data across different samples/batches using the "orig.ident" metadata
  sobj <- RunHarmony(sobj, group.by.vars = "orig.ident", reduction.use = "pca")
  
  # Step 6: Construct a k-nearest neighbor (KNN) graph
  # Uses the first 20 Harmony dimensions
  sobj <- FindNeighbors(sobj, reduction = reduction_use, dims = 1:20)
  
  # Step 7: Identify cell clusters using multiple resolutions
  # Tests resolutions from 0.1 to 1.0 in increments of 0.1
  ##sobj <- FindClusters(sobj, resolution = seq(0.1, 1, 0.1))
  
  # Step 8: Generate UMAP dimensionality reduction for visualization
  # Uses the first 20 Harmony dimensions
  sobj <- RunUMAP(sobj, dims = 1:20, reduction = reduction_use)
  
  # Return the processed Seurat object
  return(sobj)
}