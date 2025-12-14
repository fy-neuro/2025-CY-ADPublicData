# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(here)

here("src/getdir.R") %>% source()



# Set working directory relative to script location
script_dir <- get_script_path()
cat("Current script directory:", script_dir, "\n")

# Navigate to project root (assuming script is in task/20251126-1-GSE138852/code/)
project_root <- file.path(script_dir, "..", "..", "..")
setwd(project_root)
cat("Working directory set to:", getwd(), "\n")

# ...existing code...


# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)

# Set working directory to project root
# Load the count matrix from compressed CSV
cat("Loading count matrix from data/GSE138852_counts.csv.gz...\n")
counts <- read.csv("data/GSE138852_counts.csv.gz", 
                   row.names = 1, 
                   header = TRUE, 
                   check.names = FALSE)

# Convert to sparse matrix for memory efficiency
cat("Converting to sparse matrix...\n")
counts_sparse <- as(as.matrix(counts), "sparseMatrix")

# Create Seurat object
cat("Creating Seurat object...\n")
sobj <- CreateSeuratObject(
  counts = counts_sparse,
  project = "GSE138852",
  min.cells = 3,    # Include features expressed in at least 3 cells
  min.features = 200 # Include cells with at least 200 features
)

# Add metadata
sobj$dataset <- "GSE138852"

# Basic information about the dataset
cat("Dataset summary:\n")
cat("Number of cells:", ncol(sobj), "\n")
cat("Number of features:", nrow(sobj), "\n")

# Save the Seurat object
cat("Saving Seurat object...\n")
saveRDS(sobj, file = "task/20251126-1-GSE138852/GSE138852_sobject.rds")

cat("Data loading complete! Seurat object saved.\n")

# Display basic statistics
print(sobj)