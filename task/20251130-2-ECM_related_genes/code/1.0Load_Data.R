library(Seurat)
library(dplyr)
library(Matrix)
library(here)

# Get the script path using different methods
# Method 1: Using rstudioapi (works in RStudio)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getSourceEditorContext()$path
  script_dir <- dirname(script_path)
  cat("Script path (RStudio):", script_path, "\n")
  cat("Script directory:", script_dir, "\n")
}

# Method 2: Using here package for project root
project_root <- here::here()
cat("Project root (here):", project_root, "\n")

file_path <- here::here("data", "ECM_related_genes")
documents <- list.files(file_path, pattern = "*.xlsx", full.names = TRUE)
doc_names <- list.files(file_path, pattern = "*.xlsx", full.names = FALSE) %>% gsub(".xlsx$", "", .)
documents
doc_names


for (i in 1:length(documents)) {
    cat("Loading ECM related genes from:", documents[i], "\n")
    gene_data <- readxl::read_excel(documents[i], sheet = 1)
    assign(doc_names[i], gene_data)
}

sobj_path <- here::here("data", "GSE138852_sobj.rds")
cat("Loading Seurat object from:", sobj_path, "\n")
sobj <- readRDS(sobj_path)
sobj
