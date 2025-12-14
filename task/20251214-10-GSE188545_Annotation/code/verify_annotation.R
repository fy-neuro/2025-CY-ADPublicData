# Verify annotated object
library(Seurat)
library(dplyr)

cat("Verifying annotated object...\n")
annotated_file <- "GSE188545_sobj_annotated_fast.rds"
if (!file.exists(annotated_file)) {
  stop("Annotated object not found")
}
sobj <- readRDS(annotated_file)
cat("Object loaded:", sobj@project.name, "\n")
cat("Cells:", ncol(sobj), "\n")
cat("Metadata columns:", paste(colnames(sobj@meta.data), collapse = ", "), "\n")
if ("celltype" %in% colnames(sobj@meta.data)) {
  cat("Cell type distribution:\n")
  print(table(sobj$celltype))
} else {
  cat("ERROR: celltype column missing\n")
}
cat("\nUMAP reduction available:", "umap" %in% names(sobj@reductions), "\n")
cat("Annotation mapping from CSV:\n")
mapping <- read.csv("Documents/annotation_fast/cluster_annotation_mapping.csv")
print(mapping)
cat("\nVerification complete.\n")