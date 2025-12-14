
# ============================================================================
# Code to plot all genes from ECM_structure_genes table
# ============================================================================

# Load the ECM_structure_genes data
library(readxl)
library(here)

# Load ECM structure genes from Excel file
ECM_structure_genes <- read_excel(
  here::here("data", "ECM_related_genes", "ECM_structure_genes.xlsx"), 
  sheet = 1
)

# Extract gene symbols from the table
# Assuming the gene names are in a column (adjust column name/index as needed)
# If genes are in the first column, use:
genes_to_plot <- ECM_structure_genes[[1]]

# Or if genes are in a column named "Gene" or "Symbol":
# genes_to_plot <- ECM_structure_genes$Gene
# genes_to_plot <- ECM_structure_genes$Symbol

# Filter genes that exist in the Seurat object
available_genes <- genes_to_plot[genes_to_plot %in% rownames(sobj)]

# Print summary
cat("Total genes in ECM_structure_genes:", length(genes_to_plot), "\n")
cat("Genes available in Seurat object:", length(available_genes), "\n")
cat("Missing genes:", length(genes_to_plot) - length(available_genes), "\n")

# Set cell type identity (adjust "oupSample.cellType" to your actual cell type column)
celltype_column <- "oupSample.cellType"  # Change this to match your metadata

# Create the dotplot
dotplot_ECM <- dotplot_fy(
  sobj = sobj, 
  celltypes = celltype_column, 
  markers = available_genes
)

# Display the plot
print(dotplot_ECM)

# Save the plot
ggsave(
  filename = here::here("task", "20251130-2-ECM_related_genes", "plot", "ECM_structure_genes_dotplot.png"),
  plot = dotplot_ECM,
  width = 12,
  height = 10,
  dpi = 300
)

