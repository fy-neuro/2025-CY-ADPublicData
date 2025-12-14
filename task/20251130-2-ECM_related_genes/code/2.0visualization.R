library(cols4all)
library(Seurat)
library(ggplot2)

##绘制气泡图
dotplot_fy <- function(sobj, celltypes, markers){
  Idents(sobj) <- celltypes
  Dotplot_re0.5 <- DotPlot(sobj, features = markers, dot.scale = 7,
                           cluster.idents = F) +
    coord_flip() +
    scale_colour_binned_c4a_div('kovesi.bu_wh_rd2') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(
      axis.title = element_text(size = 15),   # 坐标轴标题字体大小
      axis.text = element_text(size = 15),    # 刻度标签字体大小
      legend.text = element_text(size = 15),   # 图注标签（分类名称）的字体大小
      legend.title = element_text(size = 15)    # 图注标题（例如"Celltype"）的字体大小
    ) + labs(x = "Gene Symbols", y = "Class") 
  
  return(Dotplot_re0.5)
}


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
celltype_column <- "oupSample.cellType_batchCond"  # Change this to match your metadata
cell_levels <- unique(sobj$oupSample.cellType_batchCond)
cell_levels <- cell_levels %>% sort()
cell_levels
cell_cell_levels <- c("astro_ct","astro_AD")

sobj$oupSample.cellType_batchCond <- factor(
  sobj$oupSample.cellType_batchCond, 
  levels = cell_levels
)

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
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = here::here("task", "20251130-2-ECM_related_genes", "plot", "ECM_structure_genes_dotplot.pdf"),
  plot = dotplot_ECM,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = here::here("task", "20251130-2-ECM_related_genes", "plot", "ECM_structure_genes_dotplot.eps"),
  plot = dotplot_ECM,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)
