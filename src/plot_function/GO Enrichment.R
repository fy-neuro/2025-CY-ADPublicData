perform_GO_analysis <- function(genelist,celltype) {

  # 加载所需包
  library(Seurat)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  
  # 创建输出目录（如果不存在）
  dir.create("../Documents", showWarnings = FALSE, recursive = TRUE)
  dir.create("../plot", showWarnings = FALSE, recursive = TRUE)
  
  # 寻找所有标记基因
  # sobj.markers <- FindAllMarkers(
  #   object = sobj, 
  #   only.pos = TRUE, 
  #   group.by = celltypes,
  #   min.pct = 0.25, 
  #   thresh.use = 0.25 
  # )
  
  # 筛选特定细胞类型的标记基因
  # sobj.markers_celltype <- sobj.markers[sobj.markers$cluster == celltype, ]
  
  # # 过滤掉包含"ENSG"的基因（可能是Ensembl ID）
  # sobj.markers_select <- sobj.markers_celltype %>%
  #   filter(!grepl("ENSG", gene))
  
  # 进一步筛选显著的标记基因
  # sobj.markers_filtered <- sobj.markers_celltype %>% 
  #   filter(p_val_adj < 0.01)
  
  # 检查是否有足够的基因进行分析
  # if (nrow(sobj.markers_filtered) == 0) {
  #   stop(paste("没有找到足够的标记基因用于", celltype, "的GO分析"))
  # }
  
  # 将基因名转换为ENTREZID（小鼠数据库）
  gene_entrez <- as.character(
    na.omit(
      AnnotationDbi::select(
        org.Mm.eg.db,
        keys = genelist,
        columns = 'ENTREZID',
        keytype = 'SYMBOL'
      )[, 2]
    )
  )
  
  # 检查ID转换是否成功
  if (length(gene_entrez) == 0) {
    stop(paste("基因ID转换失败，无法找到", celltype, "的有效ENTREZID"))
  }
  
  # 进行GO富集分析（注意：使用小鼠数据库org.Mm.eg.db）
  go_result <- enrichGO(
    gene = gene_entrez,
    OrgDb = org.Mm.eg.db,  
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # 保存GO分析结果表格
  go_df <- data.frame(go_result)
  result_filename <- file.path("../Documents", paste0(celltype, "_GO.csv"))
  write.table(go_df, result_filename, row.names = FALSE, col.names = TRUE, sep = ",")
  
  # 可视化GO分析结果
  p <- dotplot(go_result, label_format = 20) +
    theme(
      axis.title = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18)
    ) +
    labs(title = paste("GO analysis of", celltype))
  
  # 保存可视化结果
  plot_filename <- file.path("../plot", paste0("GO_", celltype, ".png"))
  ggsave(plot_filename, plot = p, device = "png", 
         width = 3000, height = 3500, dpi = 300, units = "px")
  
  cat(paste("GO分析结果已保存至:", result_filename, "\n"))
  # 返回GO分析结果
  # return(go_result)
}