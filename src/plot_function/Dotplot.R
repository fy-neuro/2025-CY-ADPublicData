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

