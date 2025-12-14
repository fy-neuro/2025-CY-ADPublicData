# 加载所需库
library(Seurat)

# 假设已有处理好的Seurat对象sobj
# 若没有，可使用示例数据
# sobj <- pbmc_small

# 定义要展示的marker基因


# 使用Seurat原生VlnPlot函数绘制，合并为一张图
VlnPlot(
  object = sobj,
  features = marker_genes,  # 指定多个基因
  pt.size = 0.1,            # 点的大小
  group.by = "ident",       # 按细胞类型分组
  ncol = 1,                 # 列数，设为1表示所有基因垂直排列
  # nrow = 2,               # 行数，可根据基因数量调整
  combine = TRUE,           # 合并为一张图（默认就是TRUE）
  fill.by = "ident",        # 按细胞类型填充颜色
  cols = NULL,              # 使用默认颜色，也可自定义如c("red", "blue", ...)
  split.by = NULL           # 不拆分
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5),  # 标题居中
    legend.position = "right"
  ) +
  ggtitle("Marker Genes Expression Across Cell Types")  # 添加总标题
