library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

##plot proportion dodge
plot_proportion_dodge <- function(sobj,celltypes,group){

  
  # sobj[[celltypes]]<-as.factor(sobj[[celltypes]])
  
  
  cell_number<-(prop.table(table(sobj@meta.data[[celltypes]],sobj@meta.data[[group]]),margin = 2)*100)%>%
    round(digits = 2)
  cell_number<-as.data.frame(cell_number)
  
  colnames(cell_number)<-c('celltype','Group','percentage')
  # ordered_sample <-c(
  #   "0","1","2","3","4","5","6","7","8","9","10" # 按照实际的分组顺序调整
  # )
  # cell_number$celltype <- factor(cell_number$celltype, levels = ordered_sample)
  
  
  p <- ggplot(cell_number, aes(x = celltype, 
                               y = percentage,
                               fill = Group)) +
    geom_bar(position = position_dodge(0.85), # 微调分组间距
             stat = "identity", 
             width = 0.75, # 调整柱子宽度
             color = "white", # 添加白色边框
             linewidth = 0.3) + # 边框细线
    geom_text(aes(label = percentage), # 添加数值标签
              position = position_dodge(0.75),
              vjust = -0.5, 
              size = 3, 
              color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # 优化坐标轴空白
    labs(x = "Cell Type", 
         y = "Percentage (%)",
         title = "Cell Type Proportion by Group") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "top",
      panel.grid.major.x = element_blank()
    )
  mycol <- c("#4d97cd", "#db6968", "#99cbeb", "#459943",
             "#fdc58f", "#e8c559", "#a3d393", "#f8984e")
  p <- p + scale_fill_manual(values = mycol)
  
  
  
  return(p)
}

plot_proportion_stack <- function(sobj, celltypes, group) {
  
  # 计算每个分组中各细胞类型的百分比
  cell_number <- (prop.table(table(sobj@meta.data[[celltypes]], sobj@meta.data[[group]]), margin = 2) * 100) %>%
    round(digits = 2)
  cell_number <- as.data.frame(cell_number)
  
  colnames(cell_number) <- c('celltype', 'Group', 'percentage')
  
  # 创建堆叠图
  p <- ggplot(cell_number, aes(x = Group,  # x轴改为分组变量
                               y = percentage,
                               fill = celltype)) +  # 填充改为细胞类型
    geom_bar(position = "stack",  # 使用堆叠位置
             stat = "identity", 
             width = 0.75, 
             color = "white", 
             linewidth = 0.3) +
    geom_text(aes(label = percentage),  # 调整文本位置适应堆叠
              position = position_stack(vjust = 0.5),  # 文本居中显示在每个堆叠部分
              size = 3, 
              color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                       limits = c(0, 100)) +  # 确保Y轴从0到100
    labs(x = "Group",  # 调整x轴标签
         y = "Percentage (%)",
         title = "Cell Type Proportion by Group") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "top",
      panel.grid.major.x = element_blank()
    )
  
  # 使用相同的颜色方案
  mycol <- c("#4d97cd", "#db6968", "#99cbeb", "#459943",
             "#fdc58f", "#e8c559", "#a3d393", "#f8984e")
  p <- p + scale_fill_manual(values = mycol)
  
  return(p)
}

plot_stacked_proportions <- function(sobj, sample_col, class_col, mycol1,
                                     x_label = "sample", y_label = "proportion",
                                     plot_title = NULL) {
  # 检查必要的包
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("需要安装ggplot2包才能运行此函数")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("需要安装dplyr包才能运行此函数")
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("需要安装ggpubr包以使用RotatedAxis功能")
  }
  
  # 验证输入列是否存在
  if (!sample_col %in% colnames(sobj@meta.data)) {
    stop(paste("样本列", sample_col, "不存在于meta.data中"))
  }
  if (!class_col %in% colnames(sobj@meta.data)) {
    stop(paste("分类列", class_col, "不存在于meta.data中"))
  }
  
  # 准备数据
  sample_table <- as.data.frame(table(
    sobj@meta.data[[sample_col]],
    sobj@meta.data[[class_col]]
  ))
  names(sample_table) <- c("Samples", "class", "CellNumber")
  mycol1 <- c("#4d97cd", "#db6968", "#99cbeb", "#459943",
             "#fdc58f", "#e8c559", "#a3d393", "#f8984e")
  # 绘制堆叠条形图
  p <- sample_table %>%
    ggplot2::ggplot(ggplot2::aes(x = Samples, weight = CellNumber, fill = class)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::scale_fill_manual(values = mycol1) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      axis.line.x = ggplot2::element_line(colour = "black"),
      axis.line.y = ggplot2::element_line(colour = "black"),
      plot.title = ggplot2::element_text(lineheight = .8, face = "bold", 
                                         hjust = 0.5, size = 16),
      axis.title = ggplot2::element_text(size = 15),
      axis.text = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::labs(x = x_label, y = y_label, title = plot_title) +
    RotatedAxis()
  
  return(p)
}

# celltypes <- sym("subclass") # 替换为你的细胞类型列名

# cell_ratio <- sobj@meta.data %>%
#   # 按样本(orig.ident)和细胞类型(celltypes)分组
#   group_by(orig.ident, !!celltypes) %>%
#   # 统计每个样本中每种细胞类型的细胞数量
#   summarise(cell_count = n(), .groups = "drop_last") %>%
#   # 按样本计算每种细胞类型占该样本总细胞数的比例
#   mutate(ratio = cell_count / sum(cell_count),
#          percentage = ratio * 100) %>%  # 可以转换成百分比形式
#   ungroup()


calculate_ratio <- function(sobj, celltypes,Group) {
  celltypes <- sym(celltypes)  # 将字符串转换为符号
  sample_info <- sobj@meta.data %>%
    select(orig.ident, Group) %>%  # 选择需要保留的样本信息列
    distinct() 
  #calculate ratio
  cell_ratio <- sobj@meta.data %>%
    group_by(orig.ident, !!celltypes) %>%  # 使用{{}}引用变量
    summarise(cell_count = n(), .groups = "drop_last") %>%
    mutate(ratio = cell_count / sum(cell_count),
           percentage =ratio * 100) %>%
    left_join(sample_info, by = "orig.ident") %>%
    ungroup()
  #
  
}

# library(dplyr)
# library(purrr)
# 
# # 按subclass分组，对每个分组执行Wilcoxon检验
# test_results <- cell_ratio %>%
#   group_by(subclass) %>%
#   group_split() %>%  # 按subclass拆分为列表
#   map_dfr(function(df) {
#     # 提取当前subclass名称
#     current_subclass <- unique(df$subclass)
#     
#     # 检查Group数量是否≥2（检验需要至少两组）
#     if (length(unique(df$Group)) < 2) {
#       return(tibble(
#         subclass = current_subclass,
#         statistic = NA,
#         p_value = NA,
#         method = "Wilcoxon rank sum test",
#         note = "Group数量不足，无法进行检验"
#       ))
#     }
#     
#     # 执行Wilcoxon检验（比较不同Group间的percentage）
#     test <- wilcox.test(percentage ~ Group, data = df)
#     
#     # 提取检验结果
#     tibble(
#       subclass = current_subclass,
#       statistic = test$statistic[[1]],  # 检验统计量
#       p_value = test$p.value,           # p值
#       method = test$method,             # 检验方法
#       note = "检验成功"
#     )
#   })
# 
# # 输出结果
# print(test_results, n = Inf)  # n=Inf显示所有行