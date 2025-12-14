
library(ggrepel)
library(ggvolcano)
plot_volcano <- function(data){
  p1 <- ggplot(data, aes(x = avg_log2FC, y=neg_log10_padj, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=1) +  #点的透明度、大小
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(-2, 3)) +  #调整点的颜色和x轴的取值范围
  labs(title = "volcano plot ",x="log2FoldChange", y="-log10p_val_adj") +  #x、y轴标签
  ggtitle("volcano plot ") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank()
  ) +
  theme(
    axis.title = element_text(size = 18),   # 坐标轴标题字体大小
    axis.text.x = element_text(size = 18),    # 刻度标签字体大小
    axis.text.y = element_text(size = 18), 
    legend.text = element_text(size = 18),   # 图注标签（分类名称）的字体大小
    legend.title = element_text(size = 18) ,   # 图注标题（例如"Celltype"）的字体大小
    title = element_text(size = 18)
  )+ geom_text_repel(
    data = subset(data, data$neg_log10_padj > 50 ),# 可以设置跟上面不同的阈值，用数值替换即可
    aes(label = gene), size = 4.5,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE ,color = "black")
  
  return(p1)
  }


#出图
# p1
# ggsave("D:/2024-2025/brain/GO/PV-volcano.png",device = "png",width=3000,height = 2000,dpi = 300,units = "px" )

