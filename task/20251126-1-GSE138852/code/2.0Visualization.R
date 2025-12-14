library(ggplot2)
VlnPlot(sobj, features = "SDC4",split.by  = "oupSample.batchCond") +
  theme_minimal()
ggsave("../plot/SDC4_violin_plot.png")

