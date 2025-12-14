library(here)
library(dplyr)
library(ggplot2)
here("src", "recluster.R") %>% source()

sobj <- recluster(sobj)


sobj <- RunUMAP(sobj, dims = 1:10, reduction = "pca")
DimPlot(sobj)
saveRDS(sobj, here::here("data", "GSE138852_sobj_reclustered.rds"))
