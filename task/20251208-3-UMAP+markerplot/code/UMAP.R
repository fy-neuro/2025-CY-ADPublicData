library(Seurat)
library(ggplot2)
library(dplyr)
library(here)

here("src/getdir.R") %>% source()
setwd(script_dir)


here("data/GSE138852_sobj_reclustered.rds") %>% readRDS() -> sobj
##clean data
sobj <- subset(sobj,subset = (oupSample.cellType != "doublet" & oupSample.cellType != "unID"))
##set metadata
sobj$oupSample.cellType <- sobj$oupSample.cellType %>% as.character() %>% as.factor()
levels(sobj$oupSample.cellType)
sobj$oupSample.cellType <- factor(sobj$oupSample.cellType, levels = c(
  "neuron","astro","mg","oligo","OPC","endo"
))
Idents(sobj) <- "oupSample.cellType"

##UMAP plot
library(cols4all)
c4a_gui()
DimPlot(sobj)
mycol <- c4a("brewer.pastel1")
mycol1 <- c4a("brewer.paired")
mycol
DimPlot(sobj, reduction = "umap"
        # , group.by = "seurat_clusters"
        ,label = TRUE
        ,label.size = 5
        , cols = mycol1) 
# + ggtitle("UMAP colored by clusters")
ggsave("../plot/UMAP.pdf",width = 6,height = 5)

here("src/1stAnnotation.R") %>% source()
here("src/plot_function/Dotplot.R") %>% source()


dotplot_fy(sobj = sobj, celltypes = "oupSample.cellType",
           markers = Markers_all) -> p1
ggsave("../plot/Marker_genes_dotplot.pdf",plot = p1,width = 8,height = 6)

saveRDS(sobj,file = here("data/GSE138852_sobj_cleaned.rds"))
