# 细胞类型标记基因向量（小鼠）
celltype_markers <- list(
  Microglia = c("P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Hexb"),
  Macrophages = c("Adgre1", "Mrc1", "Cd68", "Itgam", "Cd14"),
  `T cells` = c("Cd3e", "Cd4", "Cd8a", "Tcrb", "Lck"),
  Oligodendrocytes = c("Mbp", "Plp1", "Olig1", "Olig2", "Cnp"),
  `Endthelial cells` = c("Pecam1", "Cd34", "Vwf", "Esam", "Cldn5"),
  `Proliferating cells` = c("Mki67", "Pcna", "Top2a", "Ccna2", "Cdk1"),
  Neutrophils = c("S100a8", "S100a9", "Ly6g", "Elane", "Mpo"),
  Neurons = c("Snap25", "Syt1", "Map2", "NeuN", "Tubb3", "Elavl3"),
  `Epithelial Cells` = c("Epcam", "Cdh1", "Krt8", "Krt18", "Krt19")
)

# 整合所有细胞类型标记基因的单个向量（小鼠）
all_celltype_markers <- c(
  # Microglia 标记基因
  "P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Hexb",
  # Macrophages 标记基因
  "Adgre1", "Mrc1", "Cd68", "Itgam", "Cd14",
  # T cells 标记基因
  "Cd3e", "Cd4", "Cd8a", "Tcrb", "Lck",
  # Oligodendrocytes 标记基因
  "Mbp", "Plp1", "Olig1", "Olig2", "Cnp",
  # Endthelial cells 标记基因
  "Pecam1", "Cd34", "Vwf", "Esam", "Cldn5",
  # Proliferating cells 标记基因
  # "Mki67", "Pcna", "Top2a", "Ccna2", "Cdk1", "Pola1",
  # Neutrophils 标记基因
  "S100a8", "S100a9", "Ly6g", "Elane", "Mpo",
  # Neurons 标记基因
  "Snap25", "Syt1", "Map2", "NeuN", "Tubb3", "Elavl3",
  # Epithelial Cells 标记基因
  "Epcam", "Cdh1", "Krt8", "Krt18", "Krt19"
)

Markers_CNS<-unique(c(
"TH","SLC6A3"
#
,'RALYL','LDB2','NELL2'
#GlutamatergicNeuron	
,'GAD1','GAD2','SLC6A1'
  # GABAergic Neuron
  ,'PDGFRA','VCAN','OLIG1'
  # OPC
  ,'MBP','MOBP','MOG',"OPALIN"
  # OLIGO
  ,'AQP4','GFAP','FGFR3',"NHSL1","SLC25A18"
  # ASTRO
  ,  'CD2','THEMIS','CD3D'
   # T CELL
  ,'ITGAM','CD74','CX3CR1','P2RY12',"C3","CSF1R"
  # MICRO
  ,'DCN','FLT1',"LEF1","VWF"
  # VC
))


# Markers_CNS<-unique(c(
#   # "TH","SLC6A3"
#   #
#   # ,
#   'RALYL','LDB2','NELL2'
#   #GlutamatergicNeuron	
#   ,'GAD1','GAD2','SLC6A1'
#   # GABAergic Neuron
#   ,'PDGFRA','VCAN','OLIG1'
#   # OPC
#   ,'MBP','MOBP','MOG',"OPALIN"
#   # OLIGO
#   ,'AQP4','GFAP','FGFR3',"NHSL1","SLC25A18"
#   # ASTRO
#   ,  'CD2','THEMIS','CD3D'
#   # T CELL
#   ,'ITGAM','CD74','CX3CR1','P2RY12',"C3","CSF1R"
#   # MICRO
#   ,'DCN','FLT1',"LEF1","VWF"
#   # VC
# )) %>% tolower() %>%
#   stringr::str_to_title()




Markers_all<-c(
  'DCN','FLT1',"LEF1","VWF"
  # VC
  ,'PDGFRA','VCAN','OLIG1'
  # OPC
  ,'MBP','MOBP','MOG',"OPALIN"
  # OLIGO
  ,  'CD2','THEMIS','CD3D'
  # T CELL
  ,'ITGAM','CD74','CX3CR1','P2RY12',"C3","CSF1R"
  # MICRO 
  ,'AQP4','GFAP','FGFR3',"NHSL1","SLC25A18"
  # ASTRO
  ,'RALYL','LDB2','NELL2'
  #Glutamatergic Neuron	
  ,'GAD1','GRIP1','GAD2','SLC6A1'
  #GABAergic Neuron	
  
)