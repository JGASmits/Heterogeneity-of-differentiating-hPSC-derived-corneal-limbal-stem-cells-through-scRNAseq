# Used conda environment: JA_R_seuratdisk
library(Seurat)
library(SeuratDisk)
seu_ob = readRDS("differentiation_object.rds")
SaveH5Seurat(seu_ob, filename = "converted.h5Seurat")
Convert("converted.h5Seurat", dest = "h5ad")
sessionInfo()