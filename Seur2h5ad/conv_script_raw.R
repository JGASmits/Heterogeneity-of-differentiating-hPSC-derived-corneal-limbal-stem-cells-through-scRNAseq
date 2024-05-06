# Used conda environment: JA_R_seuratdisk
library(Seurat)
library(SeuratDisk)
seu_ob = readRDS("seu_unfiltered.RDS")
SaveH5Seurat(seu_ob, filename = "converted_raw.h5Seurat")
Convert("converted_raw.h5Seurat", dest = "h5ad")
sessionInfo()