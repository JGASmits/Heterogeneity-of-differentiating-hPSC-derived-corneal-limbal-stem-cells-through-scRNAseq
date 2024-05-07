# First time install:
#remotes::install_github('satijalab/seurat-wrappers') and choose "3"
#remotes::install_github('eclarke/ggbeeswarm') and choose "3"

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggbeeswarm)
seur_obj <- readRDS("/mnt/s/Radboud/data/scRNAseq_ipsc/correct/differentiation_object.rds")

# Rename and reorder idents according to the latest draft
seur_obj <- RenameIdents(object = seur_obj,  'iLSCs-diff' = 'late-epi', 'iLSCs' = 'iLSCs', 'meso-like-1' = 'meso-like-3', 'mesodermal' = 'meso-like-1')
ord_vec = c('PSC', 'PSC-like', 'meso-like-3', 'meso-like-2', 'meso-like-1', 'early-epi', 'iLSCs', 'late-epi')

seur_obj@active.ident <- factor(seur_obj@active.ident, levels = ord_vec)

seur_obj@meta.data$states <- seur_obj@active.ident

seur_obj <- SetIdent(seur_obj,value=seur_obj@meta.data$cell_type)

# Rename the cell lines
seur_obj <- RenameIdents(object = seur_obj, 'ESC' = 'hESC', 'B2HT' = 'hiPSC1', '003bC' = 'hiPSC2')
ord_vec = c('hESC', 'hiPSC1', 'hiPSC2')

seur_obj@active.ident <- factor(seur_obj@active.ident, levels = ord_vec)

seur_obj@meta.data$cell_type <- seur_obj@active.ident

# Go back to cell state ident
seur_obj <- SetIdent(seur_obj,value=seur_obj@meta.data$states)

# Full object
cds <- as.cell_data_set(seur_obj)
cds <- cluster_cells(cds, reduction_method = "UMAP") #Monocle3 Pipeline
cds <- estimate_size_factors(cds)
cds <- learn_graph(cds, use_partition = F, close_loop = T)
root_group <- colnames(cds)[pData(cds)$ident == "PSC"]
cds <- order_cells(cds, root_cells = root_group)

# Generate colors as in preprint
my_color <- c("PSC"="#EFC7EA", "meso-like-1"="#BCED8E", "early-epi"="#EEC591", "iLSCs"="#D68430", "late-epi"="#E5F370", "meso-like-2"="#D8E1A5", "meso-like-3"="#7EAFFA", "PSC-like"="#C996FC")
my_cols2 <- my_color[order(names(my_color))]

pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/full_monocle3_pseudotime.pdf') ,width=5,height=4,paper='special')
plot_cells(cds = cds,color_cells_by = "pseudotime",show_trajectory_graph = T, trajectory_graph_color = "red", cell_size = 1.5, label_roots = F, label_cell_groups=F, label_leaves=F, label_branch_points=T) + ggtitle("UMAP of Pseudotime from PSC (C1-2)")
dev.off()

pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/full_umap.pdf') ,width=5,height=4,paper='special')
DimPlot(object = seur_obj, reduction = 'umap',
        cols = my_cols2)
dev.off()

pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/cell_type_umap.pdf') ,width=5,height=4,paper='special')
DimPlot(object = seur_obj, group.by = "cell_type",reduction = 'umap') +
  labs(title = "cell line")
dev.off()

pdata_cds <- pData(cds)

pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdataframe <- as.data.frame(pdata_cds)
ord_vec = c('PSC', 'PSC-like', 'meso-like-3', 'meso-like-2', 'meso-like-1', 'early-epi', 'iLSCs', 'late-epi')
pdataframe$ident <- factor(pdataframe$ident, levels = ord_vec)

pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/time_ordered_plot_full.pdf') ,width=5,height=4,paper='special')
ggplot(pdataframe, 
       aes(x = pseudotime_monocle3, 
           y = ident, colour = ident)) +
    geom_quasirandom(groupOnX = FALSE) +
    theme_classic() + scale_color_manual(values = my_color) + 
    xlab("monocle3 pseudotime") + ylab("Cluster") +
    ggtitle("Cells ordered by monocle3 pseudotime")
dev.off()

# ESC
ESC <- subset(x = seur_obj, cells = WhichCells(seur_obj, expression = cell_type=="hESC"))
cds <- as.cell_data_set(ESC)
cds <- cluster_cells(cds, reduction_method = "UMAP") #Monocle3 Pipeline
cds <- estimate_size_factors(cds)
cds <- learn_graph(cds, use_partition = F, close_loop = T)
root_group <- colnames(cds)[pData(cds)$ident == "PSC"]
cds <- order_cells(cds, root_cells = root_group)
pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/ESC_monocle3_pseudotime.pdf') ,width=5,height=4,paper='special')
plot_cells(cds = cds,color_cells_by = "pseudotime",show_trajectory_graph = T, trajectory_graph_color = "red", cell_size = 1.5, label_roots = F, label_cell_groups=F, label_leaves=F, label_branch_points=T) + ggtitle("UMAP of Pseudotime from PSC (C1-2)")
dev.off()
pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/hESC_umap.pdf') ,width=5,height=4,paper='special')
DimPlot(object = ESC, reduction = 'umap',
        cols = my_cols2)+
  labs(title = "hESC")
dev.off()

pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdataframe <- as.data.frame(pdata_cds)
pdataframe$ident <- factor(pdataframe$ident, levels = ord_vec)
pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/time_ordered_plot_hESC.pdf') ,width=5,height=4,paper='special')
ggplot(pdataframe, 
       aes(x = pseudotime_monocle3, 
           y = ident, colour = ident)) +
    geom_quasirandom(groupOnX = FALSE) +
    theme_classic() + scale_color_manual(values = my_color) + 
    xlab("monocle3 pseudotime") + ylab("Cluster") +
    ggtitle("Cells ordered by monocle3 pseudotime")
dev.off()

# BH2T
bC <- subset(x = seur_obj, cells = WhichCells(seur_obj, expression = cell_type=="hiPSC1"))
cds <- as.cell_data_set(bC)
cds <- cluster_cells(cds, reduction_method = "UMAP") #Monocle3 Pipeline
cds <- estimate_size_factors(cds)
cds <- learn_graph(cds, use_partition = F, close_loop = T)
root_group <- colnames(cds)[pData(cds)$ident == "PSC"]
cds <- order_cells(cds, root_cells = root_group)
pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/bh2t_monocle3_pseudotime.pdf') ,width=5,height=4,paper='special')
plot_cells(cds = cds,color_cells_by = "pseudotime",show_trajectory_graph = T, trajectory_graph_color = "red", cell_size = 1.5, label_roots = F, label_cell_groups=F, label_leaves=F, label_branch_points=T) + ggtitle("UMAP of Pseudotime from PSC (C1-2)")
dev.off()
pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/hiPSC1_bh2t_umap.pdf') ,width=5,height=4,paper='special')
DimPlot(object = bC, reduction = 'umap',
        cols = my_cols2)+
  labs(title = "hiPSC1")
dev.off()

pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdataframe <- as.data.frame(pdata_cds)
pdataframe$ident <- factor(pdataframe$ident, levels = ord_vec)
pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/time_ordered_plot_bh2t.pdf') ,width=5,height=4,paper='special')
ggplot(pdataframe, 
       aes(x = pseudotime_monocle3, 
           y = ident, colour = ident)) +
    geom_quasirandom(groupOnX = FALSE) +
    theme_classic() + scale_color_manual(values = my_color) + 
    xlab("monocle3 pseudotime") + ylab("Cluster") +
    ggtitle("Cells ordered by monocle3 pseudotime")
dev.off()

# 003bc
BHT <- subset(x = seur_obj, cells = WhichCells(seur_obj, expression = cell_type=="hiPSC2"))
cds <- as.cell_data_set(BHT)
cds <- cluster_cells(cds, reduction_method = "UMAP") #Monocle3 Pipeline
cds <- estimate_size_factors(cds)
cds <- learn_graph(cds, use_partition = F, close_loop = T)
root_group <- colnames(cds)[pData(cds)$ident == "PSC"]
cds <- order_cells(cds, root_cells = root_group)
pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/003bc_monocle3_pseudotime.pdf') ,width=5,height=4,paper='special')
plot_cells(cds = cds,color_cells_by = "pseudotime",show_trajectory_graph = T, trajectory_graph_color = "red", cell_size = 1.5, label_roots = F, label_cell_groups=F, label_leaves=F, label_branch_points=T) + ggtitle("UMAP of Pseudotime from PSC (C1-2)")
pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/hipsc2_003bc_umap.pdf') ,width=5,height=4,paper='special')
DimPlot(object = BHT, reduction = 'umap',
        cols = my_cols2)+
  labs(title = "hiPSC2")
dev.off()

pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
pdataframe <- as.data.frame(pdata_cds)
pdataframe$ident <- factor(pdataframe$ident, levels = ord_vec)
pdf(paste0('/mnt/s/Radboud/data/scRNAseq_ipsc/time_ordered_plot_003bc.pdf') ,width=5,height=4,paper='special')
ggplot(pdataframe, 
       aes(x = pseudotime_monocle3, 
           y = ident, colour = ident)) +
    geom_quasirandom(groupOnX = FALSE) +
    theme_classic() + scale_color_manual(values = my_color) + 
    xlab("monocle3 pseudotime") + ylab("Cluster") +
    ggtitle("Cells ordered by monocle3 pseudotime")
dev.off()
