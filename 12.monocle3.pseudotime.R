# Hua Sun

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(dplyr)
library(data.table)

set.seed(42)

outdir <- 'out_monocle3'
dir.create(outdir)

rds <- 'seurat5.1_v6.2/multiome_integrated.plus.rds'
seu <- readRDS(rds)

DefaultAssay(seu) <- 'SCT'
Idents(seu) <- 'cell_type2'

seu[["UMAP"]] <- seu[['wnn.umap']]
cds <- as.cell_data_set(seu)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)
cds <- order_cells(cds)


plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)






