library(Seurat)


seurat_obj <- readRDS('multiome.rds')

DefaultAssay(seurat_obj) <- 'SCT'
seurat_obj <- CellCycleScoring(seurat_obj, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

saveRDS(seurat_obj, 'multiome.cc.rds')
