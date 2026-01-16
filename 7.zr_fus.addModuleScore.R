library(Seurat)


seurat_obj <- readRDS('multiome.rds')
geneset <- readLines('zr_fus.gene')

seurat_obj <- AddModuleScore(
    object = seurat_obj,
    assay = 'SCT',
    slot = 'data',
    features = list(geneset),
    name = 'zr_score'
)

saveRDS(seurat_obj, 'multiome.zrfus.rds')

