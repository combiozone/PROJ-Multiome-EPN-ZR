library(Seurat)


seurat_obj <- readRDS('multiome.rds')

DefaultAssay(seurat_obj) <- 'SCT'
seurat_obj <- PrepSCTFindMarkers(seurat_obj) 
diff_exp <- FindAllMarkers(
        object = seurat_obj,
        only.pos = TRUE,
        logfc.threshold = 0.25,
        min.pct = 0.05,
        recorrect_umi = FALSE
    )

write.table(diff_exp, file = 'findallmarker.geneExp.xls', sep = "\t", quote=FALSE, row.names = F)



