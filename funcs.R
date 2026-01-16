# Hua Sun


CallPeaksUsingMACS2 <- function(obj=NULL, macs2='MACS2', annotation=NULL)
{
    DefaultAssay(obj) <- "ATAC"
    ref <- unique(obj$ref)

    effective_genome_size <- 2.3e+09

    if (ref == 'hg38'){ effective_genome_size <- 2.7e+09 }
    
    peaks <- CallPeaks(obj, macs2.path = macs2, effective.genome.size = effective_genome_size)

    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    if (ref == 'mm10'){ peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE) }
    if (ref == 'hg38'){ peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE) }

    # quantify counts in each peak
    macs2_counts <- FeatureMatrix(
        fragments = Fragments(obj),
        features = peaks,
        cells = colnames(obj)
    )

    obj[["peaks"]] <- CreateChromatinAssay(
        counts = macs2_counts,
        fragments = Fragments(obj),
        annotation = annotation
    )
   
    return(obj)
}



ObjectListConvertBedToGRanges <- function(snatac_macs2_list)
{
    atac_peak <- c()
    for (obj in snatac_macs2_list){
        name <- unique(obj$Sample)

        # make peak bed for multiome data
        peaks <- as.data.frame(obj@assays$peaks@ranges)[,1:3]
        colnames(peaks) <- c('chr', 'start', 'end')
        
        # convert to genomic ranges
        gr <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(peaks))
        atac_peak[[name]] <- gr
    }
    print(names(atac_peak))

    # Create a unified set of peaks to quantify in each dataset
    combined.peaks <- GenomicRanges::reduce(x = unlist(GRangesList(atac_peak)))

    # Filter out bad peaks based on length
    peakwidths <- width(combined.peaks)
    combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]

    snatac_peak_list <- lapply(X = snatac_macs2_list, FUN = function(x) { x <- MakeNewATACObject(x, combined.peaks, annotation) })

    return(snatac_peak_list)
}


RenameObjectListCells <- function(objectlist)
{
    objectlist <- base::lapply(X = objectlist, FUN = function(x) { x <- RenameCells(x, new.names = paste0(as.character(unique(x$orig.ident)), "_", Cells(x=x))) })

    return(objectlist)
}



ComputeLSI <- function(obj=NULL, assay='peaks', min_cutoff=10)
{   
    DefaultAssay(obj) <- assay
    obj <- FindTopFeatures(obj, min.cutoff=min_cutoff)
    obj <- RunTFIDF(obj)
    obj <- RunSVD(obj)

    return(obj)
}


MergeObjects <- function(objectlist)
{
    samples <- c()
    for (x in objectlist){ samples <- append(samples, unique(x$Sample)) }
    combined_object <- merge(x=objectlist[[1]], y=objectlist[2:length(objectlist)], add.cell.ids = samples)

    return(combined_object)
}


RunUMAP_Plus <- function(
    obj=NULL, 
    method='uwot',
    reduc='pca',
    reduc_name='umap',
    reduc_key='UMAP_',
    min_dim=1,
    max_dim=30,
    seed=47,
    res_clus=0.4
){
    obj <- RunUMAP(obj, umap.method=method, reduction = reduc, reduction.name = reduc_name, reduction.key = reduc_key, dims = min_dim:max_dim, seed.use = seed)
    obj <- FindNeighbors(obj, reduction = reduc, dims = min_dim:max_dim)
    obj <- FindClusters(obj, resolution = res_clus)

    return(obj)
}




ATACIntegration <- function(snatac_list=NULL, assay_use='peaks', new_reduc='integrated_lsi', min_cutoff=10, max_dim=40)
{
    combined_snatac <- MergeObjects(snatac_list)
    DefaultAssay(combined_snatac) <- assay_use

    combined_snatac <- ComputeLSI(combined_snatac, assay_use, min_cutoff)
    snatac_list <- RenameObjectListCells(snatac_list)

    integration.anchors <- FindIntegrationAnchors(
            object.list = snatac_list,
            anchor.features = rownames(snatac_list[[1]]),
            reduction = "rlsi",
            dims = 2:max_dim
    )

    integrated <- IntegrateEmbeddings(
            anchorset = integration.anchors,
            reductions = combined_snatac[["lsi"]],
            new.reduction.name = new_reduc,
            dims.to.integrate = 1:max_dim
    )

    integrated[['lsi']] <- combined_snatac[['lsi']]

    return(integrated)
}


ATACIntegration <- function(snatac_list=NULL, assay_use='peaks', new_reduc='integrated_lsi', min_cutoff=10, max_dim=40, res_clus=0.4, seed=47)
{
    combined_snatac <- MergeObjects(snatac_list)
    DefaultAssay(combined_snatac) <- assay_use

    combined_snatac <- ComputeLSI(combined_snatac, assay_use, min_cutoff)
    snatac_list <- RenameObjectListCells(snatac_list)

    integration.anchors <- Seurat::FindIntegrationAnchors(
            object.list = snatac_list,
            anchor.features = rownames(snatac_list[[1]]),
            reduction = "rlsi",
            dims = 2:max_dim
    )

    integrated <- Seurat::IntegrateEmbeddings(
            anchorset = integration.anchors,
            reductions = combined_snatac[["lsi"]],
            new.reduction.name = new_reduc,
            dims.to.integrate = 1:max_dim
    )

    integrated[['lsi']] <- combined_snatac[['lsi']]
    integrated <- RunUMAP_Plus(obj=integrated, reduc=new_reduc, reduc_name='atac.umap', reduc_key='atacUMAP_', min_dim=2, max_dim=max_dim, seed=seed, res_clus=res_clus)

    return(integrated)
}




MultiomeIntegration_v1 <- function(
    obj_list=NULL, 
    snatac_obj=NULL, 
    sct_ver = "v2", 
    regress = NULL,
    var_feature=3000, 
    reduc_rna = 'pca',
    rna_umap = 'rna.umap',
    reduc_atac='integrated_lsi',
    dim_max_rna=30, 
    dim_max_atac=30, 
    weight_name = "RNA.weight",
    res_clus=0.4,
    seed=42
){
    multiome <- snatac_obj

    merged_obj <- MergeObjects(obj_list)
    multiome[['RNA']] <- merged_obj[['RNA']]
    DefaultAssay(multiome)<-"RNA"

    multiome <- SCTransform(multiome, vst.flavor = sct_ver, vars.to.regress = regress, variable.features.n = var_feature, verbose = FALSE)
    multiome <- RunPCA(multiome, seed.use = seed, verbose = TRUE)
    multiome <- RunUMAP(multiome, dims = 1:dim_max_rna, reduction = reduc_rna, reduction.name = rna_umap, seed.use = seed)
    multiome <- FindMultiModalNeighbors(
                    object = multiome, 
                    reduction.list = list(reduc_rna, reduc_atac), 
                    dims.list = list(1:dim_max_rna, 2:dim_max_atac), 
                    modality.weight.name = weight_name,
                    verbose = TRUE
                )
    multiome <- RunUMAP(
                    object = multiome, 
                    nn.name = "weighted.nn", 
                    reduction.name = "wnn.umap", 
                    reduction.key = "wnnUMAP_",
                    seed.use = seed,
                    verbose = TRUE
                )
    multiome <- FindClusters(multiome, graph.name="wsnn", resolution=res_clus)

    return(multiome)
}






