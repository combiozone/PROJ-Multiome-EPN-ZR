# Hua Sun
# Seurat v5.1


library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(stringi)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)


rdir <- 'out_cellranger_arc'
seed <- 42
ref <- 'mm10'
regress <- NULL  # set it by data
outdir <- 'out_multiome_integrated'

set.seed(seed)
source('funcs.R')

dir.create(outdir)


macs2 <- '~/software/miniconda3/bin/macs2'
ensdb <- EnsDb.Mmusculus.v79
if (ref == 'hg38'){ ensdb <- EnsDb.Hsapiens.v86 }
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- ref
genome <- BSgenome.Mmusculus.UCSC.mm10
if (ref == 'hg38'){ genome <- BSgenome.Hsapiens.UCSC.hg38 }


multiome_filtered_list <- readRDS('multiome_filtered.list.rds')

print('[INFO] macs2 ...')
snatac_macs2_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- CallPeaksUsingMACS2(x, macs2, annotation) })
saveRDS(snatac_macs2_list, file=paste0(outdir, '/snatac_macs2.list.rds'))

print('[INFO] combine peaks ...')
snatac_peak_list <- ObjectListConvertBedToGRanges(snatac_macs2_list)
snatac_peak_list <- lapply(X = snatac_peak_list, FUN = function(x) { x <- ComputeLSI(x, 'peaks', 10) })
saveRDS(snatac_peak_list, file=paste0(outdir, '/snatac_peaks.list.rds'))

integ_atac <- ATACIntegration(snatac_list=snatac_peak_list, max_dim=30, seed=seed)
saveRDS(integ_atac, file=paste0(outdir,'/snatac_integ.withumap.rds'))

print('[INFO] integrate multiome ...')
integrated_multiome <- MultiomeIntegration_v1(obj_list=multiome_filtered_list, regress=regress, snatac_obj=integ_atac)
saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.sct_lsi.rds'))








