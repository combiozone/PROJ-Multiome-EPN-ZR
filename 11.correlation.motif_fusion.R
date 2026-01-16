# Hua Sun

library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(stringr)
library(gprofiler2) # change name
library(ggpubr)


rds <- 'multiome_integrated_plus.rds'
fmeta <- 'cluster_cellType.corrected2.xls'
gene <- 'Plagl1'
motif_id <- 'MA1615.1'
fzr <- 'out_zrFusSig93/metadata_with_featuresSig.xls'
outdir <- 'out_cor.motif_fus'


dir.create(outdir)

seurat_obj <- readRDS(rds)
meta <- fread(fmeta, data.table=F)
seurat_obj$cell_type2 <- meta$cell_type2[match(seurat_obj$seurat_clusters, meta$seurat_clusters)]

# target motif-id
d_motif <- data.frame(seurat_obj@assays$chromvar@data[motif_id,])
colnames(d_motif) <- motif_id


metadata <- seurat_obj@meta.data

d_motif <- merge(d_motif, select(metadata, c('Sample', 'cell_type2')), by='row.names', all.x=TRUE)
row.names(d_motif) <- d_motif[,1]
d_motif[,1] <- NULL


# read fusion data & merge
d_fus <- read.table(fzr, sep='\t', header=T, row.names=1)
d_motif$zr_score1 <- d_fus$zr_score1[match(row.names(d_motif), row.names(d_fus))]



d_motif2 <- filter(d_motif, cell_type2 %in% c('CycProg-Like', 'Fibroblast-Like', 'Neuronal-Like', 'RGC-Like'))

p <- ggscatter(d_motif2, x = 'zr_score1', y = motif_id, 
            color='#2278B5', shape = 16, size = 0.5, alpha=0.6,
            add = "reg.line", add.params = list(color = "#E14C32", fill = "lightgray", size=0.5),
            conf.int = TRUE, cor.coef = TRUE, cor.coef.size = 2.5,
            cor.coeff.args = list(method = "pearson"), ggtheme=clean_theme()
           ) + 
        theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(
           axis.ticks = element_line(linewidth = 0.3),
           axis.ticks.length=unit(1, "mm")) +
        labs(title='', x='ZR-Fusion signal', y=paste0(gene, ' motif score')) +
        theme(plot.title = element_text(hjust = 0.5, size=10)) +
        theme(text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 6)) 

p <- p + facet_wrap(vars(cell_type2)) +
        theme(
            strip.text = element_text(size = 8, color = "black", face = "bold"),
            strip.background = element_rect(color=NA, fill=NA)
        )

ggsave(paste0(outdir, '/', motif_id, '.cor.motif_fus.perCellType.pdf'), w=2.7, h=2.5, useDingbats=F)



