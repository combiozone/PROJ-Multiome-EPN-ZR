# Hua Sun


library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)


fmd <- ''
control <- ''

cutoff_region_size <- 5
min_ratio <- 0.5
min_cnv <- 100
min_gene <- 20

cutoff_gain <- 1.2
cutoff_loss <- 0.8

rm_clu <- ''
rm_chr <- ''

icnv_exp <- ''
hmm_gene <- ''
hmm_region <- ''

outdir <- 'out_malignant_cells'

dir.create(outdir)


df_cnv_region <- read.table(hmm_region, sep='\t', header=T)
df_cnv_region$region_size <- (df_cnv_region$end - df_cnv_region$start)/1000000


df_cnv_region <- df_cnv_region %>% filter(region_size >= cutoff_region_size)
filtered_cnv_region <- unique(df_cnv_region[['cnv_name']])
#print(filtered_cnv_region)

df_cnv_gene <- read.table(hmm_gene, sep='\t', header=T)

df_cnv_gene <- df_cnv_gene %>% filter(gene_region_name %in% filtered_cnv_region)
print(dim(df_cnv_gene))


if (nchar(rm_chr) > 0){
    remove_chr <- str_split(rm_chr, ',')[[1]]
    df_cnv_gene <- df_cnv_gene %>% filter(!(chr %in% remove_chr))
}
print(unique(df_cnv_gene$chr))



filtered_cnv_gene_set <- unique(df_cnv_gene$gene)

d_icnv_exp <- readRDS(icnv_exp)
used_cells <- rownames(d_icnv_exp)

d_icnv_exp_filtered <- d_icnv_exp[,filtered_cnv_gene_set]
print(dim(d_icnv_exp_filtered))


# sample name
df_cnv_gene$sample <- gsub("\\.\\w+", "", df_cnv_gene$cell_group_name, perl=TRUE)


# read metadata
metadata <- read.table(fmd, sep='\t', header=T)
rownames(metadata) <- metadata$cell
metadata$cell <- NULL
print(dim(metadata))

metadata_filtered <- metadata[used_cells,]


if (nchar(rm_clu) > 0){
    del_clusters <- str_split(rm_clu, ',')[[1]]
    metadata_filtered <- metadata_filtered %>% filter(!(seurat_clusters %in% del_clusters))
}

print(unique(metadata_filtered$seurat_clusters))
print(nrow(metadata_filtered))


df_merged <- ''
i <- 0
for (s in unique(df_cnv_gene$sample)){
    ds <- df_cnv_gene[df_cnv_gene$sample==s,]

    ds_report <- metadata_filtered[metadata_filtered$sample==s, c('sample', 'seurat_clusters')]
    cells <- row.names(ds_report)

    s_regions <- c()
    s_cnv_gene_count <- c()
    s_cnv_ratio <- c()

    # region cnv
    for (r in unique(ds$gene_region_name)){

        genes <- ds$gene[ds$gene_region_name==r]
        
        n_gene_region <- length(genes)
        
        if (n_gene_region < min_gene){
            print(paste('skip ', r, '-', n_gene_region))
            next
        }

        ds_exp <- d_icnv_exp_filtered[cells, genes]
        total_genes <- length(genes)
        gain <- rowSums(ds_exp >= cutoff_gain)
        loss <- rowSums(ds_exp <= cutoff_loss)

        s_regions <- cbind(s_regions, r)

        cnv_genes <- gain + loss
        s_cnv_gene_count <- cbind(s_cnv_gene_count, cnv_genes)
        colnames(s_cnv_gene_count) <- as.vector(s_regions)

        cnv_ratio <- (gain+loss)/total_genes
        s_cnv_ratio <- cbind(s_cnv_ratio, cnv_ratio)
        colnames(s_cnv_ratio) <- as.vector(s_regions)
    }
    
    ds_report$total_cnv.regions <- length(as.vector(s_regions))
    ds_report$cnv_gene_count <- as.vector(unite(col="merged_cnv_genes", as.data.frame(s_cnv_gene_count), sep=","))
    ds_report$cnv_ratio <- as.vector(unite(col="merged_cnv_ratio", as.data.frame(s_cnv_ratio), sep=","))
    
    ds_report$cnv_gene_count.mean <- round(apply(s_cnv_gene_count, 1, mean, na.rm=TRUE), 0)
    ds_report$cnv_gene_count.max <- apply(s_cnv_gene_count, 1, max, na.rm=TRUE)
    ds_report$cnv_ratio.mean <- apply(s_cnv_ratio, 1, mean, na.rm=TRUE)
    ds_report$cnv_ratio.max <- apply(s_cnv_ratio, 1, max, na.rm=TRUE)

    ds_report$cnv_ratio.max <- round(ds_report$cnv_ratio.max, 4)
    ds_report$cell_status <- 'Non-malignant'
    ds_report$cell_status[ds_report$cnv_ratio.max >= min_ratio | ds_report$cnv_gene_count.max >= min_cnv] <- 'Malignant'

    i <- i + 1
    if (i==1) {
        df_merged <- ds_report
    } else {
        df_merged <- rbind(df_merged, ds_report)
    }
    
}


print(dim(df_merged))
df_merged <- df_merged %>% select(sample, seurat_clusters, total_cnv.regions, cnv_gene_count.mean, cnv_gene_count.max, cnv_ratio.mean, cnv_ratio.max, cell_status)


# cell g1 g2 g3 ...
total_genes <- ncol(d_icnv_exp)
print(total_genes)

gain <- rowSums(d_icnv_exp >= cutoff_gain)
loss <- rowSums(d_icnv_exp <= cutoff_loss)
total_cnv_genes <- gain + loss

gain_ratio.by_allgenes <- round(gain/total_genes * 100, 2)
loss_ratio.by_allgenes <- round(loss/total_genes * 100, 2)
cnv_ratio.by_allgenes <- round(total_cnv_genes/total_genes * 100, 2)

d_cnv_ratio2 <- data.frame(gain_ratio.by_allgenes, loss_ratio.by_allgenes, cnv_ratio.by_allgenes)
rownames(d_cnv_ratio2) <- row.names(d_icnv_exp)

df_merged <- merge(df_merged, d_cnv_ratio2, by='row.names', all.x=TRUE)

# output
write.table(df_merged, paste0(outdir, '/infercnv.estimate_malignant_cells.', cutoff_region_size, 'Mb.xls'), sep='\t', quote=FALSE, row.names=FALSE)



ref_clu <- stringr::str_split(control, ',')[[1]]
d_poor_ref <- df_merged %>% 
            filter(seurat_clusters %in% ref_clu) %>% 
            filter(cell_status=='Malignant')

poor_ref_cells <- d_poor_ref[,1]
writeLines(poor_ref_cells, paste0(outdir, '/poor_ref_cells.out'))




