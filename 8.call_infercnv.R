# Hua Sun

library(infercnv)

rds <- 'pre_infercnv_obj.rds'
func <- 'hmm'
sd_amp <- 2
cutoff <- 0.1

outdir <- 'out_infercnv'

dir.create(outdir)


infercnv_obj <- readRDS(rds)
infercnv_obj2 <- ''

# fast
if (func == 'default'){
    infercnv_obj2 <- infercnv::run(
        infercnv_obj,
        cutoff=cutoff, 
        out_dir=outdir,
        output_format='pdf',
        write_expr_matrix=FALSE,
        cluster_by_groups=TRUE, 
        denoise=TRUE,
        sd_amplifier = sd_amp,
        noise_logistic=TRUE,
        analysis_mode="samples",
        HMM=FALSE,
        no_plot=FALSE,
        no_prelim_plot=TRUE,
        hclust_method="ward.D2",
        num_threads = 4
    )
}



# very slow 
if (func == 'hmm'){
    infercnv_obj2 <- infercnv::run(
        infercnv_obj,
        cutoff=cutoff, 
        out_dir=outdir,
        output_format='pdf',
        write_expr_matrix=FALSE,
        cluster_by_groups=TRUE, 
        denoise=TRUE,
        sd_amplifier = sd_amp,
        noise_logistic=TRUE,
        analysis_mode="samples",
        HMM=TRUE,
        no_plot=FALSE,
        no_prelim_plot=TRUE,
        hclust_method="ward.D2",
        num_threads = 8
    )
}



# very slow 
if (func == 'hmm-subc'){
    infercnv_obj2 <- infercnv::run(
        infercnv_obj,
        cutoff=cutoff,
        out_dir=outdir,
        output_format='pdf',
        cluster_by_groups=TRUE,
        plot_steps=FALSE,
        denoise=TRUE,
        sd_amplifier=sd_amp,
        noise_logistic=TRUE,
        HMM=TRUE,
        no_prelim_plot=TRUE,
        analysis_mode='subclusters',
        hclust_method="ward.D2",
        tumor_subcluster_pval=0.05,
        tumor_subcluster_partition_method='random_trees',
        num_threads = 16
    )
}



