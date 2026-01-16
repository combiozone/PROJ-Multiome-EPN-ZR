
# https://github.com/combiozone/CellTypeEstimate

# mouse
Rscript CellTypeEstimate/cte4os.R --rds multiome.rds --db mm.brain.v2 --assay SCT --reduction 'wnn.umap' --plot --outdir out_celltype 

# human
#Rscript CellTypeEstimate/cte4os.R --rds multiome.rds --db hs.brain --assay SCT --reduction 'wnn.umap' --plot --outdir out_celltype 


