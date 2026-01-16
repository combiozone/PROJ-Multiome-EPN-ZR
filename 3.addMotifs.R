# Hua Sun

library(TFBSTools)
library(JASPAR2022)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)



multiome <- readRDS('multiome.rds')
ref <- 'hg38'

genome <- BSgenome.Mmusculus.UCSC.mm10
if (ref == 'hg38'){ genome <- BSgenome.Hsapiens.UCSC.hg38 }

pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

DefaultAssay(multiome) <- 'peaks'
multiome <- AddMotifs(object = multiome, genome = genome, pfm = pfm)
multiome <- RunChromVAR(object = multiome, genome = genome)

saveRDS(multiome, file='multiome.addMotifs.rds')

