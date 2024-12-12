####################Pando
library(Pando)
library(Seurat)
library(BSgenome.Btaurus.UCSC.bosTau9)
cattle_data <- readRDS("SeuratObject.rds")
cattle_data[['RNA']]
cattle_data[['peaks']]
cattle_data <- initiate_grn(cattle_data)
pfms <- readRDS("pfms.rds")

motif_841_name <- read.table("motif_841_name.txt",header = F )
cattle_data <- find_motifs(
    cattle_data,
    pfm = pfms,
    motif_tfs = motif_841_name[,-1],
    genome = BSgenome.Btaurus.UCSC.bosTau9
)

library(doParallel)
registerDoParallel(10)
cattle_data <- infer_grn(
    cattle_data,
    peak_to_gene_method = 'GREAT',
    parallel = T
)
coef(cattle_data)
cattle_data <- find_modules(cattle_data)
modules <- NetworkModules(cattle_data)
modules@meta
save(list=ls(),file="pando.RData")