library(cicero)
library(monocle3)
library(ggplot2)
library(patchwork)
library(stringr)
library(SeuratWrappers)
library(Signac)
library(Seurat)
set.seed(1)
#1. Load Data and create gene_anno
chromosome_length <- read.table("ARS-UCD1.2.chr30.sizes")
gene_anno <- rtracklayer::readGFF("ARS-UCD1.2.110.chr.gtf")
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

#########################
#2. Seurat to CDS Format and Run Cicero
#########################
cattle_cicero <- readRDS("cattle_cicero.rds")
conns <- run_cicero(cattle_cicero, chromosome_length)
head(conns)
#3. Finding cis-Co-accessibility Networks (CCANS)
ccans <- cicero::generate_ccans(conns)
#4. Save
saveRDS(conns, file = "./Outputs/conns.rds")
saveRDS(ccans, file = "./Outputs/ccans.rds")

all_peaks <- row.names(exprs(cds))
write.csv(x = all_peaks, file = "./Outputs/all_peaks.csv")
write.csv(x = conns, file = "./Outputs/cicero_connections.csv")

links <- ConnectionsToLinks(conns=conns, ccans=ccans)
Links(SeuratObject) <- links

save(list=ls(),file="cicero.RData")
message("all done!")
message("all done!")