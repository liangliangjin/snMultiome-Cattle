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
#2.Run Seurat to CDS Format
#########################

cattle_cicero <- readRDS("cattle_cicero.rds")

head(conns)
#3. Finding cis-Co-accessibility Networks (CCANS)
ccans <- cicero::generate_ccans(conns)
#4. Save
saveRDS(conns, file = "./Outputs/conns.rds")
saveRDS(ccans, file = "./Outputs/ccans.rds")

all_peaks <- row.names(exprs(brain_cds))
write.csv(x = all_peaks, file = "./Outputs/all_peaks.csv")
write.csv(x = conns, file = "./Outputs/cicero_connections.csv")
message("all done!")