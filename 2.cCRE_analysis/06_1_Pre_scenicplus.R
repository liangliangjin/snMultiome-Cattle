library(ArchR)
library(Matrix)
proj <- readRDS("combine_all_new_fig3.rds")
peak_matrix <- getMatrixFromProject(proj, useMatrix = 'PeakMatrix')
peak_counts <- peak_matrix@assays@data$PeakMatrix
peak_meta = as.data.frame(peak_matrix@rowRanges)
peaks = tidyr::unite(peak_meta, "peaks", seqnames, start, end)
peaks_regions = peaks$peaks
dim(peak_counts)

writeMM(peak_counts, "atac_counts.mtx")
write.table(peaks_regions, file = "peaks_regions.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(colnames(peak_counts), file = "cell_names.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

rna_matrix <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")                                                                                                       
rna_counts <- rna_matrix@assays@data$GeneExpressionMatrix
features <- rowData(rna_matrix)$name

writeMM(rna_counts, "rna_counts.mtx")
write.table(features, file = "gene_name.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

if (all(colnames(rna_counts) == colnames(peak_counts))) {
  message("The barcode names for RNA and ATAC are exactly the same.")
} else {
  write.table(colnames(rna_counts), file = "cell_names2.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  message("RNA barcodes saved to cell_names2.tsv")
}


## Add header to mtx files
#sed -i '1i%%MatrixMarket matrix coordinate integer general\n' atac_counts.mtx
#sed -i '1i%%MatrixMarket matrix coordinate integer general\n' rna_counts.mtx