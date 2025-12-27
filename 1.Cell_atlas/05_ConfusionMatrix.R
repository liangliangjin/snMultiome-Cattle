# confusionMatrix
library(pheatmap)
library(RColorBrewer)

cM_wnn_rna <- confusionMatrix(paste0(proj$Clusters_RNA), paste0(proj$WNN_Clusters))
cM_wnn_rna <- as.matrix(cM_wnn_rna / Matrix::rowSums(cM_wnn_rna))
cM_wnn_rna <- cM_wnn_rna[match(sort(rownames(cM_wnn_rna)),rownames(cM_wnn_rna)),]
max_cols <- apply(cM_wnn_rna, 1, which.max)
sorted_matrix <- cM_wnn_rna[,c(unique(max_cols),match(colnames(cM_wnn_rna)[-unique(max_cols)],colnames(cM_wnn_rna)))]


p_wnn_rna <- pheatmap::pheatmap(
  mat = as.matrix(sorted_matrix), border_color = "grey90",
  color = colorRampPalette(brewer.pal(n = 9, name ="RdPu"))(100),cluster_rows = FALSE,cluster_cols = FALSE
)


cM_wnn_atac <- confusionMatrix(paste0(proj$Clusters_ATAC), paste0(proj$WNN_Clusters))
cM_wnn_atac <- as.matrix(cM_wnn_atac / Matrix::rowSums(cM_wnn_atac))
cM_wnn_atac <- cM_wnn_atac[, colnames(sorted_matrix)[colnames(sorted_matrix) %in% colnames(cM_wnn_atac)], drop = FALSE]
sort_rows_by_values <- function(mat) {
  row_max <- apply(mat, 1, max)
  max_cols <- numeric(nrow(mat))
  for (i in 1:nrow(mat)) {
    if (row_max[i] == 0) {
      max_cols[i] <- 1
    } else {
      max_cols[i] <- which.max(mat[i, ])
    }
  }
  row_order <- order(max_cols, -row_max)
  return(mat[row_order, ])
}
sorted_matrix2 <- sort_rows_by_values(cM_wnn_atac)
p_wnn_atac <- pheatmap::pheatmap(
  mat = as.matrix(sorted_matrix2), border_color = "grey90",
  color = colorRampPalette(brewer.pal(n = 9, name ="RdPu"))(100),cluster_rows = FALSE,cluster_cols = FALSE
)

celltype2color <- setNames(colorRampPalette(ArchRPalettes$stallion)(length(unique(proj$main))), sort(unique(proj$main)))
Annotation <- read.table("Annotation_result.txt", header = T, sep = "\t")
cluster2color <- setNames(my_colors[Annotation$annotation], Annotation$celltype)
my_colors <- cluster2color[match(as.numeric(colnames(sorted_matrix)), names(cluster2color))]
my_colors[is.na(my_colors)] <- celltype2color["Unknown cells"]
image(matrix(1:length(my_colors)), col = my_colors, axes = FALSE)