setwd("/home/Jingliangliang/SC")
load("muscle-combine.RData")
#######################
#After ArchR2Seurat
#######################
library(ArchR)
library(Seurat)
library(SingleCellExperiment)
library(Signac)
library(AnnotationHub)
library(harmony)

DefaultAssay(SeuratObject) <- "RNA"
SeuratObject <- SCTransform(SeuratObject, verbose = FALSE) %>% 
	RunPCA() %>% 
	RunHarmony(group.by.vars = "Sample", reduction = "pca", reduction.save = "harmony_rna",dims.use = 1:30) %>% 
	RunUMAP(dims = 1:30, reduction = "harmony_rna", reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

colnames(SeuratObject) <- sub("_", "#", colnames(SeuratObject))
LSI_ATAC<- getReducedDims(proj, reducedDims ="Harmony_ATAC", returnMatrix = TRUE)

if(all(rownames(LSI_ATAC) %in% colnames(SeuratObject))){
	SeuratObject[["harmony_atac"]] <- CreateDimReducObject(
		embeddings = as.matrix(LSI_ATAC),
		key = "LSIATAC_",
		assay = "peaks"
	)}

set.seed(1)
options(future.seed = TRUE)
SeuratObject <- FindMultiModalNeighbors(SeuratObject, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:30, 2:30))
SeuratObject <- RunUMAP(SeuratObject, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
SeuratObject <- FindClusters(SeuratObject, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
markers <- FindAllMarkers(SeuratObject, only.pos = TRUE)
markers %>%
	group_by(cluster) %>%
	dplyr::filter(avg_log2FC > 1)

write.csv(markers, "markers-2025.csv", row.names = FALSE, quote = FALSE)


library(readr)
markers_ref <- read_csv("marker_holstein.csv")
marker_list <- strsplit(markers_ref$Markergenes, ",\\s*")
names(marker_list) <- markers_ref$Celltype
all_genes <- unique(unlist(marker_list))
ref_mat <- matrix(0, nrow = length(all_genes), ncol = length(marker_list),
                  dimnames = list(all_genes, names(marker_list)))
				  
for (ct in names(marker_list)) {
  genes <- marker_list[[ct]]
  genes <- gsub("[^A-Za-z0-9]", "_", genes)
  genes <- genes[genes != ""] 
  genes <- unique(genes)
  genes <- intersect(genes, rownames(ref_mat))
  marker_list[[ct]] <- genes
}
for(ct in names(marker_list)){
  ref_mat[marker_list[[ct]], ct] <- 1
}
se_ref <- SummarizedExperiment::SummarizedExperiment(
  assays = list(logcounts = as.matrix(ref_mat)),
  colData = DataFrame(label = colnames(ref_mat))
)
#run SingleR	
query_mat <- as.matrix(GetAssayData(SeuratObject, assay = "RNA", slot = "data"))
common_genes <- intersect(rownames(query_mat), rownames(ref_mat))
query2 <- query_mat[common_genes, , drop=FALSE]
ref2   <- ref_mat[common_genes, , drop=FALSE]

pred <- SingleR(test = query2, ref = ref2, labels = colnames(ref2), method = "single", cluster = SeuratObject@meta.data$seurat_clusters)
SeuratObject$SingleR_marker_label <- pred$labels





FeaturePlot(SeuratObject, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))

###
seurat.umap <- SeuratObject@reductions$wnn.umap
seurat.umap.df <- DataFrame(row.names=proj$cellNames, "seurat#UMAP1" = seurat.umap@cell.embeddings[proj$cellNames, "wnnUMAP_1"], "seurat#UMAP2" =  seurat.umap@cell.embeddings[proj$cellNames, "wnnUMAP_2"], check.names = FALSE)
seurat.umap.df$cn <- rownames(seurat.umap.df)

ArchR.meta <- proj@cellColData
ArchR.meta$cn <- rownames(ArchR.meta)

ordered.coords <- merge(ArchR.meta, seurat.umap.df, by = "cn")[, c("seurat#UMAP1", "seurat#UMAP2", "cn")]
rownames(ordered.coords) <- ordered.coords$cn
ordered.coords$cn <- NULL

proj@embeddings$UMAP_wnn <- SimpleList(df = seurat.umap.df, params = list())
###






DimPlot(SeuratObject, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
DimPlot(SeuratObject, reduction = "wnn.umap", group.by = "Sample", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
DimPlot(SeuratObject, reduction = "wnn.umap", cells.highlight = old2$`rownames(SeuratObject@meta.data)`[old2$`SeuratObject@meta.data$seurat_clusters`=="50"] ,label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")




