library(Seurat)
library(ArchR)
addArchRThreads(threads = 16)

proj <- readRDS("combine_all_after_filter.rds")

proj <- addImputeWeights(proj, reducedDims = "Harmony")

dir.create("0.marker", showWarnings = FALSE)
markersGS <- getMarkerFeatures(
	ArchRProj = proj,
	useMatrix = "GeneExpressionMatrix",
	groupBy = "Clusters",
	bias = c("TSSEnrichment","log10(Gex_nUMI)"),
	testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 2")
for(k in 1:length(unique(proj$Clusters))){
write.table(markerList@listData[[k]],paste0(getwd(), "/0.marker/",length(labelOld),"/combine_all_",names(markerList_rna[k]),".txt"),col.names=T,row.names=T,sep="\t");
}

markersGS_RNA <- getMarkerFeatures(
	ArchRProj = proj,
	useMatrix = "GeneExpressionMatrix",
	groupBy = "Clusters_RNA",
	bias = c("TSSEnrichment","log10(Gex_nUMI)"),
	testMethod = "wilcoxon"
)
markerList_RNA <- getMarkers(markersGS_RNA, cutOff = "FDR <= 0.01 & Log2FC >= 2")
for(k in 1:length(unique(proj$Clusters_RNA))){
write.table(markerList_RNA@listData[[k]],paste0(getwd(), "/0.marker/",length(labelOld),"/RNA_all_",names(markerList_RNA[k]),".txt"),col.names=T,row.names=T,sep="\t");
}


################All clusters of marker genes were examined and annotated manually in conjunction with the literature

################Star genes can also be visualized on UMAP
plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneExpressionMatrix",#GeneScoreMatrix GeneExpressionMatrix
  name = "GAD1",
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  size = 0.00000000001,rastr = FALSE,
  imputeWeights = getImputeWeights(proj)
)


#Rename celltypes
labelOld <- unique(proj$Clusters)
labelNew <- readLines("Annotation_result.txt")#Needs to correspond to the order in labelOld
proj$main <- mapLabels(proj$Clusters, newLabels = labelNew, oldLabels = labelOld)

#Visualization
p_main <- plotEmbedding(proj, name = "main",embedding ="UMAP",labelAsFactors=F,size = 0.000001,labelMeans=F,legendSize=10,rastr = FALSE)

p_rna <- plotEmbedding(proj, name = "main",embedding ="UMAP_RNA_0.6",labelAsFactors=F,size = 0.000001,labelMeans=F,legendSize=10,rastr = FALSE)

p_atac <- plotEmbedding(proj, name = "main",embedding ="UMAP_ATAC_0.6",labelAsFactors=F,size = 0.000001,labelMeans=F,legendSize=10,rastr = FALSE)