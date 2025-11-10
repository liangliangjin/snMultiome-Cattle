library(Seurat)
library(ArchR)
addArchRThreads(threads = 16)

proj <- readRDS("combine_all_after_filter.rds")

proj <- addImputeWeights(proj, reducedDims = "Harmony")

dir.create("./0.marker/WNN_Clusters", showWarnings = FALSE)
proj$WNN_Clusters[proj$WNN_Clusters %in% names(table(proj$WNN_Clusters)[table(proj$WNN_Clusters) < 10])] <- "Unknown cells"
markersGS <- getMarkerFeatures(
	ArchRProj = proj,
	useMatrix = "GeneExpressionMatrix",
	groupBy = "WNN_Clusters",
	bias = c("TSSEnrichment","log10(Gex_nUMI)"),
	testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
for(c in 1:length(unique(proj$WNN_Clusters))){
write.table(markerList@listData[[c]],paste0(getwd(), "/0.marker/WNN_Clusters/combine_all_",names(markerList[c]),".txt"),col.names=T,row.names=T,sep="\t");
}

################All clusters of marker genes were examined and annotated manually in conjunction with the literature

################Star genes can also be visualized on UMAP
plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneExpressionMatrix",#GeneScoreMatrix GeneExpressionMatrix
  name = "GAD1",
  embedding = "UMAP_wnn",
  quantCut = c(0.01, 0.95),
  size = 0.00000000001,rastr = FALSE,
  imputeWeights = getImputeWeights(proj)
)


#Rename celltypes
labelOld <- as.character(unique(proj$WNN_Clusters))
Annotation <- read.table("Annotation_result.txt", header = T, sep = "\t")
mapping <- setNames(as.character(Annotation$annotation), as.character(Annotation$celltype))
labelNew <- ifelse(!is.na(mapping[labelOld]), mapping[labelOld], "Unknown cells")
proj$main <- mapLabels(proj$WNN_Clusters, newLabels = labelNew, oldLabels = labelOld)

#Visualization
p_main <- plotEmbedding(proj, name = "main",embedding ="UMAP_wnn",labelAsFactors=F,size = 0.000001,labelMeans=F,legendSize=5,rastr = FALSE)+theme(legend.text = element_text(size = 10))+guides(colour  = guide_legend(override.aes = list(shape = 20,size=8)))

p_cattle <- plotEmbedding(proj, name = "Clusters3",embedding ="UMAP_wnn",labelAsFactors=F,size = 0.000001,labelMeans=F,legendSize=5,rastr = FALSE)
p_tissue <- plotEmbedding(proj, name = "Clusters4",embedding ="UMAP_wnn",labelAsFactors=F,size = 0.000001,labelMeans=F,legendSize=5,rastr = FALSE)
saveRDS(proj, file = "combine_all_new_fig1_end.rds")