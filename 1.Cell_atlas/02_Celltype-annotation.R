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
markerList <- getMarkers(markersGS, cutOff = "FDR < 0.05 & Log2FC > 1")
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
  quantCut = c(0.01, 0.95),rastr = FALSE,
  imputeWeights = getImputeWeights(proj)
)


#Rename celltypes
labelOld <- as.character(unique(proj$WNN_Clusters))
Annotation <- read.table("Annotation_result.txt", header = T, sep = "\t")
mapping <- setNames(as.character(Annotation$annotation), as.character(Annotation$celltype))
labelNew <- ifelse(!is.na(mapping[labelOld]), mapping[labelOld], "Unknown cells")
proj$main <- mapLabels(proj$WNN_Clusters, newLabels = labelNew, oldLabels = labelOld)

#Visualization
p_main <- plotEmbedding(proj, name = "main",embedding ="UMAP_wnn",labelAsFactors=F,labelMeans=F,legendSize=5,rastr = FALSE)+theme(legend.text = element_text(size = 10))+guides(colour  = guide_legend(override.aes = list(shape = 20,size=8)))

p_cattle <- plotEmbedding(proj, name = "Clusters3",embedding ="UMAP_wnn",labelAsFactors=F,labelMeans=F,legendSize=5,rastr = FALSE)
p_tissue <- plotEmbedding(proj, name = "Clusters4",embedding ="UMAP_wnn",labelAsFactors=F,labelMeans=F,legendSize=5,rastr = FALSE)
saveRDS(proj, file = "combine_all_new_fig1_end.rds")



#hclust
#############################GeneExpressionMatrix
markersGS <- getMarkerFeatures(
	ArchRProj = proj,
	useMatrix = "GeneExpressionMatrix",
	groupBy = "main",
	bias = c("log10(Gex_nUMI)"),
	testMethod = "wilcoxon"
)

heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.05 & Log2FC >= 1.25",
    #labelMarkers = markerGenes,
    transpose = TRUE
)

hclust_mat <- as.matrix(heatmapGS@matrix)[!rownames(mat) %in% c("Unknown cells"), ]
hclust_data <- hclust(dist(hclust_mat),method="ward.D2")

plot(hclust_data,hang=-1)

#order: hclust_data$order
#label: hclust_data$labels[hclust_data$order]
ratio_data <- as.data.frame(proj@cellColData) %>%
  filter(!main %in% c("Unknown cells")) %>%
  count(Clusters4, main) %>%
  as.data.frame()

names(ratio_data) <- c("Organ","Celltype","CellNumber")
ratio_data$Celltype <- factor(ratio_data$Celltype, levels = hclust_data$labels[hclust_data$order])
bar_plot_main<-ggplot(ratio_data, aes(x = Celltype, y = CellNumber)) +
    geom_col(aes(fill = Organ), color = "black", position = 'fill') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid = element_blank()) +
    scale_fill_brewer(palette = "Set3")



SeuratObject$main <- proj$main
SeuratObject_filtered <- subset(SeuratObject, subset = main != "Unknown cells")
SeuratObject_filtered$main <- factor(SeuratObject_filtered$main, levels = hclust_data$labels[hclust_data$order])

genes_to_plot <- c("PRM3", "SYCP2", "TNP1", "TNP2", "TYR", 
                   "GAD2", "LRRTM4", "NRG3", "NCAM1", "MOBP", "AQP4", "PDGFRA", 
                   "C1QA", "JCHAIN", "PTPRC", "CD3E", "P2RY12", "CD74", "STAB2", "MRC1", "CD163", "IL7R", "THEMIS", 
                   "TTN", "MYBPC3", "MYH7", "MYH1", 
                   "PECAM1", "RELN", "PTPRB", "FLT1", "FABP4", "ERG", "DMBT1", 
				   "ECM2", "SCN7A", "FSHR" ,"EBF2", "MYH11", "DCN", "SPRY1", 
				   "SLC12A1", "SLC13A3", "EGFL7", "AEBP1", "KRT17", "DSG3", "KLF5", "KRT6A", "TCHH", "KRT14", 
				   "MUCL1", "GATA3", "ALB", "CFAP299", "SFTPC", "AGER", "ADIPOQ"
				   )


DotPlot <- DotPlot(SeuratObject_filtered, features = genes_to_plot, group.by = "main", assay='RNA' )
DotPlot(SeuratObject_filtered, features = genes_to_plot, group.by = "main", assay='RNA') +
  coord_flip() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend("Percent Expression")) +
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))