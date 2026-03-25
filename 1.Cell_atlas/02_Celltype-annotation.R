library(Seurat)
library(dplyr)
library(ArchR)
addArchRThreads(threads = 16)

SeuratObject <- readRDS("SeuratObject_wnn.rds")
proj <- readRDS("combine_all_after_filter.rds")
proj <- addImputeWeights(proj, reducedDims = "Harmony")

dir.create("./0.marker/WNN_Clusters", showWarnings = FALSE)

markers_all <- FindAllMarkers(
  object = SeuratObject,
  group.by = "WNN_Clusters",
  assay = "SCT",
  layer = "data",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25,
)
markers_significant <- markers_all %>% filter(p_val_adj < 0.01 & avg_log2FC > 0.5) %>% 
  mutate(
    cluster = as.numeric(as.character(cluster)),
    gene = as.character(gene)
  ) %>% arrange(cluster, desc(avg_log2FC))
write.table(markers_significant, "markers_significant.txt", quote=F, row.names = F, sep = "\t")

# All clusters of marker genes were examined and annotated manually in conjunction with the literature
Annotation <- read.table("Annotation_result.txt", header = T, sep = "\t")

# Use large language models (LLM, From clusterProfiler::interpret) to analyze enrichment results for cross-validation of manual annotations.
library(clusterProfiler)
library(org.Bt.eg.db)

enrich_results <- list()
for(cl in unique(markers_significant$cluster)) {
	message(cl)
	enrich_results[[as.character(cl)]] <- NULL
	markers <- markers_significant %>% 
		filter(cluster == cl) %>% arrange(desc(avg_log2FC))
	geneList <- markers$avg_log2FC
	names(geneList) <- markers$gene
	geneList <- sort(geneList, decreasing = TRUE)
	if(length(geneList) >= 5){
		eG <- gseGO(gene = geneList,
			keyType = 'SYMBOL',
			OrgDb = org.Bt.eg.db,
			ont = "BP", pvalueCutoff = 0.05,
			pAdjustMethod = "BH")
		if(!is.null(eG) && nrow(eG) > 0) {
		  enrich_results[[as.character(cl)]] <- eG
		}
	}
}
saveRDS(enrich_results, file = "enrich_results.rds")

#Obtain the API and key, see https://github.com/YuLab-SMU/fanyi
#set_translate_option(appid = '', key = '', source = "dsk")
res <- list()
for (cl in names(enrich_results)) {
  tryCatch({
	ssample <- table(SeuratObject$Clusters4[SeuratObject$WNN_Clusters==cl]) %>% prop.table %>% `*`(100) %>% .[. > 5] %>% sort(decreasing = TRUE)
	tissue_string <- paste0(names(ssample), "(", round(ssample, 0), "%)", collapse = "; ")
	top50genes <- paste(head(markers_significant$gene[markers_significant$cluster == cl], 50), collapse = ", ")
    interp_result <- interpret(
      enrich_results[[cl]],
	  context = paste0("Multi-tissue single-cell sequencing of cattle, and the main tissue sources are as follows: ", tissue_string, " the top highly variable genes are as follows: ", top50genes),
      prior = Annotation$annotation[Annotation$celltype==cl],
      model = "deepseek-chat",
      task = "cell_type"
    )
    res[[cl]] <- interp_result
  }, error = function(e) {
    message("Error in cluster ", cl, ": ", e$message)
    res[[cl]] <- list(error = e$message)
  })
  message(paste0(cl, ": The cell type manually annotated as [",Annotation$annotation[Annotation$celltype==cl],"] is inferred to be [", res[[cl]]$cell_type, "] with a confidence level of [", res[[cl]]$confidence, "]."))
}

#Rename celltypes
labelOld <- as.character(unique(proj$WNN_Clusters))
Annotation <- read.table("Annotation_result_final.txt", header = T, sep = "\t") #Adjusted annotations
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