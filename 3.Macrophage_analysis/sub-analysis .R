#Macrophage across tissues
# Adipose	Macrophage	C33 
# Brain	Microglia cell C36
# Heart	Macrophage	C38 C23
# Kidney	Mesangial cell	C34
# Liver	Kupffer cell	C27
# Lung	Macrophage	C33
# Muscle	Macrophage	C38
# Ovary	Macrophage	C33
# Rumen	Macrophage	C23	C33
# Skin	Macrophage	C33
# Spleen	Macrophage	C56
# Testis	Macrophage	↓C56 C33
## c("C23","C27","C33","C34","C36","C38","C55")
sum(proj$Clusters%in%c("C23","C27","C33","C34","C36","C38","C55"))
proj2<-proj[proj$Clusters%in%c("C23","C27","C33","C34","C36","C38","C55")]


markerTest <- getMarkerFeatures(
  ArchRProj = proj_m, 
  useMatrix = "GeneExpressionMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Macrophage_hn",
  bgdGroups = "Macrophage_mg"
)
markerList <- getMarkers(markerTest, cutOff = "FDR <= 0.05 & Log2FC >= 1")
intersect(markerList@listData$Macrophage_hn$name,m1marker)
intersect(markerList@listData$Macrophage_hn$name,m2marker)

markerTest2 <- getMarkerFeatures(
  ArchRProj = proj_m, 
  useMatrix = "GeneExpressionMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Macrophage_mg",
  bgdGroups = "Macrophage_hn"
)
markerList2 <- getMarkers(markerTest2, cutOff = "FDR <= 0.05 & Log2FC >= 1")
intersect(markerList2@listData$Macrophage_mg$name,m1marker)
intersect(markerList2@listData$Macrophage_mg$name,m2marker)

pma <- plotMarkers(seMarker = markerTest, name = "Macrophage_hn", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pv <- plotMarkers(seMarker = markerTest, name = "Macrophage_hn", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pma+pv



load("macrophage.RData")
#Refer to https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html
#1-0.library packages
library(miloR)
library(Seurat)
library(SingleCellExperiment)
library(BiocParallel)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(ggrepel)
library(reshape2)
library(scater)
library(scran)
#1-1.get the SeuratObjec format: see cCRE_analysis/00_Data_format_conversion.R
#SeuratObject <- readRDS("SeuratObject_mar.rds")
#DefaultAssay(SeuratObject) <- "RNA"
#SeuratObject <- ScaleData(SeuratObject)
#SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures(object = SeuratObject))
#ElbowPlot(SeuratObject, ndims = 50)

#1-2.SeuratObjec to SingleCellExperiment
Milo <- as.SingleCellExperiment(SeuratObject)
#Add umap dimension reduction
reducedDims(Milo) <- list(
    pca = SeuratObject[["pca"]]@cell.embeddings,
    umap = SeuratObject[["umap"]]@cell.embeddings
)

#1-3.Create a Milo object
cattle_milo <- Milo(Milo)

#1-4.Construct KNN graph
#The Nhood corresponding to k should be the value of the number of samples * 5, set prop and d to the recommended values
#plotNhoodSizeHist(cattle_milo)
#set k = 20
cattle_milo <- buildGraph(cattle_milo, k = 20, d = 30, reduced.dim = "pca")

#1-5.MakeNhoods
cattle_milo <- makeNhoods(cattle_milo, prop = 0.1, k = 20, d=30, refined = TRUE, reduced_dims = "pca")

#1-6.Counting cells in neighbourhoods
cattle_milo <- countCells(cattle_milo, meta.data = as.data.frame(colData(cattle_milo)), sample="Sample")

#1-7.Defining experimental design
#"Clusters3" represents cattle breed information and "Clusters4" represents tissue information
cattle_design <- data.frame(colData(cattle_milo))[,c("Sample", "Clusters3", "Clusters4")]
cattle_design <- distinct(cattle_design) 
rownames(cattle_design) <- cattle_design$Sample 

#1-8.Computing neighbourhood connectivity
cattle_milo <- calcNhoodDistance(cattle_milo, d=30, reduced.dim = "pca")

#1-9.Testing
#cattle
model.contrasts <- c('Clusters3Hainan - Clusters3Mongolian')
da_results <- testNhoods(cattle_milo, design=~ 0 + Clusters3, design.df=cattle_design, reduced.dim="pca", fdr.weighting="graph-overlap", model.contrasts = model.contrasts)
da_results %>%                                           
  arrange(SpatialFDR) %>%
  head() 
  
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)


cattle_milo <- buildNhoodGraph(cattle_milo)
#Plot single-cell UMAP                       
umap_pl <- plotReducedDim(cattle_milo, dimred = "umap", colour_by="Clusters3", text_by = "Clusters", text_size = 8, point_size=0.5) +guides(fill="none")

nh_graph_pl <- plotNhoodGraphDA(cattle_milo, da_results, layout="umap",alpha=0.1) 

p<-umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
  
da_results <- annotateNhoods(cattle_milo, da_results, coldata_col = "Clusters")
head(da_results)
da_results$celltype <- ifelse(da_results$Clusters_fraction < 0.7, "Mixed", da_results$Clusters)

plotDAbeeswarm(da_results, group.by = "celltype")




set.seed(42)
cattle_milo <- buildNhoodGraph(cattle_milo)
da_results <- groupNhoods(cattle_milo, da_results, max.lfc.delta = 3, overlap=5)
plotNhoodGroups(cattle_milo, da_results, layout="umap")

plotDAbeeswarm(da_results, group.by = "NhoodGroup")


keep.rows <- rowSums(logcounts(cattle_milo)) != 0
cattle_milo <- cattle_milo[keep.rows, ]

## Find HVGs
dec <- modelGeneVar(cattle_milo)
hvgs <- getTopHVGs(dec, n=2000)
head(hvgs)

nhood_markers <- findNhoodGroupMarkers(cattle_milo, da_results, subset.row = hvgs, 
                                       aggregate.samples = TRUE, sample_col = "Sample",
                                       subset.groups = c("1")
                                       )

head(nhood_markers)