# cd /home/Jingliangliang/SCARlink
# conda activate scarlink-env
library(ArchR)
library(Seurat)
library(rhdf5)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(cowplot)
library(DoubletFinder)
library(clustree)
organ <- c("spleen0","spleen2")
set_resolutions <- c(0.2,0.05) # Adjust to the appropriate resolution

f <- function(a,b){
	gsub(setdiff(strsplit(b,"")[[1]],strsplit(a,"")[[1]]),"",b)
	}
	
load(paste0("/home/Jingliangliang/SC/",f(organ[1],organ[2]),"-combine.RData"))

for(i in organ){
sample<-get(i)
sample <- subset(sample, subset = nFeature_RNA > 200 & percent.mt < 5)
sample <- SCTransform(sample, verbose = F)
sample <- FindVariableFeatures(sample,selection.method = "vst", features = rownames(sample))
sample <- ScaleData(sample, features = rownames(sample))
#PCA降维
sample <- RunPCA(sample, npcs = 30, verbose = FALSE)
#细胞聚类
sample <- sample %>%
  RunUMAP(reduction = "pca", dims = 1:30,n.neighbors=30,min.dist=0.1) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0) %>%
  identity()
saveRDS(sample, file = paste0("scrna_",i,".rds"))
rm(sample)
}
# for(i in 1:2){
# scatac.object<-readRDS(paste0("/home/Jingliangliang/SC/",organ[i],".rds"))
# resolutions<-as.numeric(gsub("LSI_ATAC_","",names(scatac.object@reducedDims[grep("ATAC", names(scatac.object@reducedDims))])))
# print(resolutions)
# }
for(i in 1:length(organ)){
scatac.object<-readRDS(paste0("/home/Jingliangliang/SC/",organ[i],".rds"))
set.seed(1)
scatac.object <- addTileMatrix(input=scatac.object, binarize=FALSE, tileSize = 500, force=TRUE)
print(names(scatac.object@reducedDims))
resolutions<-as.numeric(gsub("LSI_ATAC_","",names(scatac.object@reducedDims[grep("ATAC", names(scatac.object@reducedDims))])))
set_resolution<-intersect(set_resolutions[i],resolutions)
print(set_resolution)
scatac.object <-addIterativeLSI(ArchRProj = scatac.object,
    clusterParams = list(
      resolution = set_resolution, 
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix", 
    depthCol = "nFrags",force = TRUE,
    name = "IterativeLSI")
scatac.object$celltype<-gsub("[)]","_",gsub("[(]","_",gsub(" ", "_",scatac.object$Clusters2)))
saveArchRProject(ArchRProj = scatac.object, outputDirectory = paste0("/home/Jingliangliang/2025-SCARlink/",organ[i],"_archr"), load = FALSE)
}