########################ArchR to Seurat/Signac Format
library(ArchR)
library(Seurat)
library(SingleCellExperiment)
library(Signac)
library(AnnotationHub)
library(ArchRtoSignac)

ah <- AnnotationHub()
#query(ah, "Bos_taurus.ARS-UCD1.2.105.gtf")
bos_edb <- ah[["AH109513"]]
annotations <- GetGRangesFromEnsDb(ensdb = bos_edb)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))

##If there is no peakmatrix, add it first.
#1. pathToMacs2 <- findMacs2()
#2. proj <- addGroupCoverages(ArchRProj = proj, maxCells =2000, groupBy = "Sample", force = TRUE)
#3. proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "Sample", pathToMacs2 = pathToMacs2, genomeSize = 2.7e+09, cutOff = 0.05, force = TRUE)
pkm <- getPeakMatrix(proj)

#annotations<-readRDS("~/SC/annotations_ARS.rds")
fragments_dirs <- as.list(paste0("~/SC/",unique(proj@cellColData$Sample),"/outs/"))
SeuratObject <- ArchR2Signac(
  ArchRProject = proj,
  refversion = "bosTau9",
  samples = unique(proj@cellColData$Sample),
  fragments_dir = fragments_dirs,
  pm = pkm,
  fragments_fromcellranger = "NO",
  fragments_file_extension = 'atac_fragments.tsv.gz', 
  annotation = annotations
)

##add rna
getRNAMatrix <- function(
  ArchRProject,
  SeuratObject
){
  print("In Progress:")
  print("Get RNA Matrix From ArchRProject")
  GeneExpression_matrix <- ArchR::getMatrixFromProject(ArchRProject, useMatrix='GeneExpressionMatrix')
  gem <- assays(GeneExpression_matrix)$GeneExpressionMatrix # peaks sparse martix
  print("get Gene Features From ArchRProject")
  GeneFeatures <-getFeatures(
    ArchRProj = ArchRProject,
    useMatrix = "GeneExpressionMatrix",
    select = NULL,
    ignoreCase = TRUE
  )
  colnames(gem) <- gsub("#", "_", colnames(gem))
  ix <- match(colnames(SeuratObject), colnames(gem))
  gem <- gem[,ix]
  print("Saving Gene Features From ArchRProject into Gene Score Matrix")
  rownames(gem) <- GeneFeatures
  print("Return RNA Matrix")
  gem
}
gem <- getRNAMatrix(ArchRProject = proj, SeuratObject = SeuratObject)
SeuratObject[['RNA']] <- CreateAssayObject(counts = gem)
SeuratObject <- addTwoDimRed(
  ArchRProject = proj,
  SeuratObject = SeuratObject,
  addUMAPs = "UMAP",
  reducedDims1 = "LSI_Combined",
  reducedDims2 = "Harmony" 
)
SeuratObject <- Seurat::FindVariableFeatures(SeuratObject,nfeatures=2000)
SeuratObject <- Seurat::FindVariableFeatures(SeuratObject,assay='RNA',nfeatures=2000)
saveRDS(SeuratObject,file="SeuratObject.rds")


########################Seurat to CDS Format
library(cicero)
library(monocle3)
library(ggplot2)
library(patchwork)
library(stringr)
library(SeuratWrappers)
library(Signac)
library(Seurat)
SeuratObject <- readRDS("SeuratObject.rds")
names(SeuratObject@assays)[1]<-'ATAC'
DefaultAssay(SeuratObject) <- 'ATAC'
DefaultAssay(object = SeuratObject[['umap']]) <- 'ATAC'
DefaultAssay(object = SeuratObject[['harmony']]) <- 'ATAC'
cds <- as.cell_data_set(x=SeuratObject)
cattle_cicero <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$UMAP)
saveRDS(cattle_cicero, file = "cattle_cicero.rds")
message("save cds_cicero done!")
