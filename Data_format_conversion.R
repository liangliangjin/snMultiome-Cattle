########################ArchR to Signac Format
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

pkm <- getPeakMatrix(proj)

#annotations<-readRDS("~/SC/annotations_ARS.rds")
fragments_dirs <- as.list(paste0("/home/Jingliangliang/SC/",unique(proj@cellColData$Sample),"/outs/"))
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
  GeneScore_matrix <- ArchR::getMatrixFromProject(ArchRProject, useMatrix='GeneExpressionMatrix')
  gsm <- assays(GeneScore_matrix)$GeneExpressionMatrix # peaks sparse martix
  print("get Gene Features From ArchRProject")
  GeneFeatures <-getFeatures(
    ArchRProj = ArchRProject,
    useMatrix = "GeneExpressionMatrix",
    select = NULL,
    ignoreCase = TRUE
  )
  colnames(gsm) <- gsub("#", "_", colnames(gsm))
  ix <- match(colnames(SeuratObject), colnames(gsm))
  gsm <- gsm[,ix]
  print("Saving Gene Features From ArchRProject into Gene Score Matrix")
  rownames(gsm) <- GeneFeatures
  print("Return RNA Matrix")
  gsm
}
gsm <- getRNAMatrix(ArchRProject = proj, SeuratObject = SeuratObject)
SeuratObject[['RNA']] <- CreateAssayObject(counts = gsm)
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
