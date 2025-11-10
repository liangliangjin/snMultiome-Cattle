library(ArchR)
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(org.Bt.eg.db)
addArchRThreads(threads = 16)


proj<-readRDS("combine_all_new_fig1_end.rds")

#AddCellColData
proj$Clusters3 <- ifelse(grepl("2", proj$Sample), "Mongolian", "Leiqiong")
proj <- addCellColData(ArchRProj = proj, data = paste0(proj@cellColData$main,"_",proj@cellColData$Clusters3), name = "celltype_cattle", cells = getCellNames(proj), force = TRUE)
#Callpeak using macs2
pathToMacs2 <- findMacs2()
proj <- addGroupCoverages(ArchRProj = proj, maxCells =2000, groupBy = "celltype_cattle",force = TRUE)
#Iterative overlapping peak, faCount your.fa
proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "celltype_cattle", pathToMacs2 = pathToMacs2, genomeSize = 2.7e+09, cutOff = 0.05, force = TRUE)
proj <- addPeakMatrix(proj, force = TRUE)

#Save Peak Matrix
Filter_by_proportion <- function(x, proportion){
    Binary_df = (x > 0)
    peak_sum = rowSums(Binary_df)
    x[(peak_sum >= proportion * dim(x)[2]), ]
}
Peak_mtx = getMatrixFromProject(proj, useMatrix = 'PeakMatrix')
temp_df = Peak_mtx@assays@data$PeakMatrix
peak_meta = as.data.frame(Peak_mtx@rowRanges)
peaks = tidyr::unite(peak_meta, "peaks", seqnames, start, end)
rownames(temp_df) = peaks$peaks
xsp_Count = Filter_by_proportion(temp_df, proportion = 0)
dim(xsp_Count)
Count_df = summary(xsp_Count)
colnames(Count_df) = c('Peaks', 'Cells', 'Counts')
Count_cells = colnames(xsp_Count)
Count_peaks = rownames(xsp_Count)
write.table(x = Count_df, file="./Peak_df/ArchR_cattle/Count_df.txt", sep='\t', quote = FALSE, row.names = FALSE)
write.table(x = Count_cells, file="./Peak_df/ArchR_cattle/Count_Cells.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = Count_peaks, file="./Peak_df/ArchR_cattle/Count_Peaks.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

#Get pesudobulk accessibility
meta = proj@cellColData
celltype_list = gsub("/| ","-",unique(proj$main))
write.table(celltype_list, file="celltype_list", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
Samples = unique(proj$Sample)

peaknames <- read_delim("./Peak_df/ArchR_cattle/Count_Peaks.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
cellnames <- read_delim("./Peak_df/ArchR_cattle/Count_Cells.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
datafr = sparseMatrix(i = Count_df$Peaks, 
                     j = Count_df$Cells, 
                     x = Count_df$Counts,
                     dims = c(length(peaknames[[1]]), length(cellnames[[1]])),
                     dimnames = list(peaknames[[1]], cellnames[[1]]))
#saveRDS(datafr,"datafr.rds")

#In cross-tissue analysis, filter out the small number of cells that confuse annotations
tissue_prop_by_celltype <- as.data.frame(proj@cellColData) %>%
  group_by(main) %>%
  count(tissue, name = "count") %>%
  mutate(prop = count / sum(count) * 100, keep = (prop >= 1) & (count >= 10))

tissue_celltype_reserve <- gsub(" ", "-", paste(tissue_prop_by_celltype$main, tissue_prop_by_celltype$tissue, sep = "_x_")[tissue_prop_by_celltype$keep])

#Samples and celltypes
dir.create("./Peak_df/DA_CRE", showWarnings = FALSE, recursive = TRUE)
for (temp_sample in Samples){
    for (temp_celltype in celltype_list){
		if(paste(temp_celltype, sub("\\d$", "", temp_sample), sep = "_x_") %in% tissue_celltype_reserve){
        output_file = paste0("./Peak_df/DA_CRE/" ,temp_sample, '_', temp_celltype, '.txt')
        Used_meta = meta[gsub("/| ", "-", meta$main)==temp_celltype, ]
        Used_meta = Used_meta[Used_meta$Sample==temp_sample, ]
        Used_cell = row.names(Used_meta)
		if(length(Used_cell) > 1){
        Bulk_df = rowSums(datafr[, Used_cell])
        write.table(Bulk_df, output_file, sep='\t', quote=FALSE)}
    }}
}

#Only celltypes
dir.create("./Peak_df/DA_CRE_Onlycelltypes", showWarnings = FALSE, recursive = TRUE)
for (temp_celltype in celltype_list){
	output_file = paste0("./Peak_df/DA_CRE_Onlycelltypes/",temp_celltype, '.txt')
	Used_meta = meta[gsub("/| ","-",meta$main)==temp_celltype, ]
	Used_cell = row.names(Used_meta)
	if(length(Used_cell) > 1){
	Bulk_df = rowSums(datafr[, Used_cell])
	write.table(Bulk_df, output_file, sep='\t', quote=FALSE)}
}


  
