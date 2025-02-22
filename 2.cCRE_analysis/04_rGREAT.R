##Peak set process
#awk 'FNR > 1' Res/*_leiqiong.txt > leiqiong.txt
#awk '$5 < 0.05 {print $1}' leiqiong.txt > filtered_leiqiong.txt

#awk 'FNR > 1' Res/*_mongolian.txt > mongolian.txt
#awk '$5 < 0.05 {print $1}' mongolian.txt > filtered_mongolian.txt
############################################
#############near gene
library(data.table)
library(rtracklayer)
library(rGREAT)
library(ggplot2)
library(tidyverse)
gr_leiqiong <- fread("/home/Jingliangliang/SC/Peak_df/DA_CRE/DAcCRE_celltype/Immune-cell-Microglia_leiqiong.txt",header=F)


# gr_mongolian <- fread("/home/Jingliangliang/SC/Peak_df/rGREAT/filtered_mongolian.txt",header=F)
# gr_mongolian <- unique(gr_mongolian$V1)

# gr_leiqiong_only <- setdiff(gr_leiqiong, gr_mongolian)
# gr_mongolian_only <- setdiff(gr_mongolian, gr_leiqiong)
# gr_overlap <- intersect(gr_leiqiong,gr_mongolian)

# #createGRanges
# createGRanges <- function(region_strings) {
  # if (!is.vector(region_strings) || !is.character(region_strings)) {
    # stop("Input must be a character vector.")
  # }
  # gr_parts <- tstrsplit(region_strings, "_")
  # seqnames <- gr_parts[[1]]
  # start <- as.numeric(gr_parts[[2]])
  # end <- as.numeric(gr_parts[[3]])
  # gr <- GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end))
  # return(gr)
# }

# gr_leiqiong <- createGRanges(gr_leiqiong_only)

gr_leiqiong <- GRanges(seqnames = gr_leiqiong$V1, ranges = IRanges(start = gr_leiqiong$V2, end = gr_leiqiong$V3))
res_leiqiong <- great(gr_leiqiong, "GO:BP", biomart_dataset = "btaurus_gene_ensembl")
tb_leiqiong <- getEnrichmentTable(res_leiqiong)
tb_leiqiong2 <- getEnrichmentTable(res_leiqiong) %>%
	tibble::as_tibble() %>%
	dplyr::select(id,genome_fraction,observed_region_hits,p_value) %>%
	dplyr::mutate(
		updated_pvalue = exp(pbinom(
			q = observed_region_hits - 1,
			size = attr(res_leiqiong,which = "n_total"),
			prob = genome_fraction,
			lower.tail = F,
			log.p = T
		))
	)
tb_leiqiong$p_value<-tb_leiqiong2$updated_pvalue
tb_leiqiong$log<-0
for(i in which(tb_leiqiong$p_value==0)){
tb_leiqiong$log[i]<-pbinom(q = tb_leiqiong$observed_region_hits[i] - 1,size = attr(res_leiqiong,which = "n_total"),prob = tb_leiqiong$genome_fraction[i],lower.tail = F,log.p = T) / log(10)
}
tb_leiqiong$p_adjust<-p.adjust(tb_leiqiong$p_value,"BH")

######################################
gr_mongolian <- fread("/home/Jingliangliang/SC/Peak_df/DA_CRE/DAcCRE_celltype/Immune-cell-Microglia_mongolian.txt",header=F)

gr_mongolian <- GRanges(seqnames = gr_mongolian$V1, ranges = IRanges(start = gr_mongolian$V2, end = gr_mongolian$V3))

res_mongolian <- great(gr_mongolian, "GO:BP", biomart_dataset = "btaurus_gene_ensembl")
tb_mongolian <- getEnrichmentTable(res_mongolian)
tb_mongolian2 <- getEnrichmentTable(res_mongolian) %>%
	tibble::as_tibble() %>%
	dplyr::select(id,genome_fraction,observed_region_hits,p_value) %>%
	dplyr::mutate(
		updated_pvalue = exp(pbinom(
			q = observed_region_hits - 1,
			size = attr(res_mongolian,which = "n_total"),
			prob = genome_fraction,
			lower.tail = F,
			log.p = T
		))
	)
tb_mongolian$p_value<-tb_mongolian2$updated_pvalue
tb_mongolian$log<-0
for(i in which(tb_mongolian$p_value==0)){
tb_mongolian$log[i]<-pbinom(q = tb_mongolian$observed_region_hits[i] - 1,size = attr(res_mongolian,which = "n_total"),prob = tb_mongolian$genome_fraction[i],lower.tail = F,log.p = T) / log(10)
}
tb_mongolian$p_adjust<-p.adjust(tb_mongolian$p_value,"BH")
######################################
gr <- createGRanges(gr_overlap)
res <- great(gr, "GO:BP", biomart_dataset = "btaurus_gene_ensembl")
tb <- getEnrichmentTable(res)
tb2 <- getEnrichmentTable(res) %>%
	tibble::as_tibble() %>%
	dplyr::select(id,genome_fraction,observed_region_hits,p_value) %>%
	dplyr::mutate(
		updated_pvalue = exp(pbinom(
			q = observed_region_hits - 1,
			size = attr(res,which = "n_total"),
			prob = genome_fraction,
			lower.tail = F,
			log.p = T
		))
	)
tb$p_value[which(tb$p_value==0)]<-tb2$updated_pvalue[tb2$p_value==0]
tb$p_adjust<-p.adjust(tb$p_value,"BH")
#save
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "tb_overlap")
addWorksheet(wb, "tb_leiqiong")
addWorksheet(wb, "tb_mongolian")
writeData(wb, sheet="tb_overlap", tb)
writeData(wb, sheet="tb_leiqiong", tb_leiqiong)
writeData(wb, sheet="tb_mongolian", tb_mongolian)
saveWorkbook(wb, file = "multiple_sheets.xlsx", overwrite = TRUE)
#plot
tb2<-tb[(tb$fold_enrichment>2)&(tb$p_adjust<2.2e-16),]
tb3<-tb2[order(tb2$p_adjust),]
tb4<-tb3[1:15,]
tb4 %>%
    mutate(description = fct_reorder(description, -log10(p_adjust))) %>%
    ggplot( aes(x=description, y=-log10(p_adjust))) +
    geom_bar(stat="identity", fill="#f68060") +
    coord_flip() +
    xlab("") +
    theme_bw()+theme(panel.grid = element_blank())
	
tb_leiqiong2<-tb_leiqiong[(tb_leiqiong$fold_enrichment>2)&(tb_leiqiong$p_adjust<2.2e-16),]
tb_leiqiong3<-tb_leiqiong2[order(tb_leiqiong2$p_adjust),]
tb_leiqiong4<-tb_leiqiong3[1:15,]
tb_leiqiong4 %>%
    mutate(description = fct_reorder(description, -log10(p_adjust))) %>%
    ggplot( aes(x=description, y=-log10(p_adjust))) +
    geom_bar(stat="identity", fill="#f68060") +
    coord_flip() +
    xlab("") +
    theme_bw()+theme(panel.grid = element_blank())
	
tb_mongolian2<-tb_mongolian[(tb_mongolian$fold_enrichment>2)&(tb_mongolian$p_adjust<2.2e-16),]
tb_mongolian3<-tb_mongolian2[order(tb_mongolian2$p_adjust),]
tb_mongolian4<-tb_mongolian3[1:15,]
tb_mongolian4 %>%
    mutate(description = fct_reorder(description, -log10(p_adjust))) %>%
    ggplot( aes(x=description, y=-log10(p_adjust))) +
    geom_bar(stat="identity", fill="#f68060") +
    coord_flip() +
    xlab("") +
    theme_bw()+theme(panel.grid = element_blank())