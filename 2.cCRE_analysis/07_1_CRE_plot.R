# Set2: Breed-specific CRE
library(edgeR)
library(GenomicRanges)
output_dir <- "./Peak_df/DA_CRE/DAcCRE_breed_specific/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

leiqiong_df <- read.delim("./Peak_df/DA_CRE/Settings/leiqiong_df.txt", row.names = 1, check.names = FALSE)
mongolian_df <- read.delim("./Peak_df/DA_CRE/Settings/mongolian_df.txt", row.names = 1, check.names = FALSE)
meta <- readRDS("meta.rds") 


all_counts <- cbind(leiqiong_df, mongolian_df)
group <- factor(c(rep("leiqiong", ncol(leiqiong_df)), rep("mongolian", ncol(mongolian_df))), levels = c("leiqiong", "mongolian"))

y <- DGEList(counts = all_counts, group = group)
CPM <- cpm(y)
keep <- rowSums(CPM > 2)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y)
bcv <- 0.1
et <- exactTest(y, dispersion = bcv^2)
res <- et$table
res$padj <- p.adjust(res$PValue, method = "BH")

CPM <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
binMat <- (CPM > 5) * 1
open_prop_leiqiong <- rowSums(binMat[, group == "leiqiong"]) / sum(group == "leiqiong")
open_prop_mongolian <- rowSums(binMat[, group == "mongolian"]) / sum(group == "mongolian")

res$open_prop_leiqiong <- open_prop_leiqiong[rownames(res)]
res$open_prop_mongolian <- open_prop_mongolian[rownames(res)]
res$avg_CPM <- rowMeans(CPM)[rownames(res)]

leiqiong_specific <- res[res$logFC < -1 & res$padj < 0.01 & res$open_prop_leiqiong > 0.3 & res$avg_CPM > 5, ]
leiqiong_specific <- leiqiong_specific[order(leiqiong_specific$logFC), ]

mongolian_specific <- res[res$logFC > 1 & res$padj < 0.01 & res$open_prop_mongolian > 0.3 & res$avg_CPM > 5, ]
mongolian_specific <- mongolian_specific[order(mongolian_specific$logFC, decreasing = TRUE), ] 

parse_genomic_coordinates <- function(peak_names) {
    coords <- strsplit(peak_names, "_")
    chr_index <- which(sapply(coords[[1]], function(x) grepl("chr", x)))[1]
    chr <- sapply(coords, function(x) x[chr_index])
    start <- as.numeric(sapply(coords, function(x) x[chr_index + 1]))
    end <- as.numeric(sapply(coords, function(x) x[chr_index + 2]))
    return(data.frame(chr = chr, start = start, end = end))
}
leiqiong_df <- parse_genomic_coordinates(rownames(leiqiong_specific))
leiqiong_gr <- GRanges(
  seqnames = leiqiong_df$chr,
  ranges = IRanges(start = leiqiong_df$start, end = leiqiong_df$end),
  strand = "*",
  leiqiong_specific
)
mongolian_df <- parse_genomic_coordinates(rownames(mongolian_specific))
mongolian_gr <- GRanges(
  seqnames = mongolian_df$chr,
  ranges = IRanges(start = mongolian_df$start, end = mongolian_df$end),
  strand = "*",
  mongolian_specific
)

leiqiong_specific_anno <- fastAnnoPeaks(peaks=leiqiong_gr)
mongolian_specific_anno <- fastAnnoPeaks(peaks=mongolian_gr)

write.table(leiqiong_specific_anno, paste0(output_dir, "leiqiong_specific_anno.txt"), sep="\t", row.names = F, quote=FALSE)
write.table(mongolian_specific_anno, paste0(output_dir, "mongolian_specific_anno.txt"), sep="\t", row.names = F, quote=FALSE)

message("Finished broad CRE analysis for both breeds.")






















library(data.table)

merged_data <- fread("merged_results.txt")

celltype_groups <- list(
  nerve_cells = c('Nerve-cell-Neuron', 'Nerve-cell-Granule-cell', 'Nerve-cell-GABAergic-neuron', 
                  'Nerve-cell-Excitatory-neuron', 'Nerve-cell-Astrocyte', 
                  'Nerve-cell-Oligodendrocyte-precursor-cell', 'Nerve-cell-Oligodendrocyte', 
                  'Nerve-cell-Neural-precursor-cell', 'Nerve-cell-Schwann-cell'),
  endothelial_cells = c('Endothelial-cell-Cardiac-endothelial-cell-1', 'Endothelial-cell-Cardiac-endothelial-cell-2',
                       'Endothelial-cell-Ruminal-endothelial-cell-1', 'Endothelial-cell-Ruminal-endothelial-cell-2', 
                       'Endothelial-cell-Liver-sinusoidal-endothelial-cell', 'Endothelial-cell-Ovarian-endothelial-cell', 
                       'Endothelial-cell-Other-endothelial-cell', 'Endothelial-cell-Splenic-sinusoidal-endothelial-cell'),
  immune_cells = c('Immune-cell-Microglia', 'Immune-cell-Splenic-macrophage', 'Immune-cell-Cardiac-macrophage', 
                   'Immune-cell-NK-T-cell', 'Immune-cell-T-cell', 'Immune-cell-Macrophage', 'Immune-cell-B-cell', 
                   'Immune-cell-Muscle-macrophage', 'Immune-cell-Kupffer-cell', 'Immune-cell-Testicular-T-cell'),
  muscle_cells = c('Muscle-cell-Other-smooth-muscle-cell', 'Muscle-cell-Cardiomyocyte-1', 'Muscle-cell-Cardiomyocyte-2', 
                   'Muscle-cell-Type-II-myonuclei', 'Muscle-cell-Type-I-myonuclei', 'Muscle-cell-Perivascular-smooth-muscle-like-cell', 
                   'Muscle-cell-Ruminal-smooth-muscle-cell', 'Muscle-cell-Cardiac-smooth-muscle-cell', 
                   'Muscle-cell-Hybrid-skeletal-muscle-fibers'),
  epithelial_cells = c('Epithelial-cell-Hepatocyte', 'Epithelial-cell-Follicular-cell', 'Epithelial-cell-Alveolar-epithelial-type-II-cell', 
                      'Epithelial-cell-Alveolar-epithelial-type-I-cell', 'Epithelial-cell-Keratinocyte', 
                      'Epithelial-cell-Supporting-cell', 'Epithelial-cell-Proximal-convoluted-tubule-2', 
                      'Epithelial-cell-Proximal-convoluted-tubule-1', 'Epithelial-cell-Renal-tubular-epithelial-cell', 
                      'Epithelial-cell-Loop-of-Henle---Distal-convoluted-tubule', 'Epithelial-cell-Hair-follicle-cell', 
                      'Epithelial-cell-Ruminal-epithelial-cell'),
  stromal_cells = c('Stromal-cell-Stellate-cell', 'Stromal-cell-Ruminal-fibroblast', 'Stromal-cell-Mesenchymal-cell', 
                    'Stromal-cell-Adipocyte-precursor-cell', 'Stromal-cell-Cardiac-fibroblast', 'Stromal-cell-Mesangial-cell', 
                    'Stromal-cell-Satellite-cell', 'Stromal-cell-Splenic-fibroblast', 'Stromal-cell-Adipocyte', 
                    'Stromal-cell-Lung-pericyte', 'Stromal-cell-Myotendinous-junction', 'Stromal-cell-Dermal-fibroblast'),
  germline_cells = c('Germline-cell-Spermatid', 'Germline-cell-Spermatocyte', 'Germline-cell-Spermatogonia')
)


assign_conservation_level <- function(data, celltype_groups) {
	conservation_level <- numeric(nrow(data))
  
	pb <- txtProgressBar(min = 0, max = nrow(data), style = 3)
  
	for (i in 1:nrow(data)) {
		num_ones_in_groups <- sapply(celltype_groups, function(group) sum(data[i, ..group] == 1))
		group_ones_ratio <- sapply(celltype_groups, function(group) mean(data[i, ..group] == 1))
	
		max_ratio <- max(group_ones_ratio)
		max_group_index <- which.max(group_ones_ratio)

		#wilcox.test
		group_data <- sapply(celltype_groups, function(group) data[i, ..group, drop = FALSE])
		
		num_groups <- length(celltype_groups)
		p_values <- character(choose(num_groups, 2))
		p_values_name <- NULL
		count <- 1
		for (group1_index in 1:(num_groups - 1)) {
			for (group2_index in (group1_index + 1):num_groups) {
			group1_data <- as.numeric(group_data[[group1_index]])
			group2_data <- as.numeric(group_data[[group2_index]])
			test_result <- wilcox.test(group1_data, group2_data, exact = FALSE)
			p_values[count] <- test_result$p.value
			p_values_name[count] <- paste0(group1_index,"-",group2_index)
			count <- count + 1
			}
		}
		
		p_values_name2 <- p_values_name[p_values < 0.01]
		groups_involved <- strsplit(p_values_name2, "-")
		cross <- Reduce(intersect, groups_involved)
		
		if (sum(num_ones_in_groups) == 1) {
		conservation_level[i] <- 1  # Level 1 (Only one celltype)
		} 
		else if (max_ratio < 0.5 && sum(num_ones_in_groups) > 1 && length(cross) > 0) {
			conservation_level[i] <- 2  # Level 2 (Low conservation - not across groups)
		} 
		else if (max_ratio < 0.5 && sum(num_ones_in_groups) > 1 && length(cross) == 0) {
			conservation_level[i] <- 3  # Level 3 (Low conservation - across groups)
		} 
		else if (max_ratio >= 0.5 && length(cross) > 0) {
			conservation_level[i] <- 4  # Level 4 (High conservation - not across groups)
		} 
		else if (max_ratio >= 0.5 && length(cross) == 0) {
			conservation_level[i] <- 5  # Level 5 (High conservation - across groups)
		} 
		else {
			conservation_level[i] <- 0 
		}
		
		if (sum(num_ones_in_groups) == 0) {
		conservation_level[i] <- 0
		}
		
		setTxtProgressBar(pb, i)
		cat(sprintf("\rProgress: %d/%d", i, nrow(data)))
	}
	close(pb)
	cat("\n")
	return(conservation_level)
}

conservation_level <- assign_conservation_level(merged_data, celltype_groups)

#write.table(conservation_level, "conservation_level.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

CRE_peak<-as.data.frame(fread("/storage/public/home/2021060195/SC/liftover/peakset_783227.csv",header=T))
#conservation_level <- read.table("conservation_level.txt", sep="\t", header=F)
CRE_set<-paste(CRE_peak[,1],CRE_peak[,2],CRE_peak[,3],sep = "-")#chr1-1-10 formal string
library(tidyverse)
set<-as.data.frame(list(name=CRE_set,Level=conservation_level$V1))

setwd("/storage/public/home/2021060195/SC/liftover/50-merge-nonredundant/")
library(data.table)
for(j in c("Hg38","Mm10","SusScr11","EquCab3","GCF_016772045.2")){
  assign(paste("bos", j, sep = "_"), read.table(paste("bos","_0.5_", j,".bed", sep = ""),header=F))
  assign(paste("bos", j,"unmap", sep = "_"), read.table(paste("bos","_unmap_0.5_", j,".bed", sep = ""),header=F))
  assign(paste("bos", j,"unmap_set", sep = "_"), paste(get(paste("bos", j,"unmap", sep = "_"))[,1],get(paste("bos", j,"unmap", sep = "_"))[,2],get(paste("bos", j,"unmap", sep = "_"))[,3],sep = "-"))  
}

species_unmap_sets <- list(
  hg = bos_Hg38_unmap_set,
  mm = bos_Mm10_unmap_set,
  pig = bos_SusScr11_unmap_set,
  horse = bos_EquCab3_unmap_set,
  sheep = bos_GCF_016772045.2_unmap_set
)

for (species in names(species_unmap_sets)) {
  set[[species]] <- ifelse(CRE_set %in% species_unmap_sets[[species]], 0, 1)
}

#########################################################UpSet plot
library(UpSetR)

set_plot<-set[c("name","hg","mm","pig","horse","sheep")]

for(s in names(species_unmap_sets)) {
  set_plot[[paste("unmap", s, sep = "_")]] <- (set_plot[[s]] == 0) + 0
}

upsetplot1<-upset(set_plot, nsets = 5, mb.ratio = c(0.6, 0.4),sets=c("sheep","horse","pig","mm","hg"),keep.order = TRUE,order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE),sets.bar.color = c("#20498D", "#20498D", "#20498D", "#20498D","#20498D"),matrix.color = '#225EA8',point.size = 3.3,line.size = 0.8,shade.color = 'grey',shade.alpha = 0.2,matrix.dot.alpha = 0.7)


#################
#############################################################
TSS<-read.table("gene_TSSup1k_down1k.txt",header=F)
TSS$V1<-paste0("chr",TSS$V1)

library("bedtoolsr")#options(bedtools.path = "/path/to/")
TSS_int<-bt.intersect(a = CRE_peak, b = TSS,wa=T,wb=T)
###distal as 0 and 1 when there is overlap on the proximal promoter
#set$peaktype <- CRE_peak$peakType
set$TSS <- 0
set$TSS[TSS_int$V17] <- 1 #
##mammal
for(level in 1:5) {
  name <- paste("mam_Level", level, sep = "")
  subset_data <- set[set$hg == 1 & set$mm == 1 & set$pig == 1 & set$horse == 1 & set$sheep == 1 & set$Level == level, ]
  subset_data$n <- 1:nrow(subset_data)
  assign(name, subset_data)
}

##unique to cattle
for(level in 1:5) {
  name <- paste("Level", level, sep = "")
  subset_data <- set[set$hg == 0 & set$mm == 0 & set$pig == 0 & set$horse == 0 & set$sheep == 0 & set$Level == level, ]
  subset_data$n <- 1:nrow(subset_data)
  assign(name, subset_data)
}

##Livestock
for(level in 1:5) {
  name <- paste("livestock_Level", level, sep = "")
  subset_data <- set[set$hg == 0 & set$mm == 0 & set$pig == 1 & set$horse == 1 & set$sheep == 1 & set$Level == level, ]
  subset_data$n <- 1:nrow(subset_data)
  assign(name, subset_data)
}

##bovidae
for (level in 1:5) {
  name <- paste("bovidae_Level", level, sep = "")
  subset_data <- set[set$hg == 0 & set$mm == 0 & set$pig == 0 & set$horse == 0 & set$sheep == 1 & set$Level == level, ]
  subset_data$n <- 1:nrow(subset_data)
  assign(name, subset_data)
}


# Create canvas layout
par(mfrow = c(4, 5))
for (name in c("mam_Level","livestock_Level","bovidae_Level","Level")){
	for(level in 1:5){
		data <- get(paste0(name, level))
		pie(table(data$TSS), 
			labels = paste0(names(table(data$TSS)), " (", round(prop.table(table(data$TSS)) * 100, 2), "%)"), 
			main = paste0(name, level, " (", nrow(data), ")"), 
			col = rainbow(length(table(data$TSS))))
		}
}
# Save plot
# saveRDS(set, "set.rds")