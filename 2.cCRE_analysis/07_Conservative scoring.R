###################################phastCons score
#BiocManager::install("phastCons100way.UCSC.hg38")
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(phastCons100way.UCSC.hg38)
phast <- phastCons100way.UCSC.hg38

# set <- readRDS("set.rds") 
# bos_Hg38 <- read.table("bos_0.5_Hg38.bed", sep = "",header=F)

# Writes the sequence after liftover conversion
set[set$hg==1,(ncol(set)+1):(ncol(set)+3)] <- bos_Hg38
names(set)[(ncol(set)-2):ncol(set)] <- c("hg_chr", "hg_start", "hg_end")

mam_index<-(set$hg==1&set$mm==1&set$pig==1&set$horse==1&set$sheep==1)
mam<-paste(set$hg_chr[mam_index],paste(set$hg_start[mam_index],set$hg_end[mam_index],sep="-"),sep=":")

GRanges_mam<-GRanges(mam)
mam_score<-gscores(phast, GRanges_mam)

set_mm<-set[mam_index,]
set_mm$mam_score<-mam_score$default

mam_liftover_list <- list()
for (i in 1:5) {
  mam_liftover_list[[i]] <- set_mm[set_mm$Level == i, ]
  mam_liftover_list[[i]]$n <- 1:nrow(mam_liftover_list[[i]])
}

# Density plot
plot_list <- list()
for (i in 1:5) {
  d <- data.frame(score = mam_liftover_list[[i]]$mam_score)
  plot_list[[i]] <- ggplot(d, aes(x = score)) +
    geom_density(fill = "yellow", adjust = 0.5) +
    theme(panel.background = element_rect(fill = "white"))
}
phastCons_plot <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[5]]
print(phastCons_plot)

###################################phyloP score
set_mm2<-set_mm[, c("hg_chr", "hg_start", "hg_end")]
write.table(set_mm2,file="/storage/public/home/2021060195/02.phyloP/region_file.bed",sep="\t",col.names = F,row.names = F,quote = FALSE)


# run 01.bw2bed.sh
	# Step 1: BigWig to bed file.
	# Step 2: Chromosomal resolution of phyloP file.
	# Step 3: Bedtools intersect.
	# Step 4: Add line numbers and sort.
	# Step 5: Chromosomal resolution of region file.
# run 02.cal.R
	# Calculate the phyloP score for each CRE.
##################
setwd("/storage/public/home/2021060195/02.phyloP")
library(data.table)
chr_files <- list.files(pattern = "^chr.*\\.bed$")
bed_data_list <- lapply(chr_files, function(file) {
fread(file, header = FALSE, col.names = c('chr', 'start', 'end', 'score'))
})
names(bed_data_list) <- chr_files
region_data <- fread('/storage/public/home/2021060195/02.phyloP/region_file.bed', header = FALSE, col.names = c('chr', 'start', 'end'))
library(parallel)
num_cores <- detectCores()
process_region <- function(i) {
region_chr <- region_data$chr[i]
region_start <- region_data$start[i]
region_end <- region_data$end[i]
match_chr <- match(paste0(region_chr, ".bed"), chr_files)
bed_data <- bed_data_list[[match_chr]]
scores <- bed_data[start >= region_start & end <= region_end, "score"]
return(scores)
}
region_scores <- mclapply(seq_len(nrow(region_data)), process_region, mc.cores = num_cores)
means <- sapply(region_scores, function(df) mean(df$score, na.rm = TRUE))
write.table(means,file="phyloP_means.txt",sep="\t",col.names = F,row.names = F,quote = FALSE)