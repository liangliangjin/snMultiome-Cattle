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

mam_index <- (set$hg==1&set$mm==1&set$pig==1&set$horse==1&set$sheep==1)
mam <- paste(set$hg_chr[mam_index],paste(set$hg_start[mam_index],set$hg_end[mam_index],sep="-"),sep=":")

GRanges_mam <- GRanges(mam)
mam_score <- gscores(phast, GRanges_mam)

set_mm <- set[mam_index,]
set_mm$mam_score <- mam_score$default

mam_list <- list()
for (i in 1:5) {
  mam_list[[i]] <- set_mm[set_mm$Level == i, ]
  mam_list[[i]]$n <- 1:nrow(mam_list[[i]])
}

# Density plot
plot_list <- list()
for (i in 1:5) {
  d <- data.frame(score = mam_list[[i]]$mam_score)
  plot_list[[i]] <- ggplot(d, aes(x = score)) +
    geom_density(fill = "#FADAA8", adjust = 0.5) +
    theme(panel.background = element_rect(fill = "white"))
}
phastCons_plot <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[5]]
print(phastCons_plot)

###################################phyloP score
set_mm2 <- set_mm[, c("hg_chr", "hg_start", "hg_end")]
write.table(set_mm2,file="/storage/public/home/2021060195/02.phyloP/region_file.bed",sep="\t",col.names = F,row.names = F,quote = FALSE)

# run 01.bw2bed.sh
	# Step 1: BigWig to bed file.
	# Step 2: Chromosomal resolution of phyloP file.
	# Step 3: Bedtools intersect.
	# Step 4: Add line numbers and sort.
	# Step 5: Chromosomal resolution of region file.
# run 02.cal.R
	# Calculate the phyloP score for each CRE.
		# Multitask by chromosome
##################
library(data.table)
chr_files <- list.files(pattern = "_scores.txt")

combined_data <- chr_files %>%
  map(~ fread(.x, sep = "\t", header = TRUE)) %>%
  rbindlist()
  
if (!"phyloP" %in% colnames(set_mm)) {
  set_mm$phyloP <- NA
}

set_mm$phyloP[combined_data$index] <- combined_data$mean_score

mam_list <- list()
for (i in 1:5) {
  mam_list[[i]] <- set_mm[set_mm$Level == i, ]
  mam_list[[i]]$n <- 1:nrow(mam_list[[i]])
}

# Density plot
plot_list2 <- list()
for (i in 1:5) {
  d <- data.frame(score = mam_list[[i]]$phyloP)
  plot_list2[[i]] <- ggplot(d, aes(x = score)) +
    geom_density(fill = "#ADD8E6", adjust = 0.5) +
    theme(panel.background = element_rect(fill = "white")) +
	xlim(-1, 1) + ylim(0, 1.5) +
	geom_vline(xintercept = 0, linetype = "dashed", color = "black")
}
phyloP_plot <- plot_list2[[1]] + plot_list2[[2]] + plot_list2[[3]] + plot_list2[[4]] + plot_list2[[5]]
print(phyloP_plot)