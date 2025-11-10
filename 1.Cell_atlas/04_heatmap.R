#
library(Seurat)
library(ArchR)
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(org.Bt.eg.db)
#
proj <- readRDS("all_sample.rds")

markersGS_RNA <- getMarkerFeatures(
	ArchRProj = proj,
	useMatrix = "GeneExpressionMatrix",
	groupBy = "main",
	bias = c("TSSEnrichment","log10(Gex_nUMI)"),
	testMethod = "wilcoxon"
)

heatmapRNA <- plotMarkerHeatmap(
    seMarker = markersGS_RNA,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    #labelMarkers = markerGenes,
    transpose = TRUE
)
exp_rna<-heatmapRNA@matrix

#hclust plot 
hclust_rna <- hclust(dist(as.matrix(exp_rna)),method="ward.D2")
plot(hclust_rna,hang=-1)

#bar plot 
ratio_data<-as.data.frame(table(proj$Clusters4,proj$main))
names(ratio_data) <- c("Organ","Celltype","CellNumber")
ratio_data$Celltype <- factor(ratio_data$Celltype, levels = hclust_rna2$labels[hclust_rna2$order])
bar_plot_RNA<-ggplot(ratio_data, aes(x = Celltype, y = CellNumber)) +
    geom_col(aes(fill = Organ), color = "black", position = 'fill') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid = element_blank()) +
    scale_fill_brewer(palette = "Set3")

#Enrichment
#Cell types that do not belong to this group are used as the background
proj$Clusters_GO<-as.character(sapply(proj$main, function(x) strsplit(x, "-")[[1]][1]))
for (i in unique(proj$main)){
print(paste0("start ",i,":"));
bg<-unique(proj$main[!proj$Clusters_GO==unique(proj$Clusters_GO[which(proj$main==i)])])
markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneExpressionMatrix",
  groupBy = "main",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment","log10(nFrags)","log10(Gex_nUMI)"),
  useGroups = i,
  bgdGroups = bg
)
a <- getMarkers(
    seMarker = markers,
    n = round(nrow(assay(markers))*0.01),
    returnGR = FALSE
)
gene <- a[[i]]$name
if(length(gene)>0){
write.table(gene,paste0("./marker/gene_",gsub("/","_",gsub(" ","_",i)),".txt"),col.names=T,row.names=F,sep="\t");
print(paste0("Output success!    {",i,"}    ",length(gene)))
}
else{print("error")}
rm(list=c("bg","markers","a","gene"))
}
#Annotate using the clusterProfiler package
library(AnnotationHub)
library(clusterProfiler)
library(stringr)
library(msigdbr)
library(org.Bt.eg.db)
library(KEGG.db)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringdist)
library("viridis")
files<-list.files("./marker/",pattern = "txt",full.names = TRUE)
for(i in files){
geneid<-read.table(file=i,header=T,sep="\t")
gene<-geneid$x
eG <- enrichGO(gene = gene, 
			   keyType = 'SYMBOL',
               OrgDb = org.Bt.eg.db,
               ont="BP", 
			   pAdjustMethod = "BH",pvalueCutoff = 0.1,qvalueCutoff = 0.1,
               readable =F)
write.csv(eG,file=paste0("./marker/",strsplit(i,".txt")[[1]][1],".csv"))
}

plot_files <- list.files(pattern = "\\.csv$", full.names = TRUE)
all_results <- list()
top_pathways <- list()
for (file in plot_files) {
    group_name <- tools::file_path_sans_ext(basename(file))
    data <- read.csv(file) %>% arrange(p.adjust)
	top_pathways[[group_name]] <- data %>% head(3) %>% mutate(group = group_name)
    all_results[[group_name]] <- data %>% mutate(group = group_name)
}

all_results_df <- bind_rows(all_results)
top_pathways_df <- bind_rows(top_pathways)

top_ids <- unique(top_pathways_df$ID)
selected_df <- all_results_df %>%
  filter(ID %in% top_ids) %>%
  mutate(log_p = -log10(p.adjust + 1e-10))


order <- gsub("/","_",gsub(" ","_",paste0("gene_",hclust_rna2$labels[hclust_rna2$order])))
selected_df$group <- factor(selected_df$group, levels = hclust_rna2$labels[hclust_rna2$order])

pathway_order <- selected_df %>%
  filter(!is.na(p.adjust) & p.adjust > 0) %>%
  mutate(
    logp = -log10(p.adjust),
    weight = logp * Count
  ) %>%
  group_by(Description) %>%
  summarize(
    weighted_group_pos = weighted.mean(as.numeric(group), weight + 1e-5),
    .groups = "drop"
  ) %>%
  arrange(weighted_group_pos) %>%
  pull(Description)


selected_df$Description <- factor(selected_df$Description, levels = rev(pathway_order)) 


ggplot(selected_df, aes(x = group, y = Description)) +
  geom_point(aes(size = Count, color = log_p)) +
  scale_color_viridis(option = "viridis") +
  labs(
    title = "Enrichment map",
    x = "Cell types",
    y = "GO pathway",
    color = "-log10(p.adjust)",
    size = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  )