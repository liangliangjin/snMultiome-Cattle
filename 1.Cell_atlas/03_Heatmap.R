#
library(Seurat)
library(ArchR)
library(Signac)
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(org.Bt.eg.db)
#
proj <- readRDS("all_sample.rds")

##Fig2.Peak2gene
proj <- readRDS("all_sample.rds")
proj <- addPeak2GeneLinks(ArchRProj = proj, reducedDims = "Harmony", useMatrix = "GeneExpressionMatrix")
p <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "main")


SeuratObject <- readRDS("SeuratObject_wnn.rds")
DefaultAssay(SeuratObject) <- "peaks"
linkedpeaks <- list()
SeuratObject <- RegionStats(SeuratObject, BSgenome.Btaurus.UCSC.bosTau9)
Idents(SeuratObject) <- "main"
for(i in levels(SeuratObject)){
  celltype <- subset(SeuratObject, idents = i)
  celltype <- LinkPeaks(celltype, 
                        peak.assay = "peaks", 
                        expression.assay = "RNA", 
                        distance = 500000)
  linkedpeaks[[i]] <- celltype@assays[["peaks"]]@links
}
saveRDS(linkedpeaks, "linkpeaks.rds")

link_stats <- data.frame(
  celltype = names(linkedpeaks),
  n_links = sapply(linkedpeaks, nrow)
)

ggplot(link_stats, aes(x = celltype, y = n_links)) +
  geom_col() +
  theme_classic() +
  coord_flip()


#corplot
DefaultAssay(SeuratObject) <- "SCT"
rna.avg <- AverageExpression(
    SeuratObject,
    assays = "SCT",
    features = VariableFeatures(SeuratObject),
	group.by = "main"
)
celltypes <- unique(colnames(rna.avg$RNA))
corM <- cor(as.matrix(rna.avg$SCT), method = "pearson")
ordered_indices <- corrplot::corrMatOrder(
	corM[celltypes, celltypes],
	order = "hclust",
	hclust.method = "ward.D2"
	)
ordered_celltypes <- celltypes[ordered_indices]
pdf("RNA.correaltion.pdf", width = 9, height = 9)
corrplot(
    corM[ordered_celltypes, ordered_celltypes],
    method = "square",
    type = "upper",
    tl.col = "black",
    tl.cex = 0.6,
    is.corr = F,
    col = rev(COL2("RdBu", 100)),
    order = "original",col.lim = c(-1, 1)
)
dev.off()
#
DefaultAssay(SeuratObject) <- "peaks"
SeuratObject <- RunTFIDF(SeuratObject)
SeuratObject <- FindTopFeatures(SeuratObject, min.cutoff = "q0")
SeuratObject <- RunSVD(SeuratObject)
SeuratObject.atac.markers <-
    parallel::mclapply(unique(SeuratObject$main), function(x) {
        xx <- FindMarkers(
            SeuratObject,
            ident.1 = x,
            only.pos = T,
            test.use = "LR",
            max.cells.per.ident = 300,
            latent.vars = "nCount_peaks"
        )
        return(data.frame(xx, gene = rownames(xx), cluster = x))
    }, mc.cores = 20)
peaks <- lapply(SeuratObject.atac.markers, function(x) {
    if (!is.null(nrow(x))) {
            return(x)
        }
    })
peaksID <- do.call(rbind, peaks) %>%
    group_by(cluster) %>%
	filter(p_val_adj < 0.001) %>%
    top_n(500, avg_log2FC) %>%
    pull(gene) %>%
    unique()

SeuratObject.atac.avg <- AverageExpression(SeuratObject,
    assays = "peaks",
    features = peaksID,
	slot = "data",
	group.by = "main"
)
corATAC <- cor(as.matrix(SeuratObject.atac.avg$peaks), method = "pearson")
pdf("ATAC.correaltion.pdf", width = 9, height = 9)
corrplot(
    corATAC[ordered_celltypes, ordered_celltypes],
    method = "square",
    type = "lower",
    tl.col = "black",
    tl.cex = 0.6,
    is.corr = F,
    col = rev(COL2("RdBu", 100)),
    order = "original", col.lim = c(-1, 1)
)
dev.off()
	
#GO enrichment
#Cell types that do not belong to this group are used as the background
proj$Clusters_GO <- as.character(sapply(proj$main, function(x) strsplit(x, "-")[[1]][1]))
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
    seMarker = markers, n = round(nrow(assay(markers))*0.01),
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
library(clusterProfiler)
library(AnnotationHub)
library(stringr)
library(msigdbr)
library(org.Bt.eg.db)
library(KEGG.db)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringdist)
library(viridis)
library(patchwork)
files<-list.files("./marker/", pattern = "txt", full.names = TRUE)
for(i in files){
geneid<-read.table(file=i, header=T, sep="\t")
gene<-geneid$x
eG <- enrichGO(gene = gene, 
			   keyType = 'SYMBOL',
               OrgDb = org.Bt.eg.db,
               ont="BP", 
			   pAdjustMethod = "BH",pvalueCutoff = 0.01,qvalueCutoff = 0.01,
               readable =T)
			   
ego_simp <- simplify(eG, cutoff=0.7, by="p.adjust", select_fun=min)
write.csv(ego_simp, file = paste0(strsplit(i,".txt")[[1]][1],".csv"))
}

plot_files <- list.files("./marker/", pattern = "\\.csv$", full.names = TRUE)
all_results <- list()
top_pathways <- list()
for (file in plot_files) {
    group_name <- tools::file_path_sans_ext(basename(file))
    data <- read.csv(file) %>% arrange(p.adjust)
	if (nrow(data) == 0) {
    next
	}
    all_results[[group_name]] <- data %>% mutate(group = group_name)
}

all_results_df <- bind_rows(all_results)

selected_df <- all_results_df %>%
  mutate(log_p = -log10(p.adjust + 1e-10))

exclude_keywords <- c(
  "cell adhesion",
  "cell junction",
  "cell-cell",
  "cell migration",
  "cell morphogenesis",
  "extracellular matrix",
  "metabolic process",
  "cell cycle",
  "DNA repair", 
  "ribosome",
  "translation",
  "mitochondrial", 
  "immune response",
  "apoptotic process",
  "cell motility",
  "vesicle transport",
  "actin filament-based process",
  "actin filament organization"
)
exclude_pattern <- str_to_lower(paste(exclude_keywords, collapse = "|"))

selected_df_clean <- selected_df %>%
  filter(!str_detect(str_to_lower(str_squish(Description)), exclude_pattern)) %>%
  distinct(group, Description, .keep_all = TRUE)


anchor_df <- selected_df_clean %>%
  group_by(group) %>%
  slice_max(log_p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    major_type = sub("gene_(.*?)-.*", "\\1", group)
  )

anchor_df <- anchor_df %>%
  arrange(major_type, desc(log_p))

group_levels <- unique(anchor_df$group)
go_levels <- unique(anchor_df$Description)

selected_df_final <- selected_df_clean %>% 
  mutate(
    group = factor(group, levels = group_levels),
    Description = factor(Description, levels = rev(go_levels)) 
  )

selected_df_final <- selected_df_clean %>% 
  mutate(
    group = factor(group, levels = group_levels),
    Description = factor(Description, levels = rev(go_levels))
  ) %>% filter(!is.na(Description) & Description != "<NA>" & !is.na(group)) 

  
ggplot(selected_df_final, aes(x = group, y = Description)) +
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
  
#Differences among cattle breeds 
SeuratObject$main <- proj$main
SeuratObject_filtered <- subset(SeuratObject, subset = main != "Unknown cells")
SeuratObject_filtered_LQ <- subset(SeuratObject_filtered, subset = Clusters3 == "Leiqiong")
SeuratObject_filtered_MG <- subset(SeuratObject_filtered, subset = Clusters3 == "Mongolian")

genes_to_plot <- c(
  # Endothelial
  "PECAM1", "FLT1", "CLDN5", "VWF", 
  # Epithelial
  "EPCAM", "KRT8", "SFTPC", "KLF5", 
  # Germline 
  "DDX4",
  # Immune
  "CD3E", "CD79A","LYZ", "C1QA", 
  # Muscle
  "MYH7", "MYH1", "MYH11", 
  # Nerve
  "NCAM1", "GAD2", "MOBP", "AQP4", 
  # Stromal
  "PDGFRA", "PDGFRB", "DCN" 
)

p1 <- DotPlot(SeuratObject_filtered_LQ, features = genes_to_plot, group.by = "main", assay='RNA') +
  coord_flip() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend("Percent Expression")) + scale_color_viridis(option = "plasma", direction = -1)
  
p2 <- DotPlot(SeuratObject_filtered_MG, features = genes_to_plot, group.by = "main", assay='RNA') +
  coord_flip() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend("Percent Expression")) + scale_color_viridis(option = "plasma", direction = -1)
  
combined_dotplot <- p1 / p2 +
  plot_layout(heights = c(1, 11),
              guides = "collect")