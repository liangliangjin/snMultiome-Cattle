##
df <- as.data.frame(getCellColData(proj))

###cellnumber
df_counts <- df %>%
  group_by(Clusters4, Clusters3) %>%
  summarise(count = n())
total_counts <- aggregate(count ~ Clusters4, data = df_counts, sum)
df_counts <- merge(df_counts, total_counts, by = "Clusters4", suffixes = c("", "_total"))
df_counts$proportion <- df_counts$count / df_counts$count_total


p_cellnumber<- ggplot(df_counts, aes(x = Clusters4, y = count, fill = Clusters3)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label = count_total), 
              position = position_stack(vjust = 0.5), 
              color = "white") + 
    labs(
        title = "Bar Plot of Samples by Clusters3",
        x = "Clusters4",
        y = "count"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


p_cellnumber<- ggplot(df_counts, aes(x = Clusters4, y = count, fill = Clusters3)) + 
    geom_bar(stat = "identity") +
    geom_text(aes(label = count_total), 
              position = position_stack(vjust = 0.5), 
              color = "black", size = 4.5) +
    labs(
        title = "Bar Plot of Samples by Clusters3",
        x = "Clusters4",
        y = "Count"
    ) +
    theme_light(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    ) +
    scale_fill_brewer(palette = "Set3")

##TSS
mean_values <- df %>%
  group_by(Clusters4, Clusters3) %>%
  summarise(mean_TSS_Enrichment = median(TSSEnrichment, na.rm = TRUE))
p_tss<-ggplot(mean_values, aes(x = Clusters4, y = mean_TSS_Enrichment)) +
  geom_point(size = 3, color = "blue") +
  geom_line(group = 1, color = "blue") +
  labs(
    title = "Median TSS Enrichment by Clusters3",
    x = "Clusters3",
    y = "Median TSS Enrichment"
  ) +
  theme_minimal() +  coord_cartesian(ylim = c(0, NA))+ facet_wrap(~Clusters3)

nFrags_values <- df %>%
  group_by(Clusters4, Clusters3) %>%
  summarise(nFrags = median(nFrags, na.rm = TRUE))
p_nFrags<-ggplot(nFrags_values, aes(x = Clusters4, y = log10(nFrags),fill = Clusters3)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Median nFrags by Clusters3",
    x = "Clusters3",
    y = "Median TSS Enrichment"
  ) +
  theme_minimal() +  coord_cartesian(ylim = c(0, NA))+ facet_wrap(~Clusters3)
p_cellnumber+p_tss+p_nFrags



p_qc_rna<-df %>% 
  	ggplot(aes(x=Gex_nUMI, y=Gex_nGenes, color=Gex_MitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500,linetype = "dashed") +
  	geom_hline(yintercept = 250,linetype = "dashed") +
  	facet_wrap(~Clusters3)


								


# Cell Cycle
seurat_phase <- lapply(seurat_objects, function(obj) {
    NormalizeData(obj)
	ScaleData(obj)
	})
	
for (i in seq_along(seurat_phase)) {
  s.genes <- intersect(cc.genes$s.genes, rownames(seurat_phase[[i]]))
  g2m.genes <- intersect(cc.genes$g2m.genes, rownames(seurat_phase[[i]]))
  seurat_phase[[i]] <- CellCycleScoring(seurat_phase[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
}
combined_meta <- NULL
for (i in seq_along(seurat_phase)) {
  sce <- seurat_phase[[i]]
  meta_data_with_id <- cbind(sce@meta.data, Source = paste("Object", i))
  combined_meta <- rbind(combined_meta, meta_data_with_id)
}
print(combined_meta)
all(proj$cellNames%in%rownames(combined_meta))
combined_meta2<-combined_meta[rownames(combined_meta)%in%proj$cellNames,]
combined_meta2$celltype<-proj$main[match(rownames(combined_meta2),proj$cellNames)]

pS <- ggplot(combined_meta, aes(x = Phase, y = S.Score, fill = Phase)) +
  geom_violin(trim = FALSE) +
  labs(title = "S.Score Distribution", x = "Phase", y = "S.Score") +
  theme_minimal()

pG2M <- ggplot(combined_meta, aes(x = Phase, y = G2M.Score, fill = Phase)) +
  geom_violin(trim = FALSE) +
  labs(title = "G2M.Score Distribution", x = "Phase", y = "G2M.Score") +
  theme_minimal()

