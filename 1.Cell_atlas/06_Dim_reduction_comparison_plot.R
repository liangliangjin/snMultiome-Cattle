library(ggplot2)
library(ArchR)

archr_umap <- as.matrix(proj@embeddings$UMAP$df[, 1:2])
archr_clusters <- as.numeric(factor(proj@cellColData$Clusters))
archr_clusters_anno <- as.numeric(factor(proj@cellColData$main_old))

wnn_umap <- as.matrix(proj@embeddings$UMAP_wnn$df[, 1:2])
wnn_clusters <- as.numeric(factor(proj@cellColData$WNN_Clusters))
wnn_clusters_anno <- as.numeric(factor(proj@cellColData$main))

ch_archr <- calinhara(archr_umap, archr_clusters, cn = length(unique(archr_clusters)))
ch_wnn   <- calinhara(wnn_umap, wnn_clusters, cn = length(unique(wnn_clusters)))
ch_archr_anno <- calinhara(archr_umap, archr_clusters_anno, cn = length(unique(archr_clusters_anno)))
ch_wnn_anno   <- calinhara(wnn_umap, wnn_clusters_anno, cn = length(unique(wnn_clusters_anno)))

df <- data.frame(
  Method = factor(rep(c("ArchR", "WNN"), each = 2), levels = c("ArchR", "WNN")),
  Group  = factor(rep(c("Cluster", "Annotation"), times = 2), levels = c("Cluster", "Annotation")),
  CH     = c(ch_archr, ch_archr_anno, ch_wnn, ch_wnn_anno)
)

ggplot(df, aes(x = Group, y = CH, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), color = "black", width = 0.65) +
  geom_text(aes(label = round(CH, 0)), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c(
    "ArchR" = "#8B73BA", 
    "WNN"   = "#D6BC36" 
  )) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  ) +
  ylab("Calinski–Harabasz Index") +
  ggtitle("CH Index Comparison: ArchR vs Seurat WNN") +
  theme(plot.title = element_text(hjust = 0.5))



library(lisi)
library(ggplot2)
library(dplyr)
library(tidyr)
lsi_meta <- data.frame(
  breed   = as.factor(proj@cellColData$Clusters3),
  sample  = as.factor(proj@cellColData$Sample)
)

archr_lisi_c <- compute_lisi(archr_umap, lsi_meta, c('breed', 'sample'))
wnn_lisi_c   <- compute_lisi(wnn_umap,lsi_meta, c('breed', 'sample'))


lisi_df <- bind_rows(as.data.frame(archr_lisi_c) %>% mutate(method = "ArchR-UMAP", cell = rownames(.)), as.data.frame(wnn_lisi_c) %>%
  mutate(method = "WNN-UMAP", cell = rownames(.))) %>% mutate(log_LISI = log(`breed`))

lisi_long <- lisi_df %>%
  pivot_longer(
    cols = c(breed, sample),
    names_to = "lisi_type",
    values_to = "lisi_value"
  ) %>%
  mutate(
    log_lisi_value = log(lisi_value),
    lisi_type = case_when(
      lisi_type == "breed" ~ "iLISI_breed",
      lisi_type == "sample" ~ "iLISI_sample"
    )
  )

compare_df <- lisi_long %>%
  group_by(lisi_type) %>%
  summarise(
    p_value = wilcox.test(
      log_lisi_value[method == "ArchR-UMAP"], 
      log_lisi_value[method == "WNN-UMAP"], 
      paired = TRUE
    )$p.value,
    median_archr = median(log_lisi_value[method == "ArchR-UMAP"]),
    median_wnn = median(log_lisi_value[method == "WNN-UMAP"]),
    median_diff = median_wnn - median_archr,
    .groups = "drop"
  ) %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

ggplot(lisi_long, aes(x = lisi_type, y = log_lisi_value, fill = method)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_text(data = compare_df, 
            aes(x = lisi_type, y = max(lisi_long$log_lisi_value) * 1.05, label = label),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  theme_classic() +
  ylab("log(LISI)") +
  xlab("") +
  scale_fill_manual(values = c("ArchR-UMAP" = "#1f77b4", "WNN-UMAP" = "#ff7f0e")) +
  scale_x_discrete(labels = c("iLISI_breed" = "iLISI\n(Breed)", 
                             "iLISI_sample" = "iLISI\n(Sample)")) +
  ggtitle("Comparison of LISI: ArchR vs WNN") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "top", panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )
