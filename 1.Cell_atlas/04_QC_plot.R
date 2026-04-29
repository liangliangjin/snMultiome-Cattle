#Figure S1b tissue origin
#Load packages
suppressPackageStartupMessages({
library(ggplot2)
library(tidyr)
library(patchwork)
library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
})
#celltype <- strsplit(as.character(proj@cellColData$main), "-")
#cellgroup <- sapply(celltype, function(x) x[1])
#addCellColData(proj, data = cellgroup, name = "Clusters_GO", cells = getCellNames(proj))
df <- as.data.frame(getCellColData(proj))
metadata <- getCellColData(proj, select = c("Clusters_GO", "tissue")) %>%
  as.data.frame()

prop_df <- metadata %>%
  dplyr::group_by(Clusters_GO, tissue) %>% filter(Clusters_GO != "Unknown cells") %>%
  dplyr::summarise(count = n(), .groups = "drop_last") %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  dplyr::ungroup() %>% mutate(Clusters_GO = forcats::fct_rev(factor(Clusters_GO)))

ggplot(prop_df, aes(x = proportion, y = Clusters_GO, fill = tissue)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_x_continuous(
    name = "Proportion",
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    y = NULL,
    title = "Tissue Origin Proportion"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.box = "horizontal",
    axis.text.y = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(15, 20, 15, 20),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    panel.background = element_rect(fill = "white")
  )

#Figure S1c cattle origin
df <- as.data.frame(getCellColData(proj))
df_counts <- df %>%
  dplyr::group_by(Clusters_GO, Clusters3) %>%
  dplyr::summarise(count = n())  %>% filter(Clusters_GO != "Unknown cells") %>%
  mutate(Clusters_GO = forcats::fct_rev(factor(Clusters_GO)))
total_counts <- aggregate(count ~ Clusters_GO, data = df_counts, sum)
df_counts <- merge(df_counts, total_counts, by = "Clusters_GO", suffixes = c("", "_total"))
df_counts$proportion <- df_counts$count / df_counts$count_total


ggplot(df_counts, aes(x = proportion, y = Clusters_GO, fill = Clusters3)) +
  geom_bar(stat = "identity", width = 0.8) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.box = "horizontal",
    axis.text.y = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(15, 20, 15, 20),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    panel.background = element_rect(fill = "white")
  ) + scale_fill_manual(values = c("Leiqiong" = "#656565", "Mongolian" = "#CCCCCA"))



# Figure S1d
metadata <- getCellColData(
    proj,
    select = c("Sample", "Gex_nUMI", "Gex_nGenes", "nFrags")
) %>% as.data.frame()

mycolor <- paletteDiscrete(values=unique(proj$Sample))

plotvin <- function(data, x_var = "Sample", y_var, title = NULL, 
                   custom_colors = NULL) {
    data[[x_var]] <- factor(data[[x_var]], levels = sort(unique(data[[x_var]])))
    p <- ggplot(data, aes_string(x = x_var, y = y_var, fill = x_var)) +
        geom_violin(scale = "width", trim = TRUE) +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 12),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = "black", size = 0.5),
            axis.line.y = element_line(color = "black", size = 0.5),
            axis.line.x = element_blank(),
            
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 0.8)
        ) +
        labs(
            x = NULL,
            y = y_var,
            title = ifelse(is.null(title), y_var, title)
        )
    if (!is.null(custom_colors)) {
        p <- p + scale_fill_manual(values = custom_colors)
    } else {
        p <- p + scale_fill_viridis_d(option = "D", begin = 0.1, end = 0.9)
    }
    
    return(p)
}
p1 <- plotvin(metadata, y_var = "Gex_nUMI", title = "Gex_nUMI", custom_colors = mycolor)
p2 <- plotvin(metadata, y_var = "Gex_nGenes", title = "Gex_nGenes", custom_colors = mycolor)
p3 <- plotvin(metadata, y_var = "nFrags", title = "nFrags", custom_colors = mycolor)
p1 / p2 / p3

# Figure S1e,f TSS & frag
p_TSS <- plotTSSEnrichment(ArchRProj = proj, groupBy = "Sample") + theme(legend.position = "none")
p_frag <- plotFragmentSizes(ArchRProj = proj, groupBy = "Sample")
(p_TSS + theme(legend.position = "none")) / p_frag


# Figure S1g
# load SeuratObject
SeuratObject <- readRDS("SeuratObject_wnn.rds")
meta <- SeuratObject@meta.data
rna <- GetAssayData(SeuratObject, "RNA", "counts")
atac <- GetAssayData(SeuratObject, "peaks", "counts")
ind1 <- sort(unique(meta$Clusters3))[1]
ind2 <- sort(unique(meta$Clusters3))[2]
tissues <- unique(meta$Clusters4)
rna_plots <- list()
atac_plots <- list()
for(tis in tissues) {
  cells1 <- colnames(SeuratObject)[meta$Clusters4 == tis & meta$Clusters3 == ind1]
  cells2 <- colnames(SeuratObject)[meta$Clusters4 == tis & meta$Clusters3 == ind2]
  if(length(cells1) == 0 | length(cells2) == 0) next
  # snRNA
  x <- log10(1 + rowMeans(rna[, cells1]))
  y <- log10(1 + rowMeans(rna[, cells2]))
  keep <- x > 0.1 | y > 0.1
  df <- data.frame(x = x[keep], y = y[keep])
  if(nrow(df) > 10) {
    ct <- cor.test(df$x, df$y)
    rna_plots[[tis]] <- ggplot(df, aes(x, y)) +
      geom_point(alpha = 0.2, size = 0.5) +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      annotate("text", x = 0.1, y = 2.9, 
               label = paste0(round(ct$estimate, 3), "\n", 
                              if(ct$p.value < 0.0001) "p<0.0001" else paste0("p=", round(ct$p.value, 4))), 
               hjust = 0, size = 2.5) +
      labs(title = tis) + xlim(0, 3) + ylim(0, 3) + 
      theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA),
                              axis.line = element_blank(), axis.ticks = element_blank(),
                              axis.text = element_blank(), axis.title = element_blank(),
                              plot.title = element_text(hjust = 0.5, size = 9)) +
      coord_fixed()
  }
  # snATAC
  x <- log10(1 + Matrix::rowMeans(atac[, cells1]))
  y <- log10(1 + Matrix::rowMeans(atac[, cells2]))
  keep <- x > 0.1 | y > 0.1
  df <- data.frame(x = x[keep], y = y[keep])
  if(nrow(df) > 10) {
    ct <- cor.test(df$x, df$y)
    atac_plots[[tis]] <- ggplot(df, aes(x, y)) +
      geom_point(alpha = 0.2, size = 0.5, color = "darkgreen") +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      annotate("text", x = 0.05, y = 0.95, 
               label = paste0(round(ct$estimate, 3), "\n", 
                              if(ct$p.value < 0.0001) "p<0.0001" else paste0("p=", round(ct$p.value, 4))), 
               hjust = 0, size = 2.5) +
      labs(title = tis) + xlim(0, 1) + ylim(0, 1) + 
      theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA),
                              axis.line = element_blank(), axis.ticks = element_blank(),
                              axis.text = element_blank(), axis.title = element_blank(),
                              plot.title = element_text(hjust = 0.5, size = 9)) +
      coord_fixed()
  }
}
common <- intersect(names(rna_plots), names(atac_plots))
ordered <- sort(common)
plots <- c()
for(tis in ordered) plots <- c(plots, rna_plots[tis], atac_plots[tis])
wrap_plots(plots, ncol = length(ordered), byrow = FALSE)



# Figure S1g corr_plot
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


								
com_tissues <- intersect(unique(proj$tissue[proj$Clusters3=="Mongolian"]), unique(proj$tissue[proj$Clusters3=="Leiqiong"]))
df$log_umi <- log10(df$Gex_nUMI + 1)

plot_results <- list()
cor_summary <- data.frame()
for(tissue in com_tissues){
    cat("\nProcessing:", tissue, "...\n")
    result <- plot_breed_correlation(tissue, df)
    if(!is.null(result)){
        plot_results[[tissue]] <- result
        cor_summary <- rbind(cor_summary, data.frame(
            Tissue = tissue,
            Correlation = result$correlation,
            R_squared = result$r_squared,
            Mongolian_cells = result$n_mongolian,
            Leiqiong_cells = result$n_leiqiong,
            Cells_used = result$n_used
        ))
        print(result$plot)
    }
}

# Figure S1h Cell Cycle
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
combined_meta2 <- combined_meta[rownames(combined_meta)%in%proj$cellNames,]
combined_meta2$celltype <- proj$main[match(rownames(combined_meta2),proj$cellNames)]

pS <- ggplot(combined_meta, aes(x = Phase, y = S.Score, fill = Phase)) +
  geom_violin(trim = FALSE) +
  labs(title = "S.Score Distribution", x = "Phase", y = "S.Score") +
  theme_minimal()

pG2M <- ggplot(combined_meta, aes(x = Phase, y = G2M.Score, fill = Phase)) +
  geom_violin(trim = FALSE) +
  labs(title = "G2M.Score Distribution", x = "Phase", y = "G2M.Score") +
  theme_minimal()


# Figure S1i
# load SeuratObject
SeuratObject <- readRDS("SeuratObject_wnn.rds")
meta <- SeuratObject@meta.data
rna <- GetAssayData(SeuratObject, "RNA", "data")
atac <- GetAssayData(SeuratObject, "peaks", "data")
ind1 <- sort(unique(meta$Clusters3))[1]
ind2 <- sort(unique(meta$Clusters3))[2]
tissues <- unique(meta$Clusters4)
rna_plots <- list()
atac_plots <- list()
for(tis in tissues) {
  cells1 <- colnames(SeuratObject)[meta$Clusters4 == tis & meta$Clusters3 == ind1]
  cells2 <- colnames(SeuratObject)[meta$Clusters4 == tis & meta$Clusters3 == ind2]
  if(length(cells1) == 0 | length(cells2) == 0) next
  # snRNA
  x <- log10(1 + rowMeans(rna[, cells1]))
  y <- log10(1 + rowMeans(rna[, cells2]))
  keep <- x > 0.1 | y > 0.1
  df <- data.frame(x = x[keep], y = y[keep])
  if(nrow(df) > 10) {
    ct <- cor.test(df$x, df$y)
    rna_plots[[tis]] <- ggplot(df, aes(x, y)) +
      geom_point(alpha = 0.2, size = 0.5) +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      annotate("text", x = 0.1, y = 2.9, 
               label = paste0(round(ct$estimate, 2), "\n", 
                              if(ct$p.value < 0.0001) "p<0.0001" else paste0("p=", round(ct$p.value, 4))), 
               hjust = 0, size = 2.5) +
      labs(title = tis) + xlim(0, 3) + ylim(0, 3) + 
      theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA),
                              axis.line = element_blank(), axis.ticks = element_blank(),
                              axis.text = element_blank(), axis.title = element_blank(),
                              plot.title = element_text(hjust = 0.5, size = 9)) +
      coord_fixed()
  }
  # snATAC
  x <- log10(1 + Matrix::rowMeans(atac[, cells1]))
  y <- log10(1 + Matrix::rowMeans(atac[, cells2]))
  keep <- x > 0.1 | y > 0.1
  df <- data.frame(x = x[keep], y = y[keep])
  if(nrow(df) > 10) {
    ct <- cor.test(df$x, df$y)
    atac_plots[[tis]] <- ggplot(df, aes(x, y)) +
      geom_point(alpha = 0.2, size = 0.5, color = "darkgreen") +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      annotate("text", x = 0.05, y = 0.95, 
               label = paste0(round(ct$estimate, 2), "\n", 
                              if(ct$p.value < 0.0001) "p<0.0001" else paste0("p=", round(ct$p.value, 4))), 
               hjust = 0, size = 2.5) +
      labs(title = tis) + xlim(0, 1) + ylim(0, 1) + 
      theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA),
                              axis.line = element_blank(), axis.ticks = element_blank(),
                              axis.text = element_blank(), axis.title = element_blank(),
                              plot.title = element_text(hjust = 0.5, size = 9)) +
      coord_fixed()
  }
}
common <- intersect(names(rna_plots), names(atac_plots))
ordered <- sort(common)
plots <- c()
for(tis in ordered) plots <- c(plots, rna_plots[tis], atac_plots[tis])
wrap_plots(plots, ncol = length(ordered), byrow = FALSE)