# Dotplot of Activator/Repressor
library(dplyr)
library(ggplot2)
eregulon <- read.csv('~/SC/scenicplus/outs/eRegulon_direct.tsv', sep = '\t')
eregulon <- eregulon %>%
    mutate(
        class = case_when(
            rho_TF2G >= 0.05 ~ "Activator",
            rho_TF2G <= -0.05 ~ "Repressor",
            TRUE ~ "Unclear"
        ),
        R2G_group = cut(importance_R2G,
                        breaks = quantile(importance_R2G, probs = c(0, 0.25, 0.5, 0.75, 1)),
                        labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
                        include.lowest = TRUE)
    )

ggplot(eregulon, aes(x = rho_TF2G, y = importance_TF2G, 
                          color = class, size = R2G_group)) +
    geom_point(alpha = 0.5) +
    geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", color = "gray", size = 0.5) +
    scale_color_manual(values = c("Activator" = "#ff7f0e", "Repressor" = "#1f77b4", "Unclear" = "#7f7f7f")) +
    scale_size_manual(
        name = "importance_R2G",
        values = c(1, 1.5, 2, 2.5)
    ) +
    labs(
        title = "Region-TF-gene relationships",
        x = "rho_TF2G",
        y = "importance_TF2G",
        color = "Regulation Type"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "right",
        legend.box = "vertical",
        panel.border = element_rect(color = "black", fill = NA, size = 0.8)
    )


# Heatmap_dotplot
library(tidyverse)
library(RColorBrewer)

rel_data <- read.csv("~/SC/scenicplus/scenicplus_TF_expression_AUC_RSS_long_table.csv") %>%
  filter(Repressor_activator == "Activators") %>%
  filter(grepl("_\\+/\\+_", Regulon))

cell_type_order <- rel_data %>% pull(Cell_type) %>% unique()
cell_type_order <- cell_type_order[cell_type_order != "Unknown cells"]
cell_labels <- sub("^[^-]+-", "", cell_type_order)

rel_data <- rel_data %>%
  filter(Cell_type %in% cell_type_order) %>%
  mutate(
    Cell_type = factor(Cell_type, levels = cell_type_order),
    Cell_label = sub("^[^-]+-", "", as.character(Cell_type))
  )

top_regs <- rel_data %>%
  group_by(Cell_type) %>%
  slice_max(RSS_scale, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  pull(Regulon) %>%
  unique()

rel_data <- rel_data %>%
  filter(Regulon %in% top_regs) %>%
  mutate(ct_rank = as.numeric(Cell_type))

ord <- rel_data %>%
  mutate(expr_pos = pmax(Expression_scale, 0)) %>%
  group_by(Regulon) %>%
  summarise(
    bucket = ct_rank[which.max(RSS_scale)][1],
    max_rss = max(RSS_scale, na.rm = TRUE),
    expr_center = ifelse(sum(expr_pos) > 0, weighted.mean(ct_rank, expr_pos), bucket),
    .groups = "drop"
  ) %>%
  mutate(diag_dist = abs(expr_center - bucket)) %>%
  arrange(bucket, desc(max_rss), diag_dist)

regulon_order <- rev(ord$Regulon)

rel_data <- rel_data %>%
  mutate(
    Regulon = factor(Regulon, levels = regulon_order),
    TF_label = sub("_.*$", "", as.character(Regulon)),
    x = as.numeric(Cell_type),
    y = as.numeric(Regulon)
  )

tf_labels <- rel_data %>%
  distinct(Regulon, TF_label) %>%
  arrange(as.numeric(Regulon)) %>%
  pull(TF_label)

tri_expr <- rel_data %>%
  transmute(id = row_number(), fill = Expression_scale,
            x1 = x - .5, y1 = y - .5,
            x2 = x + .5, y2 = y - .5,
            x3 = x - .5, y3 = y + .5) %>%
  pivot_longer(-c(id, fill), names_to = c(".value", "pt"), names_pattern = "(.)(.)")

tri_auc <- rel_data %>%
  transmute(id = row_number(), fill = AUC_scale,
            x1 = x + .5, y1 = y + .5,
            x2 = x + .5, y2 = y - .5,
            x3 = x - .5, y3 = y + .5) %>%
  pivot_longer(-c(id, fill), names_to = c(".value", "pt"), names_pattern = "(.)(.)")

p <- ggplot() +
  geom_polygon(data = tri_expr, aes(x, y, group = id, fill = fill), linewidth = 0) +
  geom_polygon(data = tri_auc, aes(x, y, group = id, fill = fill), linewidth = 0) +
  geom_point(data = rel_data, aes(x, y, size = RSS_scale), shape = 23, fill = "white", color = "grey85") +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limits = c(-2.5, 2.5),
                       name = "Scaled TF expression / AUC") +
  scale_radius(range = c(0.8, 3.8), limits = c(0, 1), name = "Scaled RSS") +
  scale_x_continuous(breaks = seq_along(cell_type_order), labels = cell_labels, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(regulon_order), labels = tf_labels, expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 5),
    panel.grid = element_blank()
  )

ggsave("~/SC/scenicplus/TF_expression_AUC_RSS_heatmap_dotplot_top15.pdf", p, width = 15, height = 20)
p
