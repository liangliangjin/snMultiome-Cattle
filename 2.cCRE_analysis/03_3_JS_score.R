suppressPackageStartupMessages({
library(Matrix)
library(dplyr)
library(data.table)
library(tidyr)
library(purrr)
})

datafr <- readRDS("datafr.rds") # CRE x Cell
#proj <- readRDS("combine_all_new_fig1.rds")
#meta <- as.data.frame(proj@cellColData)
meta <- readRDS("meta.rds")
meta <- as.data.frame(meta)

common_cells <- intersect(colnames(datafr), rownames(meta))
datafr <- datafr[, common_cells]
meta <- meta[common_cells, , drop = FALSE]

calculate_raw_proportions <- function(datafr, meta) {
  celltypes <- unique(meta$main)
  raw_prop_matrix <- matrix(
    0,
    nrow = nrow(datafr),
    ncol = length(celltypes),
    dimnames = list(rownames(datafr), celltypes)
  )
  for (ct in celltypes) {
    ct_cells <- rownames(meta)[meta$main == ct]
    if (length(ct_cells) > 0) {
      raw_prop_matrix[, ct] <- Matrix::rowMeans(datafr[, ct_cells, drop = FALSE] > 0)
    }
  }
  raw_prop_matrix
}

raw_prop_matrix <- calculate_raw_proportions(datafr, meta)

bias_by_main <- meta %>%
  dplyr::group_by(main) %>%
  dplyr::summarise(
    TSS_median = median(TSSEnrichment, na.rm = TRUE),
    log10_nFrags_median = median(log10(nFrags), na.rm = TRUE),
    .groups = "drop"
  )

bias_vec_tss <- bias_by_main$TSS_median
names(bias_vec_tss) <- bias_by_main$main

bias_vec_frag <- bias_by_main$log10_nFrags_median
names(bias_vec_frag) <- bias_by_main$main

bias_vec_tss <- bias_vec_tss[colnames(raw_prop_matrix)]
bias_vec_frag <- bias_vec_frag[colnames(raw_prop_matrix)]

tss_norm <- mean(log10(bias_vec_tss), na.rm = TRUE) / log10(bias_vec_tss)
frag_norm <- mean(bias_vec_frag, na.rm = TRUE) / bias_vec_frag

tss_norm <- pmin(tss_norm, 2)
frag_norm <- pmin(frag_norm, 2)

combined_norm <- tss_norm * frag_norm
combined_norm <- pmin(combined_norm, 3)

prop_norm_matrix <- t(
  t(raw_prop_matrix) * combined_norm[colnames(raw_prop_matrix)]
)

makeprobsvec <- function(p) {
  if (sum(p) <= 0) return(rep(0, length(p)))
  p / sum(p)
}

shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <= 0) return(0)
  p.norm <- p[p > 0] / sum(p)
  -sum(log2(p.norm) * p.norm)
}

JSdistVec <- function(p, q) {
  JSdiv <- shannon.entropy((p + q) / 2) -
    0.5 * (shannon.entropy(p) + shannon.entropy(q))
  JSdiv <- max(JSdiv, 0)
  sqrt(JSdiv)
}

specificity_score <- function(normpropmat) {
  marker_gene_specificities <- lapply(seq_len(ncol(normpropmat)), function(cell_type_i) {
    perfect_specificity <- rep(0, ncol(normpropmat))
    perfect_specificity[cell_type_i] <- 1
    apply(normpropmat, 1, function(x) {
      if (sum(x) > 0) {
        1 - JSdistVec(makeprobsvec(x), perfect_specificity)
      } else {
        0
      }
    })
  })
  out <- do.call(cbind, marker_gene_specificities)
  rownames(out) <- rownames(normpropmat)
  colnames(out) <- colnames(normpropmat)
  out
}

JS_score <- specificity_score(prop_norm_matrix)

colnames(JS_score) <- gsub(" ", "-", colnames(JS_score))
colnames(JS_score) <- gsub("_", "-", colnames(JS_score))
saveRDS(JS_score, "JS_score_binomial_biasadj.rds")
