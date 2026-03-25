####Generate JS score matrix
datafr <- readRDS("datafr.rds") #CRE x Cell
meta = as.data.frame(proj@cellColData)
common_cells <- intersect(colnames(datafr), rownames(meta))
datafr <- datafr[, common_cells]
meta <- meta[common_cells, , drop = FALSE]


calculate_raw_proportions <- function(datafr, meta){
  celltypes <- unique(meta$main)
  raw_prop_matrix <- matrix(0, nrow=nrow(datafr), ncol=length(celltypes))
  rownames(raw_prop_matrix) <- rownames(datafr)
  colnames(raw_prop_matrix) <- celltypes
  for(ct in celltypes){
    ct_cells <- rownames(meta)[meta$main == ct]
    if(length(ct_cells) > 0){
      ct_data <- datafr[, ct_cells, drop=FALSE]
      raw_prop_matrix[, ct] <- Matrix::rowMeans(ct_data > 0)
    }
  }
  return(raw_prop_matrix)
}
datafr_bin <- (datafr > 0) * 1
raw_prop_matrix <- calculate_raw_proportions(datafr_bin, meta)

calculate_median_sites_per_cluster <- function(datafr, meta) {
  clusters <- unique(meta$main)
  median_sites <- numeric(length(clusters))
  names(median_sites) <- clusters
  for(cluster in clusters) {
    cluster_cells <- rownames(meta)[meta$main == cluster]
    cluster_data <- datafr[, cluster_cells, drop = FALSE]
    sites_per_cell <- Matrix::colSums(cluster_data > 0)
    median_sites[cluster] <- median(sites_per_cell)
  }
  return(median_sites)
}

median_accessibility <- calculate_median_sites_per_cluster(datafr, meta)
depthmat <- data.frame(
  cluster = names(median_accessibility),
  depth = as.numeric(median_accessibility)
)

logdepthnorm <- mean(log10(depthmat$depth)) / log10(depthmat$depth)
names(logdepthnorm) <- depthmat$cluster

cell_counts <- table(meta$main)
cell_counts <- cell_counts[names(median_accessibility)]
weights <- log10(cell_counts[colnames(raw_prop_matrix)] + 1) / log10(max(cell_counts) + 1)

propmatnormbylogdepth <- t(t(raw_prop_matrix) * logdepthnorm[colnames(raw_prop_matrix)] * as.numeric(weights))

makeprobsvec <- function(p){
  phat <- p / sum(p)
  phat[is.na(phat)] <- 0
  phat
}
#entropy
shannon.entropy <- function(p){
  if(min(p)<0 || sum(p)<=0) return(Inf)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}
#calculate JS specificity score
JSdistVec <- function(p,q){
  JSdiv <- shannon.entropy((p+q)/2) - 0.5*(shannon.entropy(p)+shannon.entropy(q))
  JSdiv[is.infinite(JSdiv)] <- 1
  JSdiv[JSdiv<0] <- 0
  sqrt(JSdiv)
}
specificity_score <- function(normpropmat){
  marker_gene_specificities <- lapply(1:ncol(normpropmat), function(cell_type_i){
    perfect_specificity <- rep(0.0,ncol(normpropmat))
    perfect_specificity[cell_type_i] <- 1.0
    apply(normpropmat, 1, function(x){
      if(sum(x) > 0) 1 - JSdistVec(makeprobsvec(x), perfect_specificity)
      else 0
    })
  })
  do.call(cbind, marker_gene_specificities)
}

markerlistmaker <- function(markering, daps){
  markerlist <- list(); betamarkerlist <- list(); nonmarkerlist <- list(); markerback <- list()
  for(i in 1:ncol(markering)){
  if(colnames(markering)[i] %in% names(daps)){
    commoners <- intersect(rownames(daps[[match(colnames(markering)[i], names(daps))]]), rownames(markering))
    markerlist[[i]] <- markering[match(commoners, rownames(markering)), i]
	names(markerlist[[i]]) <- rownames(markering)[match(commoners, rownames(markering))]
    betamarkerlist[[i]] <- daps[[match(colnames(markering)[i], names(daps))]][match(commoners, rownames(daps[[match(colnames(markering)[i], names(daps))]])), c("binomial_beta")]
    nonmarkerlist[[i]] <- markering[-match(commoners, rownames(markering)), i]
    names(markerlist)[i] <- colnames(markering)[i]
    names(betamarkerlist)[i] <- colnames(markering)[i]
    names(nonmarkerlist)[i] <- colnames(markering)[i]
  }}
  markerback$markerlist <- markerlist
  markerback$betamarkerlist <- betamarkerlist
  markerback$nonmarkerlist <- nonmarkerlist
  markerback
}

JS_score <- specificity_score(propmatnormbylogdepth)
colnames(JS_score) <- gsub(" ", "_", colnames(propmatnormbylogdepth))
saveRDS(JS_score, file = "JS_score.rds")
