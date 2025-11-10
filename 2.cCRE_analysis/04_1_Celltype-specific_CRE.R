# Set1: Celltype-specific CRE
#Filter edgeR test files
library(data.table)
library(purrr)
meta <- readRDS("meta.rds")

file_list <- list.files("./Peak_df/DA_CRE_Onlycelltypes/Res", pattern = "\\.txt$", full.names = TRUE)
edgeR_list <- map(file_list, function(file) {
  data <- fread(file)
  setnames(data, 1, "CRE")
  data[padj < 0.01 & logFC < 0]
})
names(edgeR_list) <- gsub("\\.txt$", "", basename(file_list))

#Filter binomialff test
rds_files <- list.files(path = "./results", pattern = "^diff_test_cluster_\\d+\\.rds$", full.names = TRUE)
daps_list <- map(rds_files, function(f) {
  df <- readRDS(f)
  if ("qval" %in% colnames(df)) {
    df[df$qval < 0.01, , drop = FALSE]
  } else {
    data.frame()
  }
})
names(daps_list) <- gsub("^.*diff_test_cluster_|\\.rds$", "", rds_files)
daps_list <- daps_list[order(as.numeric(names(daps_list)))]
names(daps_list) <- gsub(" ", "-", unique(meta$main)[as.numeric(names(daps_list))])

#
common_celltypes <- intersect(names(edgeR_list), names(daps_list))
intersection_results <- list()
for (celltype in common_celltypes) {
  edgeR_df <- edgeR_list[[celltype]]
  binomial_df <- daps_list[[celltype]]
  if (!"CRE" %in% colnames(binomial_df) && nrow(binomial_df) > 0) {
    binomial_df$CRE <- rownames(binomial_df)
  }
  if (!inherits(edgeR_df, "data.table")) edgeR_df <- as.data.table(edgeR_df)
  if (!inherits(binomial_df, "data.table")) binomial_df <- as.data.table(binomial_df)
  common_cre <- intersect(edgeR_df$CRE, binomial_df$CRE)
  if (length(common_cre) > 0) {
    edgeR_sub <- edgeR_df[CRE %in% common_cre] 
    binomial_sub <- binomial_df[CRE %in% common_cre]
    result <- data.table(
      CRE = common_cre,
      celltype = celltype,
      edgeR_logFC = edgeR_sub$logFC,
      edgeR_logCPM = edgeR_sub$logCPM,
      edgeR_PValue = edgeR_sub$PValue,
      edgeR_padj = edgeR_sub$padj,
      binomial_pval = binomial_sub$pval,
      binomial_beta = binomial_sub$beta,
      binomial_qval = binomial_sub$qval
    )
    intersection_results[[celltype]] <- result
	rownames(intersection_results[[celltype]]) <- result$CRE
    cat("Celltype", celltype, ":", length(common_cre), "intersection CREs\n")
  }
}
final_intersection <- rbindlist(intersection_results, fill = TRUE)

saveRDS(final_intersection, "final_intersection.rds")

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
specificity_scorer <- function(normpropmat){
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

marker_specificities_out <- specificity_scorer(propmatnormbylogdepth)

markerdup <- marker_specificities_out^2
markering <- markerdup * propmatnormbylogdepth
rownames(markering) <- rownames(propmatnormbylogdepth)
colnames(markering) <- gsub(" ", "-", colnames(propmatnormbylogdepth))

markerlists <- markerlistmaker(markering, intersection_results)

nulldist = log10(as.numeric(unlist(markerlists$nonmarkerlist)))
truedist = log10(as.numeric(unlist(markerlists$markerlist)))
specificity_threshold <- quantile(nulldist, probs = 1-0.01, na.rm = TRUE)

wbmat <- data.frame()
markerlist <- markerlists$markerlist
for(i in 1:length(markerlist)){
if(is.null(markerlist[[i]])) next
specsites <- markerlist[[i]][which(log10(markerlist[[i]]) > specificity_threshold)]
currmatout <- matrix(names(specsites))
currcluster <- names(markerlist)[i]
currmatout = cbind(currmatout,intersection_results[[match(currcluster,names(intersection_results))]][match(names(specsites),rownames(intersection_results[[match(currcluster,names(intersection_results))]])),])
currmatout <- cbind(currmatout, specsites)
currmatout <- currmatout[order(currmatout[,"specsites"], decreasing=T), ]
clusterer <- strsplit(currcluster, "-(?=[^-]+$)", perl = TRUE)[[1]]
clusternow <- gsub("clusters_","",clusterer[1])
subnow <- gsub("cluster_","",clusterer[2])
assigns <- cbind(rep(clusternow,nrow(currmatout)),rep(subnow,nrow(currmatout)))
currmatforwbmat <- cbind(currmatout, assigns)
rownames(currmatforwbmat) <- NULL
colnames(currmatforwbmat)[(ncol(currmatforwbmat)-1):ncol(currmatforwbmat)] <- c("group","subset_cluster")
wbmat <- rbind(wbmat, currmatforwbmat)
}
wbmat[["V1"]] <- NULL 
set1 <- wbmat[wbmat$edgeR_logFC < -1 & wbmat$edgeR_logCPM > 1]
CRE_split <- do.call(rbind, strsplit(as.character(set1$CRE), "_"))
set1_gr <- GRanges(
  seqnames = CRE_split[,1],
  ranges = IRanges(start = as.numeric(CRE_split[,2]), end = as.numeric(CRE_split[,3])),
  strand = "*",
  set1
)
set1_gr_anno <- fastAnnoPeaks(peaks=set1_gr)

#Filter wilcoxon test form ArchRProj
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "main",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markersArchR <- getMarkers(markersPeaks)
names(markersArchR) <- gsub(" ", "-", names(markersArchR))

set1_final_anno <- GRanges()
for (ct in intersect(names(markersArchR), unique(set1_gr_anno$celltype))) {
  df_archr <- as.data.frame(markersArchR[[ct]])
  if (nrow(df_archr) == 0) next
  gr_archr <- GRanges(seqnames = df_archr$seqnames, ranges = IRanges(df_archr$start, df_archr$end))
  gr_custom <- set1_gr_anno[set1_gr_anno$celltype == ct]
  overlap <- subsetByOverlaps(gr_custom, gr_archr, minoverlap = 1)
  if (length(overlap) > 0) {
    overlap$source <- "Both"
    set1_final_anno <- c(set1_final_anno, overlap)
    cat(ct, ":", length(overlap), "common peaks\n")
  }
}

write.table(set1_final_anno, "set1_celltype_specificity_anno.txt", row.names=F, quote=F, sep="\t")
saveRDS(set1_final_anno, "set1_celltype_specificity_anno.rds")
