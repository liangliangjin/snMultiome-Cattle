#The current script requires a large amount of computing resources, so we submit it separately to the computing server for each cell type.
#!/bin/bash
# #jsub_submit.sh
# cells_per_sample=120
# cores=32
# outdir="results"
# mkdir -p log $outdir
# #Loop 1:61 cluster index, you need to modify the corresponding value according to your own data.
# for i in $(seq 1 61); do
    # echo "Submitting cluster index: $i"
    # jsub -R "span[hosts=1]" -q normal -J "Cluster_$i" -n $cores \
         # -o "./log/Cluster_${i}.o" -e "./log/Cluster_${i}.e" \
         # Rscript run_monocle_glm.R \
         # --cluster_index $i \
         # --cells_per_sample $cells_per_sample \
         # --cores $cores \
         # --outdir $outdir
# done
# #Usage: bash jsub_submit.sh

#Here is script run_monocle_glm.R
suppressPackageStartupMessages({
	library(data.table)
	library(ArchR)
	library(matrixStats)
	library(reshape2)
	library(Matrix)
	library(optparse)
	library(monocle)
	library(dplyr)
})

option_list = list(
  make_option(c("--cluster_index"), type="integer"),
  make_option(c("--cells_per_sample"), type="integer", default=120),
  make_option(c("--cores"), type="integer", default=64),
  make_option(c("--outdir"), type="character", default="results")
)
opt <- parse_args(OptionParser(option_list=option_list))
cluster_idx <- opt$cluster_index
cells_per_sample <- opt$cells_per_sample
cores <- opt$cores
outdir <- opt$outdir

#Reading data
#proj <- readRDS("proj.rds") #ArchR proj
datafr <- readRDS("datafr.rds") #CRE x Cell
meta <- readRDS("meta.rds")
#meta = as.data.frame(proj@cellColData
common_cells <- intersect(colnames(datafr), rownames(meta))
datafr <- datafr[, common_cells]
meta <- meta[common_cells, , drop = FALSE]

cluster_name <- unique(meta$main)[cluster_idx]
cells_in_cluster <- rownames(meta)[meta$main == cluster_name]
if (length(cells_in_cluster) < 10) {
	message("Cluster ", cluster_name, " has too few cells, skipped.")
	quit(save = "no", status = 0)
}
set.seed(1)

ref_cells <- unlist(lapply(unique(meta$Sample), function(smp){
	cells_in_smp <- rownames(meta)[meta$Sample == smp & meta$main != cluster_name]
	if(length(cells_in_smp) <= cells_per_sample) return(cells_in_smp)
	sample(cells_in_smp, cells_per_sample)
}))
sample_cells <- c(cells_in_cluster, ref_cells)

submat <- (datafr[, sample_cells] > 0) * 1
submeta <- meta[sample_cells, , drop=FALSE]

pda <- data.frame(CellID = sample_cells,
	ReadDepth = submeta$nFrags,
	CellCluster = factor(ifelse(submeta$main == cluster_name, cluster_name, "Reference")))
rownames(pda) <- pda$CellID
pda <- new("AnnotatedDataFrame", data = pda)

fda <- data.frame(Peak = rownames(submat))
rownames(fda) <- fda$Peak
fda <- new("AnnotatedDataFrame", data = fda)

cds <- newCellDataSet(submat,
                      featureData = fda,
                      phenoData = pda,
                      expressionFamily = binomialff(),
                      lowerDetectionLimit = 1)
pData(cds)$Size_Factor <- 1
cds@expressionFamily@vfamily <- "binomialff"
message("cds finished......")

diff_test_helperBeta <- function(x, 
                             fullModelFormulaStr, 
                             reducedModelFormulaStr, 
                             expressionFamily, 
                             relative_expr,
                             weights,
                             disp_func=NULL,
                             verbose=FALSE
                             ){ 
  
  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, sep="")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, sep="")
  
  x_orig <- x
  disp_guess <- 0
  
  if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
    if (relative_expr == TRUE)
    {
      x <- x / Size_Factor
    }
    f_expression <- round(x)
    if (is.null(disp_func) == FALSE){
      disp_guess <- calculate_NB_dispersion_hint(disp_func, round(x_orig))
      if (is.null(disp_guess) == FALSE && disp_guess > 0 && is.na(disp_guess) == FALSE  ) {
        # FIXME: In theory, we could lose some user-provided parameters here
        # e.g. if users supply zero=NULL or something. 
        if (expressionFamily@vfamily == "negbinomial")
          expressionFamily <- negbinomial(isize=1/disp_guess)
        else
          expressionFamily <- negbinomial.size(size=1/disp_guess)
      }
    }
  }else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal")){
    f_expression <- x
  }else if (expressionFamily@vfamily %in% c("binomialff")){
    f_expression <- x/Size_Factor
    #f_expression[f_expression > 1] <- 1
  }else{
    f_expression <- log10(x)
  }
  
  test_res <- tryCatch({
    if (expressionFamily@vfamily %in% c("binomialff")){
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)                         
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))                    
      }
    }else{
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)                         
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))                    
      }
    }

    #print(full_model_fit)
    #print(coef(reduced_model_fit))
    compareModelsBeta <- function(full_models, reduced_models){
      stopifnot(length(full_models) == length(reduced_models))
      test_res <- mapply(function(x,y) { 
        if (is.null(x) == FALSE && is.null(y) == FALSE) {
          lrt <- VGAM::lrtest(x,y) 
          pval=lrt@Body["Pr(>Chisq)"][2,]
          family = x@family@vfamily
          if (length(family) > 1)
            family = family[1]
          beta = x@coefficients[2]
          data.frame(status = "OK", family=family, pval=pval,beta=beta)
        } else { data.frame(status = "FAIL", family=NA, pval=1.0,beta=0) } 
      } , full_models, reduced_models, SIMPLIFY=FALSE, USE.NAMES=TRUE)
      
      test_res <- do.call(rbind.data.frame, test_res)
      test_res$qval <- p.adjust(test_res$pval, method="BH")
      test_res
    }
    
    compareModelsBeta(list(full_model_fit), list(reduced_model_fit))
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { 
    if(verbose)
      print (e);
      data.frame(status = "FAIL", family=expressionFamily@vfamily, pval=1.0, qval=1.0, beta=0)
    #data.frame(status = "FAIL", pval=1.0) 
  }
  )
  test_res
}

differentialGeneTest <- function(cds, 
                                 fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
                                 reducedModelFormulaStr="~1", 
                                 relative_expr=TRUE,
                                 cores=1, 
                                 verbose=FALSE
){
  status <- NA
  if (relative_expr && cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
    if (is.null(sizeFactors(cds)) || sum(is.na(sizeFactors(cds)))){
      stop("Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet.")
    }
  }
  
  if (cores > 1){
    diff_test_res<-mcesApply(cds, 1, diff_test_helperBeta, 
                             c("BiocGenerics", "VGAM", "Matrix"), 
                             cores=cores, 
                             fullModelFormulaStr=fullModelFormulaStr,
                             reducedModelFormulaStr=reducedModelFormulaStr,
                             expressionFamily=cds@expressionFamily,
                             relative_expr=relative_expr,
                             disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                             verbose=verbose
                             #       ,
                             # backup_method = backup_method, 
                             # use_epislon = use_epislon, 
                             # stepsize = stepsize
    )
    diff_test_res
  }else{
    diff_test_res<-smartEsApply(cds,1,diff_test_helperBeta, 
                                convert_to_dense=TRUE,
                                fullModelFormulaStr=fullModelFormulaStr,
                                reducedModelFormulaStr=reducedModelFormulaStr, 
                                expressionFamily=cds@expressionFamily, 
                                relative_expr=relative_expr,
                                disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                                verbose=verbose
                                #          ,
                                # backup_method = backup_method, 
                                # use_epislon = use_epislon,
                                # stepsize = stepsize
                                
    )
    diff_test_res
  }
  
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)
  
  diff_test_res$qval <- 1
  diff_test_res$qval[which(diff_test_res$status == 'OK')] <- p.adjust(subset(diff_test_res, status == 'OK')[, 'pval'], method="BH")
  
  diff_test_res <- merge(diff_test_res, fData(cds), by="row.names")
  row.names(diff_test_res) <- diff_test_res[, 1] #remove the first column and set the row names to the first column
  diff_test_res[, 1] <- NULL 
  
  diff_test_res
}


diff_test <- differentialGeneTest(cds,
                                  fullModelFormulaStr = "~CellCluster + ReadDepth",
                                  reducedModelFormulaStr = "~ReadDepth",
                                  cores = cores-1)

saveRDS(diff_test, file = file.path(outdir, paste0("diff_test_cluster_", cluster_idx, ".rds")))
message("Finished cluster: ", cluster_name)
#Script run_monocle_glm.R ends.