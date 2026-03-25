## Generate the fragments files with the sample prefix barcode.
#
# for f in ~/SC/*/outs/atac_fragments.tsv.gz; do
  # tissue=$(basename $(dirname $(dirname $f)))
  # output="${f%.tsv.gz}_prefixed.tsv.gz"
  # zcat "$f" | awk -v t="$tissue" '{if(NR>2) $4=t"#"$4; print}' OFS="\t" | bgzip > "$output"
  # tabix -p bed "$output"
  # echo "finished: $tissue"
# done


SeuratObject <- readRDS("SeuratObject_wnn.rds")
DefaultAssay(SeuratObject) <- "peaks"

samples <- unique(SeuratObject$Sample)
frag_paths <- paste0("~/SC/", samples, "/outs/atac_fragments_prefixed.tsv.gz")
names(frag_paths) <- samples
valid_paths <- frag_paths[file.exists(frag_paths)]
fragments_list <- list()
for (samp in names(valid_paths)) {
  path <- valid_paths[samp]
  cells_samp <- colnames(SeuratObject)[SeuratObject$Sample == samp]
  frag <- CreateFragmentObject(
    path = path,
    cells = cells_samp
  )
  fragments_list[[samp]] <- frag
}
SeuratObject@assays$peaks@fragments <- fragments_list
print(Fragments(SeuratObject[["peaks"]]))

dir.create("~/SC/chrombpnet/celltype_fragments/", recursive = TRUE, showWarnings = FALSE)
# prepare fragment file for bias model training
# Train a bias model only in a tissue. This bias model should generalize across tissues (assuming they follow a similar sequencing protocol).

SplitFragments(
  SeuratObject,
  assay = "peaks",
  group.by = "Sample",
  idents = "spleen0" ,
  outdir = "chrombpnet/celltype_fragments/",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE, 
  file.suffix = "_input_bias"
)

# split fragment by celltype
Idents(SeuratObject) <- "main"
SplitFragments(
  SeuratObject,
  assay = "peaks",
  group.by = "main",
  outdir = "chrombpnet/celltype_fragments/",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE
)

# split fragment by celltype and cattle
table(SeuratObject$celltype_cattle)
Idents(SeuratObject) <- "celltype_cattle"
Idents(SeuratObject) %>% unique()
dir.create("~/SC/chrombpnet/celltype_fragments_cattle/", recursive = TRUE, showWarnings = FALSE)
SplitFragments(
  SeuratObject,
  assay = "peaks",
  group.by = "celltype_cattle",
  outdir = "chrombpnet/celltype_fragments_cattle/",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE
)


# prepare bed files for chrombpnet interpretation (e.g., cell type- and injury-specific regions)
celltype_injury_dars <- readRDS("markers/markers_Uvsothers_peaks_all.rds")
for(i in names(celltype_injury_dars)){
  celltype_injury_dars[[i]]$peak <- sapply(strsplit(rownames(celltype_injury_dars[[i]]), "[.]"), "[[", 2)
  celltype_injury_dars[[i]] <- celltype_injury_dars[[i]] %>%
    dplyr::filter(p_val < 0.05 & avg_log2FC > 0.5)
}

for(i in names(celltype_injury_dars)){
  chr <- sapply(strsplit(celltype_injury_dars[[i]]$peak, "-"), "[", 1) 
  start <- sapply(strsplit(celltype_injury_dars[[i]]$peak, "-"), "[", 2) 
  end <- sapply(strsplit(celltype_injury_dars[[i]]$peak, "-"), "[", 3) 
  bed <- cbind(chr, start, end)
  bed_file <- paste0(i, "_I_dars.bed")
  write.table(bed, bed_file,row.names = F,col.names = F, sep="\t", quote=FALSE)
}