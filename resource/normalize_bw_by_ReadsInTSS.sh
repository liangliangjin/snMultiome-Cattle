module add ucscTools

scale_table=celltype_bw_norm_factors.tsv

##scale_table calculation (Rscript):
# library(dplyr)
# norm_df <- as.data.frame(SeuratObject@meta.data) %>%
  # group_by(celltype_cattle) %>%
  # summarise(
    # nCells = n(),
    # sum_nFrags = sum(nFrags, na.rm = TRUE),
    # sum_ReadsInTSS = sum(ReadsInTSS, na.rm = TRUE),
    # median_TSSEnrichment = median(TSSEnrichment, na.rm = TRUE)
  # ) %>%
  # mutate(
    # ReadsInTSS_scale = 10000 / sum_ReadsInTSS,
    # CPM_scale = 1000000 / sum_nFrags
  # )
  
# colnames(norm_df)
# [1] "celltype_cattle"      "nCells"               "sum_nFrags"           "sum_ReadsInTSS"       "median_TSSEnrichment" "ReadsInTSS_scale"     "CPM_scale"  

chrom_size=ARSUCD1.2.chr.sizes
bw_dir=.
out_dir=bw_ArchR_ReadsInTSS_norm
mkdir -p ${out_dir}

# Adjustable: Minimum cell count threshold
min_cells=20

tail -n +2 ${scale_table} | while IFS=$'\t' read celltype nCells sum_nFrags sum_ReadsInTSS median_TSSEnrichment ReadsInTSS_scale CPM_scale
do
    in_bw="${bw_dir}/${celltype}.bw"
    out_bg="${out_dir}/${celltype}.ArchRnorm.bedGraph"
    out_bw="${out_dir}/${celltype}.ArchRnorm.bw"
    if [ "${nCells}" -lt "${min_cells}" ]; then
        echo "Skip ${celltype}: nCells=${nCells} < ${min_cells}"
        continue
    fi
    if [ ! -f "${in_bw}" ]; then
        echo "Warning: ${in_bw} not found, skip"
        continue
    fi
    echo "Processing ${celltype}"
    echo "  nCells=${nCells}, ReadsInTSS=${sum_ReadsInTSS}, scale=${ReadsInTSS_scale}"
    bigWigToBedGraph "${in_bw}" stdout | \
    awk -v s="${ReadsInTSS_scale}" 'BEGIN{OFS="\t"}{$4=$4*s; print}' > "${out_bg}"
    bedGraphToBigWig "${out_bg}" "${chrom_size}" "${out_bw}"
    rm "${out_bg}"
done
