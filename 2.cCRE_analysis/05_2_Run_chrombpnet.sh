#conda activate chrombpnet
source ~/conda_envs/chrombpnet/bin/activate
sort -k 1,1V -k 2,2n all_input_bias.bed > sorted_fragments_bias.tsv
bgzip -@ 20 sorted_fragments_bias.tsv
tabix -f -p bed sorted_fragments_bias.tsv.gz

macs2 callpeak -f AUTO -t sorted_fragments_bias.tsv.gz -g 2.7e+09 -B -p 0.01 --nomodel --extsize 200 --outdir macs2_from_frag -n bias_fragments

# remove blacklist
bedtools slop -i ARSUCD1.2_blacklist.bed -g ARSUCD1.2.chr30.size -b 1057 > temp.bed
bedtools intersect -v -a ./macs2_from_frag/bias_fragments.narrowPeak -b temp.bed  > bias_peaks_no_blacklist.bed


sed 's/^>/>chr/' /storage/public/home/2021060195/reference/01.ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa > /storage/public/home/2021060195/03.snATAC/07.chrombpnet/cattle_genome_chr.fa

chrombpnet prep splits -c ARSUCD1.2.chr30.size -tcr chr1 chr3 -vcr chr10 chr20 -op splits/fold_0

chrombpnet prep nonpeaks -g cattle_genome_chr.fa -p bias_peaks_no_blacklist.bed -c ARSUCD1.2.chr30.size -fl splits/fold_0.json -br ARSUCD1.2_blacklist.bed -o output_nonpeaks

# train bias model
chrombpnet bias pipeline \
       -ifrag sorted_fragments_bias.tsv.gz \
       -d "ATAC" \
       -g cattle_genome_chr.fa \
       -c ARSUCD1.2.chr30.size \
       -p bias_peaks_no_blacklist.bed \
       -n output_nonpeaks_negatives.bed \
       -fl splits/fold_0.json \
       -b 0.5 \
       -o bias/ \
       -fp bias_model


# train cell type models (same for cell type-, injury-specific models) 
#creat celltype.list
while read i; do
 grep '^chr' chrombpnet/celltype_fragments/${i}.bed > celltype_fragments_multiome_only/${i}_chr.bed
 macs2 callpeak -f AUTO \
 -t celltype_fragments_multiome_only/${i}_chr.bed \
 -g 2.7e+09 -B -p 0.01 --nomodel --extsize 200 \
 --outdir macs2_celltype_peaks -n ${i}_fragments
done < celltype.list

mkdir -p filtered_peaks
mkdir -p negatives
bedtools slop -i ARSUCD1.2_blacklist.bed -g ARSUCD1.2.chr30.size -b 1057 > temp.bed

# create peak sets without blacklist regions and prepare non-peak regions
while read i; do
  bedtools intersect -v -a macs2_celltype_peaks/${i}_fragments_peaks.narrowPeak -b temp.bed  > filtered_peaks/${i}_peaks_no_blacklist.bed
  grep '^chr' filtered_peaks/${i}_peaks_no_blacklist.bed > filtered_peaks/${i}_peaks_no_blacklist_chr.bed
  chrombpnet prep nonpeaks \
    -g cattle_genome_chr.fa \
    -p filtered_peaks/${i}_peaks_no_blacklist_chr.bed \
    -c ARSUCD1.2.chr30.size \
    -fl splits/fold_0.json \
	-br ARSUCD1.2_blacklist.bed \
    -o negatives/${i}_nonpeaks
done < celltype.list

while read i; do
  chrombpnet pipeline \
    -ifrag celltype_fragments_multiome_only/${i}_chr.bed \
    -d "ATAC" \
    -g cattle_genome_chr.fa \
    -c ARSUCD1.2.chr30.size \
    -p filtered_peaks/${i}_peaks_no_blacklist_chr.bed \
    -n negatives/${i}_nonpeaks_negatives.bed \
    -fl splits/fold_0.json \
    -b bias_fold_0/models/bias_model_bias.h5 \
    -o models/${i} \
    > logs/${i}_logs.txt 2>&1
done < celltype.list


# predict other sets of peaks (e.g., DARs enriched in specific cell types or in catlle subspecies)
# obtain count and profile predictions
# compute base pair-level contribution scores
# identify important motifs and link to TFs
dars_intersect_peaks_path = "./dars_edgeR_filter/"
while read i; do
  chrombpnet pred_bw -bm models/${i}_fl0/models/bias_model_scaled.h5 \
    -cm models_fold_0/${i}/models/chrombpnet.h5 \
    -r ${dars_intersect_peaks_path}/${i}_fragments_peaks.narrowPeak \
    -c ARSUCD1.2.chr30.size -g cattle_genome_chr.fa \
    -op interpret_out/${i}/${i}_dars_bigwig
    
  chrombpnet contribs_bw -m models_fold_0/${i}/models/chrombpnet.h5 \
    -r ${dars_intersect_peaks_path}/${i}_fragments_peaks.narrowPeak \
    -c ARSUCD1.2.chr30.size -g cattle_genome_chr.fa \
    -op interpret_out/${i}/${i}_dars_contrib
 
  modisco motifs -i interpret_out/${i}/${i}_dars_contrib.counts_scores.h5 -n 1000000 -o interpret_out/${i}/${i}_modisco.h5 -v
  modisco report -i interpret_out/${i}/${i}_modisco.h5 -o interpret_out/${i}/modisco_report/ -s interpret_out/${i}/modisco_report/ -m motifs.meme.txt
done < celltype.list