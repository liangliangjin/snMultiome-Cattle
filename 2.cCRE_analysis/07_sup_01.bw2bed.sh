# Step 1: BigWig to bed file
## Download hg38.phyloP100way.bw
bigWigToBedGraph hg38.phyloP100way.bw hg38.phyloP100way.bed

# Step 2: Chromosomal resolution
awk '{print > $1".bed"}' hg38.phyloP100way.bed

# Step 3: Bedtools intersect
for chr in {1..22} X Y; do
  bedtools intersect -a chr${chr}.bed -b region_file.bed -wa > chr${chr}_filter.bed
done

# Step 4: Add line numbers and sort
awk '{print $0 "\t" NR}' region_file.bed > region_with_id.bed
bedtools sort -i region_with_id.bed > sorted_with_id.bed

# Step 5: Chromosomal resolution
awk '{print > $1"_region.bed"}' sorted_with_id.bed