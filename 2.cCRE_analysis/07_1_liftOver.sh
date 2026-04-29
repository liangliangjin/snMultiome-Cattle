RAW_INPUT="peakset_478422.bed"
#0-based bed to 1-based bed
INPUT="peakset_478422.with_id.bed"
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$1":"$2+1"-"$3,1}' "$RAW_INPUT" > "$INPUT"

declare -A forward_chains=(
    ["Hg38"]="chain/bosTau9ToHg38.over.chain.gz"
    ["Mm10"]="chain/bosTau9ToMm10.over.chain.gz"
    ["SusScr11"]="chain/bosTau9ToSusScr11.over.chain.gz"
    ["EquCab3"]="chain/bosTau9ToEquCab3.over.chain.gz"
    ["Oar_rambouillet"]="chain/bosTau9ToGCF_016772045.2.over.chain.gz"
)

declare -A reverse_chains=(
    ["Hg38"]="chain/hg38ToBosTau9.over.chain.gz"
    ["Mm10"]="chain/mm10ToBosTau9.over.chain.gz"
    ["SusScr11"]="chain/susScr11ToBosTau9.over.chain.gz"
    ["EquCab3"]="chain/equCab3ToBosTau9.over.chain.gz"
    ["Oar_rambouillet"]="chain/GCF_016772045.2ToBosTau9.over.chain.gz"
)
for species in Hg38 Mm10 EquCab3 SusScr11 Oar_rambouillet; do
    echo "Processing $species ..."
    forward_chain="${forward_chains[$species]}"
    reverse_chain="${reverse_chains[$species]}"
    forward_out="bos_${species}_0.5.bed"
    forward_unmap="bos_unmap_0.5_${species}.bed"
    reciprocal_out="reciprocal_${species}.bed"
    reciprocal_unmap="reciprocal_unmap_${species}.bed"
    # bosTau9 to target species
    $LIFTOVER -minMatch=0.5 -multiple $INPUT $forward_chain $forward_out $forward_unmap
    # target species back to bosTau9
    $LIFTOVER -minMatch=0.5 -multiple $forward_out $reverse_chain $reciprocal_out $reciprocal_unmap
    # Keep only those records that return to the original peak
	bos_conserved="bos_conserved_${species}_0.5.bed"
    conserved_ids="conserved_${species}_0.5.ids"
    target_conserved="${species}_conserved_0.5.bed"
    awk 'BEGIN{FS=OFS="\t"}
    {
        split($4,a,/[:-]/)
        if ($1 != a[1]) next
        s = ($2 > a[2] ? $2 : a[2])
        e = ($3 < a[3] ? $3 : a[3])
        if (e <= s) next
        overlap = e - s
        len_ret = $3 - $2
        len_org = a[3] - a[2]
        if (overlap/len_ret >= 0.5 && overlap/len_org >= 0.5)
            print $4
    }' "$reciprocal_out" | sort -u > "$conserved_ids"
    awk 'NR==FNR{keep[$1]=1; next} ($4 in keep)' "$conserved_ids" "$INPUT" > "$bos_conserved"
    awk 'NR==FNR{keep[$1]=1; next} ($4 in keep)' "$conserved_ids" "$forward_out" > "$target_conserved"
done