conda activate create_cistarget_databases

create_cistarget_databases_dir=/storage/public/home/2021060195/03.snATAC/06.scenicplus/create_cisTarget_databases
reference_dir=/storage/public/home/2021060195/03.snATAC/06.scenicplus/create_cisTarget_databases

${create_cistarget_databases_dir}/create_fasta_with_padded_bg_from_bed.sh \
    ${reference_dir}/ARS-UCD1.2_Btau5.0.1Y.with_chr.fa \
    ${reference_dir}/ARS-UCD1.2.chr30.with_chr.sizes \
    consensus.peakset.bed \
    cattle_consensus.with_1kb_bg.fa \
    1000 \
    yes

OUT_DIR=""${PWD}""
FASTA_FILE="${OUT_DIR}/cattle_consensus.with_1kb_bg.fa"
CBDIR="${OUT_DIR}/create_cisTarget_databases/aertslab_motif_colleciton/v10nr_clust_public/singletons"
MOTIF_LIST="${OUT_DIR}/create_cisTarget_databases/motifs.txt"
DATABASE_PREFIX="10x_cattle_1kb_bg_with_mask"

python ${create_cistarget_databases_dir}/create_cistarget_motif_databases.py \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 60
