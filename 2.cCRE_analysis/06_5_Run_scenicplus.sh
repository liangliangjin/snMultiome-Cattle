conda activate scenicplus
cd ~/03.snATAC/06.scenicplus/scplus_pipeline/Snakemake
snakemake --cores 64


# # The following is the config.yaml file for the snakemake pipeline.
# input_data:
  # cisTopic_obj_fname: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/cistopic.pkl"
  # GEX_anndata_fname: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/adata.h5ad"
  # region_set_folder: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/region_sets"
  # ctx_db_fname: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/10x_cattle_1kb_bg_with_mask.regions_vs_motifs.rankings.feather"
  # dem_db_fname: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/10x_cattle_1kb_bg_with_mask.regions_vs_motifs.scores.feather"
  # path_to_motif_annotations: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/motifs-v10nr_cattle.tbl"

# output_data:
  # # output for prepare_GEX_ACC .h5mu
  # combined_GEX_ACC_mudata: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/ACC_GEX.h5mu"
  # # output for motif enrichment results .hdf5
  # dem_result_fname: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/dem_results.hdf5"
  # ctx_result_fname: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/ctx_results.hdf5"
  # # output html for motif enrichment results .html
  # output_fname_dem_html: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/dem_results.html"
  # output_fname_ctx_html: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/ctx_results.html"
  # # output for prepare_menr .h5ad
  # cistromes_direct: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/cistromes_direct.h5ad"
  # cistromes_extended: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/cistromes_extended.h5ad"
  # # output tf names .txt
  # tf_names: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/tf_names.txt"
  # # output for download_genome_annotations .tsv
  # genome_annotation: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/genome_annotation.tsv"
  # chromsizes: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/chromsizes.tsv"
  # # output for search_space .tsb
  # search_space: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/search_space.tsv"
  # # output tf_to_gene .tsv
  # tf_to_gene_adjacencies: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/tf_to_gene_adj.tsv"
  # # output region_to_gene .tsv
  # region_to_gene_adjacencies: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/region_to_gene_adj.tsv"
  # # output eGRN .tsv
  # eRegulons_direct: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/eRegulon_direct.tsv"
  # eRegulons_extended: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/eRegulons_extended.tsv"
  # # output AUCell .h5mu
  # AUCell_direct: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/AUCell_direct.h5mu"
  # AUCell_extended: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/AUCell_extended.h5mu"
  # # output scplus mudata .h5mu
  # scplus_mdata: "/storage/public/home/2021060195/03.snATAC/06.scenicplus/outs/scplusmdata.h5mu"

# params_general:
  # temp_dir: "/storage/public/home/2021060195/temp"
  # n_cpu: 40
  # seed: 666

# params_data_preparation:
  # # Params for prepare_GEX_ACC
  # bc_transform_func: "\"lambda x: f'{x}'\""
  # is_multiome: True
  # key_to_group_by: ""
  # nr_cells_per_metacells: 10
  # # Params for prepare_menr
  # direct_annotation: "Direct_annot"
  # extended_annotation: "Orthology_annot"
  # # Params for download_genome_annotations
  # species: "cattle"
  # biomart_host: "http://www.ensembl.org"
  # # Params for search_space
  # search_space_upstream: "1000 150000"
  # search_space_downstream: "1000 150000"
  # search_space_extend_tss: "10 10"

# params_motif_enrichment:
  # species: "cattle"
  # annotation_version: "v10nr_clust"
  # motif_similarity_fdr: 0.001
  # orthologous_identity_threshold: 0.0
  # annotations_to_use: "Direct_annot Orthology_annot"
  # fraction_overlap_w_dem_database: 0.4
  # dem_max_bg_regions: 500
  # dem_balance_number_of_promoters: True
  # dem_promoter_space: 1_000
  # dem_adj_pval_thr: 0.05
  # dem_log2fc_thr: 1.0
  # dem_mean_fg_thr: 0.0
  # dem_motif_hit_thr: 3.0
  # fraction_overlap_w_ctx_database: 0.4
  # ctx_auc_threshold: 0.005
  # ctx_nes_threshold: 3.0
  # ctx_rank_threshold: 0.05




# params_inference:
  # # Params for tf_to_gene
  # tf_to_gene_importance_method: "GBM"
  # # Params regions_to_gene
  # region_to_gene_importance_method: "GBM"
  # region_to_gene_correlation_method: "SR"
  # # Params for eGRN inference
  # order_regions_to_genes_by: "importance"
  # order_TFs_to_genes_by: "importance"
  # gsea_n_perm: 1000
  # quantile_thresholds_region_to_gene: "0.85 0.90 0.95"
  # top_n_regionTogenes_per_gene: "5 10 15"
  # top_n_regionTogenes_per_region: ""
  # min_regions_per_gene: 0
  # rho_threshold: 0.05
  # min_target_genes: 10

