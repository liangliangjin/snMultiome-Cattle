library(biomaRt)
library(dplyr)
orth_tf <- read.table("Bos_taurus_ortholog_to_human.txt",sep="\t", header=TRUE)
colnames(orth_tf) <- c("Species", "bovine_ens", "cov_bov", "id_bov", "human_ens", "cov_hum", "id_hum", "Ortholog_Species")
unique_human_ensg <- unique(orth_tf$human_ens)
unique_bovine_ensbt <- unique(orth_tf$bovine_ens)
human_mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
bovine_mart <- useEnsembl(biomart = "genes", dataset = "btaurus_gene_ensembl")
human_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique_human_ensg,
  mart = human_mart
) %>%
  filter(hgnc_symbol != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  rename(human_ens = ensembl_gene_id, human_symbol = hgnc_symbol)
  
bovine_symbols <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = unique_bovine_ensbt,
  mart = bovine_mart
) %>%
  filter(external_gene_name != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  rename(bovine_ens = ensembl_gene_id, bovine_symbol = external_gene_name)
  
orth_tf_with_sym <- orth_tf %>%
  left_join(human_symbols, by = "human_ens") %>%
  left_join(bovine_symbols, by = "bovine_ens") %>%
  mutate(
    bovine_tf_final = case_when(
      !is.na(bovine_symbol) & bovine_symbol != "" ~ bovine_symbol,
      !is.na(bovine_ens) ~ bovine_ens,
      TRUE ~ NA_character_
    )
  )


original_tbl <- read.table("/home/Jingliangliang/SC/scenicplus/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl", sep="\t", header = T, comment.char = "", check.names = FALSE)
mapping <- setNames(orth_tf_with_sym$bovine_tf_final, orth_tf_with_sym$human_symbol)
cattle_tbl <- original_tbl %>%
  mutate(
    bovine_gene = mapping[gene_name],
    matched     = !is.na(bovine_gene)
  ) %>%
  filter(matched) %>%  
  mutate(gene_name = bovine_gene) %>%
  select(-bovine_gene, -matched) 
write.table(cattle_tbl,file = "motifs-v10nr_cattle.tbl", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

