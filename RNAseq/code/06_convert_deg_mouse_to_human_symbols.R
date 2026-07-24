library(here)
library(tidyverse)
library(babelgene)

here::i_am("code/06_convert_deg_mouse_to_human_symbols.R")


# ========== 0.0 - Analysis settings ==========
# Positive log2FoldChange means higher in old ID than old CD.
# The nominal P cutoff matches the existing DEG overlap workflow in
# code/04_deg_overlap_and_youth_shift.R.

deg_pvalue_cutoff = 0.05
min_ortholog_support = 3

dir.create(here("data/processed/orthologs"), recursive = TRUE, showWarnings = FALSE)


# ========== 1.0 - Read old ID vs old CD DESeq2 results ==========
de_results = read.csv(
  here("data/processed/deseq2/deseq2 ID vs CD old only cohort adjusted.csv")
) %>%
  as_tibble()


# ========== 2.0 - Define significant up/down mouse DEG sets ==========
deg_mouse = de_results %>%
  dplyr::filter(
    !is.na(pvalue),
    pvalue <= deg_pvalue_cutoff,
    !is.na(log2FoldChange),
    log2FoldChange != 0
  ) %>%
  mutate(
    direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")
  ) %>%
  arrange(direction, pvalue, desc(abs(log2FoldChange)))


# ========== 3.0 - Convert mouse symbols to human ortholog symbols ==========
mouse_to_human = babelgene::orthologs(
  genes = unique(deg_mouse$gene.id),
  species = "mouse",
  human = FALSE,
  min_support = min_ortholog_support,
  top = TRUE
) %>%
  as_tibble() %>%
  transmute(
    gene.id = symbol,
    human_symbol,
    mouse_entrez = entrez,
    human_entrez,
    mouse_ensembl = ensembl,
    human_ensembl,
    ortholog_support_n = support_n,
    ortholog_support_sources = support
  )

deg_human_mapping = deg_mouse %>%
  left_join(mouse_to_human, by = "gene.id") %>%
  arrange(direction, pvalue, desc(abs(log2FoldChange)), gene.id)

unmapped_mouse_genes = deg_human_mapping %>%
  dplyr::filter(is.na(human_symbol)) %>%
  select(direction, gene.id, baseMean, log2FoldChange, pvalue, padj)


# ========== 4.0 - Write mapping tables ==========
write.csv(
  deg_human_mapping,
  here("data/processed/orthologs/old ID vs old CD nominal DEG mouse to human orthologs.csv"),
  row.names = FALSE
)

write.csv(
  unmapped_mouse_genes,
  here("data/processed/orthologs/old ID vs old CD nominal DEG mouse genes without human ortholog.csv"),
  row.names = FALSE
)


# ========== 5.0 - Write copy-pasteable human symbol lists ==========
human_up_genes = deg_human_mapping %>%
  dplyr::filter(direction == "Upregulated", !is.na(human_symbol)) %>%
  pull(human_symbol) %>%
  unique()

human_down_genes = deg_human_mapping %>%
  dplyr::filter(direction == "Downregulated", !is.na(human_symbol)) %>%
  pull(human_symbol) %>%
  unique()

writeLines(
  human_up_genes,
  here("data/processed/orthologs/old ID vs old CD upregulated human symbols.txt")
)

writeLines(
  human_down_genes,
  here("data/processed/orthologs/old ID vs old CD downregulated human symbols.txt")
)


# ========== 6.0 - Console summary ==========
summary_tbl = deg_human_mapping %>%
  group_by(direction) %>%
  summarise(
    mouse_deg_n = n_distinct(gene.id),
    mapped_human_symbol_n = n_distinct(human_symbol[!is.na(human_symbol)]),
    unmapped_mouse_gene_n = n_distinct(gene.id[is.na(human_symbol)]),
    .groups = "drop"
  )

print(summary_tbl)
