library(here)
library(tidyverse)
library(babelgene)

here::i_am("code/06_convert_deg_mouse_to_human_symbols.R")


# ========== 0.0 - Analysis settings ==========
# Positive log2FoldChange means higher in the first group named in the DESeq2
# contrast:
#   - old ID vs old CD: higher in old ID.
#   - young CD vs old CD: higher in young CD.
#
# The nominal P cutoff matches the existing DEG overlap workflow in
# code/04_deg_overlap_and_youth_shift.R.

deg_pvalue_cutoff = 0.05
min_ortholog_support = 3

ortholog_output_dir = here("data/processed/webCSEA/orthologs")
dir.create(ortholog_output_dir, recursive = TRUE, showWarnings = FALSE)

contrast_inputs = tribble(
  ~contrast_label, ~contrast_slug, ~deseq2_file, ~positive_direction_label, ~negative_direction_label,
  "old ID vs old CD", "old ID vs old CD",
  here("data/processed/deseq2/deseq2 ID vs CD old only cohort adjusted.csv"),
  "Upregulated", "Downregulated",
  "young CD vs old CD", "young CD vs old CD",
  here("data/processed/deseq2/deseq2 young CD vs old CD.csv"),
  "Upregulated", "Downregulated"
)


# ========== 1.0 - Convert one DEG contrast to human symbols ==========
convert_contrast_to_human_symbols = function(contrast_label,
                                             contrast_slug,
                                             deseq2_file,
                                             positive_direction_label,
                                             negative_direction_label) {
  de_results = read.csv(deseq2_file) %>%
    as_tibble()
  
  deg_mouse = de_results %>%
    dplyr::filter(
      !is.na(pvalue),
      pvalue <= deg_pvalue_cutoff,
      !is.na(log2FoldChange),
      log2FoldChange != 0
    ) %>%
    mutate(
      contrast = contrast_label,
      direction = ifelse(
        log2FoldChange > 0,
        positive_direction_label,
        negative_direction_label
      )
    ) %>%
    arrange(direction, pvalue, desc(abs(log2FoldChange)))
  
  if (nrow(deg_mouse) == 0) {
    warning("No nominal DEGs found for contrast: ", contrast_label)
    return(tibble())
  }
  
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
    select(contrast, direction, gene.id, baseMean, log2FoldChange, pvalue, padj)
  
  write.csv(
    deg_human_mapping,
    file.path(
      ortholog_output_dir,
      paste0(contrast_slug, " nominal DEG mouse to human orthologs.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    unmapped_mouse_genes,
    file.path(
      ortholog_output_dir,
      paste0(contrast_slug, " nominal DEG mouse genes without human ortholog.csv")
    ),
    row.names = FALSE
  )
  
  human_positive_genes = deg_human_mapping %>%
    dplyr::filter(direction == positive_direction_label, !is.na(human_symbol)) %>%
    pull(human_symbol) %>%
    unique()
  
  human_negative_genes = deg_human_mapping %>%
    dplyr::filter(direction == negative_direction_label, !is.na(human_symbol)) %>%
    pull(human_symbol) %>%
    unique()
  
  writeLines(
    human_positive_genes,
    file.path(
      ortholog_output_dir,
      paste0(contrast_slug, " upregulated human symbols.txt")
    )
  )
  
  writeLines(
    human_negative_genes,
    file.path(
      ortholog_output_dir,
      paste0(contrast_slug, " downregulated human symbols.txt")
    )
  )
  
  deg_human_mapping
}


# ========== 2.0 - Run all configured contrasts ==========
all_deg_human_mapping = pmap_dfr(contrast_inputs, convert_contrast_to_human_symbols)

write.csv(
  all_deg_human_mapping,
  file.path(ortholog_output_dir, "all nominal DEG mouse to human orthologs.csv"),
  row.names = FALSE
)


# ========== 3.0 - Console summary ==========
summary_tbl = all_deg_human_mapping %>%
  group_by(contrast, direction) %>%
  summarise(
    mouse_deg_n = n_distinct(gene.id),
    mapped_human_symbol_n = n_distinct(human_symbol[!is.na(human_symbol)]),
    unmapped_mouse_gene_n = n_distinct(gene.id[is.na(human_symbol)]),
    .groups = "drop"
  )

print(summary_tbl)
