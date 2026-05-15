library(here)
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(cowplot)

here::i_am("code/02 Young old CD DESeq2 GSEA overlap.R")


# ========== 0.0 - Analysis settings ==========
# This script compares two pathway-level contrasts:
#   1. young CD vs old CD: aging signal among control-diet mice.
#   2. old ID vs old CD: inulin diet signal among old mice.
#
# Positive log2FoldChange/stat/NES values mean higher in the first group listed:
#   - young CD vs old CD: higher in young CD.
#   - old ID vs old CD: higher in old ID.
#
# The old ID vs old CD contrast is re-run here so the overlap figure is built
# from one self-contained script rather than depending on previously saved files.

gsea_rank_column = "stat"
gsea_padj_cutoff = 0.10
top_overlap_pathways_per_panel = 20
gsea_seed = 20260514


# ========== 1.0 - Read counts and metadata ==========
meta = read.csv(here("data/metadata/metadata.csv")) %>%
  mutate(
    sample = paste0("r", sample),
    cohort = factor(cohort),
    age.grp = factor(age.grp, levels = c("young", "old")),
    diet = factor(diet, levels = c("CD", "ID"))
  )

cnts = read.csv(here("data/processed/counts matrix with mouse gene symbols.csv")) %>%
  column_to_rownames(var = "gene_id")

cnts = cnts[, meta$sample]
stopifnot(identical(colnames(cnts), meta$sample))


# ========== 2.0 - DESeq2 helper for the two contrasts ==========
run_deseq2_contrast = function(samples, design_formula, contrast_vector, comparison_label, output_file) {
  sample_ids = samples$sample
  
  x = cnts[, sample_ids]
  x[is.na(x)] = 0
  x = round(as.matrix(x))
  storage.mode(x) = "integer"
  
  # Keep genes with visible expression in at least three samples. This is a
  # light prefilter; DESeq2 still performs its own independent filtering.
  keep_genes = rowSums(x >= 10) >= 3
  x = x[keep_genes, ]
  
  col_data = samples %>%
    column_to_rownames(var = "sample")
  
  stopifnot(identical(colnames(x), rownames(col_data)))
  
  dds = DESeqDataSetFromMatrix(
    countData = x,
    colData = col_data,
    design = design_formula
  )
  
  dds = DESeq(dds)
  res = results(dds, contrast = contrast_vector)
  
  res_tbl = data.frame(res) %>%
    rownames_to_column(var = "gene.id") %>%
    as_tibble() %>%
    arrange(padj, pvalue) %>%
    mutate(
      design = as.character(design_formula)[2],
      comparison = comparison_label
    )
  
  write.csv(res_tbl, output_file, row.names = FALSE)
  
  list(
    dds = dds,
    res = res,
    res_tbl = res_tbl,
    output_file = output_file
  )
}


# ========== 3.0 - DESeq2: young CD vs old CD ==========
# Cohort is fully confounded with age in the control-diet subset:
# old CD samples are cohorts 1/2 and young CD samples are cohorts 4/5.
# Therefore the model is age only.
samples_age_cd = meta %>%
  dplyr::filter(diet == "CD") %>%
  arrange(age.grp, cohort, sample)

table(samples_age_cd$age.grp, samples_age_cd$cohort)

age_cd_deseq = run_deseq2_contrast(
  samples = samples_age_cd,
  design_formula = ~ age.grp,
  contrast_vector = c("age.grp", "young", "old"),
  comparison_label = "young CD vs old CD",
  output_file = here("data/processed/deseq2 young CD vs old CD.csv")
)

age_cd_res = age_cd_deseq$res_tbl


# ========== 4.0 - DESeq2: old ID vs old CD ==========
# This repeats the diet contrast from the earlier script, adjusted for cohort.
samples_diet_old = meta %>%
  dplyr::filter(age.grp == "old") %>%
  arrange(cohort, diet, sample)

table(samples_diet_old$diet, samples_diet_old$cohort)

diet_old_deseq = run_deseq2_contrast(
  samples = samples_diet_old,
  design_formula = ~ cohort + diet,
  contrast_vector = c("diet", "ID", "CD"),
  comparison_label = "old ID vs old CD",
  output_file = here("data/processed/deseq2 old ID vs old CD for age overlap.csv")
)

diet_old_res = diet_old_deseq$res_tbl


# ========== 5.0 - GO BP and GO CC GSEA ==========
gene_symbol_to_entrez = AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(c(age_cd_res$gene.id, diet_old_res$gene.id)),
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
) %>%
  as_tibble() %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, ENTREZID)

run_go_gsea = function(de_tbl, ontology, comparison_slug, comparison_label) {
  ranked_genes = de_tbl %>%
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(pvalue),
      !is.na(.data[[gsea_rank_column]])
    ) %>%
    mutate(rank_value = .data[[gsea_rank_column]]) %>%
    inner_join(gene_symbol_to_entrez, by = c("gene.id" = "SYMBOL")) %>%
    arrange(desc(abs(rank_value))) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  gene_ranks = ranked_genes$rank_value
  names(gene_ranks) = ranked_genes$ENTREZID
  gene_ranks = sort(gene_ranks, decreasing = TRUE)
  
  set.seed(gsea_seed)
  gsea = gseGO(
    geneList = gene_ranks,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    seed = TRUE,
    verbose = FALSE
  )
  
  gsea_tbl = as.data.frame(gsea) %>%
    as_tibble() %>%
    arrange(p.adjust) %>%
    mutate(
      ontology = ontology,
      comparison = comparison_label,
      direction = ifelse(NES > 0, "Up", "Down")
    )
  
  ontology_label = ifelse(ontology == "BP", "bp", "cc")
  write.csv(
    gsea_tbl,
    here(
      "data/processed",
      paste0("go ", ontology_label, " gsea ", comparison_slug, " stat ranking.csv")
    ),
    row.names = FALSE
  )
  
  list(
    gsea = gsea,
    gsea_tbl = gsea_tbl
  )
}

age_gobp = run_go_gsea(age_cd_res, "BP", "young CD vs old CD", "young CD vs old CD")
age_gocc = run_go_gsea(age_cd_res, "CC", "young CD vs old CD", "young CD vs old CD")
diet_gobp = run_go_gsea(diet_old_res, "BP", "old ID vs old CD for age overlap", "old ID vs old CD")
diet_gocc = run_go_gsea(diet_old_res, "CC", "old ID vs old CD for age overlap", "old ID vs old CD")

age_gsea_tbl = bind_rows(age_gobp$gsea_tbl, age_gocc$gsea_tbl)
diet_gsea_tbl = bind_rows(diet_gobp$gsea_tbl, diet_gocc$gsea_tbl)


# ========== 6.0 - Plot GO BP/CC GSEA for young CD vs old CD ==========
age_gsea_plot_data = age_gsea_tbl %>%
  dplyr::filter(p.adjust <= gsea_padj_cutoff) %>%
  mutate(
    ontology = recode(ontology, BP = "GO Biological Process", CC = "GO Cellular Component"),
    direction_label = ifelse(NES > 0, "Higher in young CD", "Higher in old CD"),
    pathway = str_wrap(Description, width = 58)
  ) %>%
  group_by(ontology, direction_label) %>%
  slice_min(p.adjust, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(ontology, direction_label, NES) %>%
  mutate(
    pathway_for_plot = paste(pathway, ontology, direction_label, sep = "___"),
    pathway_for_plot = factor(pathway_for_plot, levels = unique(pathway_for_plot)),
    neg_log10_padj = -log10(pmax(p.adjust, .Machine$double.xmin))
  )

if (nrow(age_gsea_plot_data) > 0) {
  age_nes_limit = max(abs(age_gsea_plot_data$NES), na.rm = TRUE)
  age_nes_limit = ceiling(age_nes_limit * 2) / 2
  
  age_gsea_plot = ggplot(age_gsea_plot_data, aes(NES, pathway_for_plot)) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_point(aes(size = setSize, color = neg_log10_padj), alpha = 0.9) +
    facet_grid(ontology + direction_label ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_color_viridis_c(name = expression(-log[10]~adjusted~italic(P))) +
    coord_cartesian(xlim = c(-age_nes_limit, age_nes_limit)) +
    labs(
      title = "GO GSEA: young CD vs old CD",
      x = "Normalized enrichment score",
      y = NULL,
      size = "Gene set size"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8, lineheight = 0.9),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 28),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
} else {
  age_gsea_plot = ggdraw() +
    draw_label("No young CD vs old CD GO GSEA terms passed the adjusted P cutoff", size = 12)
}

ggsave(
  here("plots/go bp cc gsea young CD vs old CD stat ranking.pdf"),
  age_gsea_plot,
  width = 11,
  height = 12
)


# ========== 7.0 - Compare overlapping up and down pathways ==========
age_sig_pathways = age_gsea_tbl %>%
  dplyr::filter(p.adjust <= gsea_padj_cutoff) %>%
  dplyr::select(
    ontology,
    ID,
    Description,
    direction,
    age_NES = NES,
    age_p.adjust = p.adjust,
    age_setSize = setSize
  )

diet_sig_pathways = diet_gsea_tbl %>%
  dplyr::filter(p.adjust <= gsea_padj_cutoff) %>%
  dplyr::select(
    ontology,
    ID,
    Description,
    direction,
    diet_NES = NES,
    diet_p.adjust = p.adjust,
    diet_setSize = setSize
  )

overlap_pathways = age_sig_pathways %>%
  inner_join(
    diet_sig_pathways,
    by = c("ontology", "ID", "Description", "direction")
  ) %>%
  mutate(
    direction_label = case_when(
      direction == "Up" ~ "Up in young CD and old ID",
      direction == "Down" ~ "Down in young CD and old ID",
      TRUE ~ direction
    ),
    best_p.adjust = pmin(age_p.adjust, diet_p.adjust, na.rm = TRUE),
    worst_p.adjust = pmax(age_p.adjust, diet_p.adjust, na.rm = TRUE)
  ) %>%
  arrange(ontology, direction, worst_p.adjust)

write.csv(
  overlap_pathways,
  here("data/processed/go gsea overlapping pathways young CD vs old CD and old ID vs old CD.csv"),
  row.names = FALSE
)

comparison_outline_colors = c(
  "young CD vs old CD" = "#2C7FB8",
  "old ID vs old CD" = "#238B45"
)

comparison_shapes = c(
  "young CD vs old CD" = 21,
  "old ID vs old CD" = 24
)

make_overlap_plot = function(ontology_code, ontology_label, output_file) {
  overlap_plot_terms = overlap_pathways %>%
    dplyr::filter(ontology == ontology_code) %>%
    group_by(direction_label) %>%
    slice_min(worst_p.adjust, n = top_overlap_pathways_per_panel, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      pathway = str_wrap(Description, width = 54),
      pathway_for_plot = paste(pathway, direction_label, sep = "___")
    )
  
  overlap_plot_data = bind_rows(
    overlap_plot_terms %>%
      transmute(
        direction_label,
        pathway_for_plot,
        pathway,
        comparison = "young CD vs old CD",
        NES = age_NES,
        p.adjust = age_p.adjust,
        setSize = age_setSize
      ),
    overlap_plot_terms %>%
      transmute(
        direction_label,
        pathway_for_plot,
        pathway,
        comparison = "old ID vs old CD",
        NES = diet_NES,
        p.adjust = diet_p.adjust,
        setSize = diet_setSize
      )
  ) %>%
    mutate(
      comparison = factor(comparison, levels = names(comparison_outline_colors)),
      neg_log10_padj = -log10(pmax(p.adjust, .Machine$double.xmin))
    ) %>%
    arrange(direction_label, NES) %>%
    mutate(
      pathway_for_plot = factor(pathway_for_plot, levels = unique(pathway_for_plot))
    )
  
  if (nrow(overlap_plot_data) > 0) {
    overlap_nes_limit = max(abs(overlap_plot_data$NES), na.rm = TRUE)
    overlap_nes_limit = ceiling(overlap_nes_limit * 2) / 2
    
    overlap_plot = ggplot(overlap_plot_data, aes(NES, pathway_for_plot)) +
      geom_vline(xintercept = 0, color = "grey78", linewidth = 0.35) +
      geom_segment(
        data = overlap_plot_terms,
        aes(
          x = age_NES,
          xend = diet_NES,
          y = pathway_for_plot,
          yend = pathway_for_plot
        ),
        inherit.aes = FALSE,
        color = "grey72",
        linewidth = 0.45
      ) +
      geom_point(
        aes(
          shape = comparison,
          color = comparison,
          fill = neg_log10_padj,
          size = setSize
        ),
        alpha = 0.95,
        stroke = 1
      ) +
      facet_grid(direction_label ~ ., scales = "free_y") +
      scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
      scale_shape_manual(values = comparison_shapes, name = NULL) +
      scale_color_manual(values = comparison_outline_colors, name = NULL) +
      scale_fill_viridis_c(
        name = expression(-log[10]~adjusted~italic(P)),
        option = "magma"
      ) +
      coord_cartesian(xlim = c(-overlap_nes_limit, overlap_nes_limit)) +
      labs(
        title = paste0("Overlapping ", ontology_label, " pathways"),
        subtitle = paste0("Pathways shown pass adjusted P <= ", gsea_padj_cutoff, " in both contrasts and the same direction"),
        x = "Normalized enrichment score",
        y = NULL,
        size = "Gene set size"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8, lineheight = 0.9),
        panel.grid.major.y = element_blank(),
        legend.position = "top",
        legend.box = "horizontal",
        legend.justification = "left",
        plot.margin = margin(5.5, 5.5, 5.5, 28)
      )
  } else {
    overlap_plot = ggdraw() +
      draw_label(
        paste0("No same-direction ", ontology_label, " overlaps at adjusted P <= ", gsea_padj_cutoff),
        size = 12
      )
  }
  
  n_direction_panels = max(1, n_distinct(overlap_plot_terms$direction_label))
  output_height = max(7, 2.4 * n_direction_panels + 0.18 * nrow(overlap_plot_terms))
  
  ggsave(
    output_file,
    overlap_plot,
    width = 12,
    height = output_height
  )
  
  overlap_plot
}

overlap_bp_plot = make_overlap_plot(
  ontology_code = "BP",
  ontology_label = "GO Biological Process",
  output_file = here("plots/go bp gsea overlapping pathways young CD vs old CD and old ID vs old CD.pdf")
)

overlap_cc_plot = make_overlap_plot(
  ontology_code = "CC",
  ontology_label = "GO Cellular Component",
  output_file = here("plots/go cc gsea overlapping pathways young CD vs old CD and old ID vs old CD.pdf")
)

overlap_bp_plot
overlap_cc_plot
