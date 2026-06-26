library(here)
library(tidyverse)
library(cowplot)
library(ggforce)
library(edgeR)
library(ComplexHeatmap)

here::i_am("code/03_deg_overlap_and_youth_shift.R")


# ========== 0.0 - Analysis settings ==========
# This script compares significant DEGs from:
#   1. old ID vs old CD: positive log2FoldChange means higher in old ID.
#   2. young CD vs old CD: positive log2FoldChange means higher in young CD.
#
# The nominal P cutoff matches the existing old ID vs CD DEG workflow in
# code/01_prepare_counts_and_diet_deseq2.R.

deg_pvalue_cutoff = 0.05
min_abs_age_log2fc_for_rescue_fraction = 0.10
top_youth_shift_genes_to_table = 100
top_youth_shift_genes_to_plot = 50


# ========== 1.0 - Read saved DESeq2 results ==========
id_cd_results = read.csv(here("data/processed/deseq2/deseq2 ID vs CD old only cohort adjusted.csv")) %>%
  as_tibble()

age_cd_results = read.csv(here("data/processed/deseq2/deseq2 young CD vs old CD.csv")) %>%
  as_tibble()


# ========== 2.0 - Define significant up/down DEG sets ==========
id_cd_sig = id_cd_results %>%
  dplyr::filter(
    !is.na(pvalue),
    pvalue <= deg_pvalue_cutoff,
    !is.na(log2FoldChange),
    log2FoldChange != 0
  ) %>%
  mutate(
    contrast = "old ID vs old CD",
    direction = ifelse(log2FoldChange > 0, "Up", "Down")
  )

age_cd_sig = age_cd_results %>%
  dplyr::filter(
    !is.na(pvalue),
    pvalue <= deg_pvalue_cutoff,
    !is.na(log2FoldChange),
    log2FoldChange != 0
  ) %>%
  mutate(
    contrast = "young CD vs old CD",
    direction = ifelse(log2FoldChange > 0, "Up", "Down")
  )

id_cd_up_genes = id_cd_sig %>%
  dplyr::filter(direction == "Up") %>%
  pull(gene.id) %>%
  unique()

id_cd_down_genes = id_cd_sig %>%
  dplyr::filter(direction == "Down") %>%
  pull(gene.id) %>%
  unique()

age_cd_up_genes = age_cd_sig %>%
  dplyr::filter(direction == "Up") %>%
  pull(gene.id) %>%
  unique()

age_cd_down_genes = age_cd_sig %>%
  dplyr::filter(direction == "Down") %>%
  pull(gene.id) %>%
  unique()


# ========== 3.0 - Count shared and distinct DEG sets ==========
deg_venn_summary = tibble(
  direction = c("Upregulated", "Downregulated"),
  id_cd_label = c("Old ID > old CD", "Old ID < old CD"),
  age_cd_label = c("Young CD > old CD", "Young CD < old CD"),
  id_cd_total = c(length(id_cd_up_genes), length(id_cd_down_genes)),
  age_cd_total = c(length(age_cd_up_genes), length(age_cd_down_genes)),
  id_cd_only = c(
    length(setdiff(id_cd_up_genes, age_cd_up_genes)),
    length(setdiff(id_cd_down_genes, age_cd_down_genes))
  ),
  shared = c(
    length(intersect(id_cd_up_genes, age_cd_up_genes)),
    length(intersect(id_cd_down_genes, age_cd_down_genes))
  ),
  age_cd_only = c(
    length(setdiff(age_cd_up_genes, id_cd_up_genes)),
    length(setdiff(age_cd_down_genes, id_cd_down_genes))
  )
)

deg_venn_gene_lists = bind_rows(
  tibble(
    direction = "Upregulated",
    overlap_group = "Old ID vs old CD only",
    gene.id = setdiff(id_cd_up_genes, age_cd_up_genes)
  ),
  tibble(
    direction = "Upregulated",
    overlap_group = "Shared",
    gene.id = intersect(id_cd_up_genes, age_cd_up_genes)
  ),
  tibble(
    direction = "Upregulated",
    overlap_group = "Young CD vs old CD only",
    gene.id = setdiff(age_cd_up_genes, id_cd_up_genes)
  ),
  tibble(
    direction = "Downregulated",
    overlap_group = "Old ID vs old CD only",
    gene.id = setdiff(id_cd_down_genes, age_cd_down_genes)
  ),
  tibble(
    direction = "Downregulated",
    overlap_group = "Shared",
    gene.id = intersect(id_cd_down_genes, age_cd_down_genes)
  ),
  tibble(
    direction = "Downregulated",
    overlap_group = "Young CD vs old CD only",
    gene.id = setdiff(age_cd_down_genes, id_cd_down_genes)
  )
) %>%
  arrange(
    factor(direction, levels = c("Upregulated", "Downregulated")),
    factor(
      overlap_group,
      levels = c("Old ID vs old CD only", "Shared", "Young CD vs old CD only")
    ),
    gene.id
  )

write.csv(
  deg_venn_summary,
  here("data/processed/overlaps/deseq2 DEG venn overlap summary.csv"),
  row.names = FALSE
)

write.csv(
  deg_venn_gene_lists,
  here("data/processed/overlaps/deseq2 DEG venn overlap gene lists.csv"),
  row.names = FALSE
)


# ========== 4.0 - Figure-ready Venn diagrams ==========
venn_colors = c(
  id_cd = "#2C7FB8",
  age_cd = "#238B45"
)

make_venn_panel = function(summary_row) {
  circle_data = tibble(
    set = c(summary_row$id_cd_label, summary_row$age_cd_label),
    x0 = c(-0.72, 0.72),
    y0 = c(0, 0),
    r = c(1.12, 1.12),
    fill = c(venn_colors["id_cd"], venn_colors["age_cd"])
  )
  
  count_labels = tibble(
    label = c(summary_row$id_cd_only, summary_row$shared, summary_row$age_cd_only),
    x = c(-1.05, 0, 1.05),
    y = c(0, 0, 0)
  )
  
  set_labels = tibble(
    label = c(
      paste0(summary_row$id_cd_label, "\n", "n = ", summary_row$id_cd_total),
      paste0(summary_row$age_cd_label, "\n", "n = ", summary_row$age_cd_total)
    ),
    x = c(-1.28, 1.28),
    y = c(-1.42, -1.42),
    color = c(venn_colors["id_cd"], venn_colors["age_cd"])
  )
  
  ggplot() +
    geom_circle(
      data = circle_data,
      aes(x0 = x0, y0 = y0, r = r, fill = fill),
      color = "grey20",
      linewidth = 0.5,
      alpha = 0.24,
      show.legend = FALSE
    ) +
    geom_text(
      data = count_labels,
      aes(x, y, label = label),
      size = 7.2,
      fontface = "bold",
      color = "black",
      family = "Helvetica"
    ) +
    geom_text(
      data = set_labels,
      aes(x, y, label = label, color = color),
      size = 4.4,
      fontface = "bold",
      lineheight = 0.95,
      family = "Helvetica",
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    coord_equal(xlim = c(-2.05, 2.05), ylim = c(-1.75, 1.28), clip = "off") +
    labs(
      title = summary_row$direction
    ) +
  theme_void(base_size = 12, base_family = "Helvetica") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5,
        color = "black",
        margin = margin(0, 0, 6, 0)
      ),
      plot.margin = margin(4, 8, 8, 8)
    )
}

up_venn_plot = make_venn_panel(deg_venn_summary[1, ])
down_venn_plot = make_venn_panel(deg_venn_summary[2, ])

deg_venn_title = ggdraw() +
  draw_label(
    "Overlap of significant DEGs",
    x = 0.5,
    y = 0.70,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 18,
    fontfamily = "Helvetica"
  ) +
  draw_label(
    paste0("Nominal P <= ", deg_pvalue_cutoff),
    x = 0.5,
    y = 0.25,
    hjust = 0.5,
    vjust = 0.5,
    size = 11,
    fontfamily = "Helvetica"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

deg_venn_panels = plot_grid(
  up_venn_plot,
  down_venn_plot,
  nrow = 1
)

deg_venn_figure = plot_grid(
  deg_venn_title,
  deg_venn_panels,
  ncol = 1,
  rel_heights = c(0.15, 1)
) +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  here("plots/overlaps/deseq2 DEG venn overlap up down old ID vs CD and young CD vs old CD.pdf"),
  deg_venn_figure,
  width = 9.4,
  height = 4.9,
  device = cairo_pdf,
  bg = "white"
)

deg_venn_summary


# ========== 5.0 - Youth-like old ID shifts relative to old CD ==========
# This section asks which genes move in the same direction in:
#   1. young CD vs old CD: the youth-associated expression difference.
#   2. old ID vs old CD: the inulin-associated shift in old mice.
#
# Same-direction genes are interpreted as old ID shifting toward the young CD
# state relative to old CD. Opposite-direction genes move away from young CD.

age_contrast = age_cd_results %>%
  dplyr::select(
    gene.id,
    age_baseMean = baseMean,
    age_log2FoldChange = log2FoldChange,
    age_lfcSE = lfcSE,
    age_stat = stat,
    age_pvalue = pvalue,
    age_padj = padj
  )

id_contrast = id_cd_results %>%
  dplyr::select(
    gene.id,
    id_baseMean = baseMean,
    id_log2FoldChange = log2FoldChange,
    id_lfcSE = lfcSE,
    id_stat = stat,
    id_pvalue = pvalue,
    id_padj = padj
  )

youth_shift_gene_scores = age_contrast %>%
  inner_join(id_contrast, by = "gene.id") %>%
  dplyr::filter(
    !is.na(age_log2FoldChange),
    !is.na(id_log2FoldChange),
    !is.na(age_pvalue),
    !is.na(id_pvalue)
  ) %>%
  mutate(
    age_direction = case_when(
      age_log2FoldChange > 0 ~ "Higher in young CD",
      age_log2FoldChange < 0 ~ "Lower in young CD",
      TRUE ~ "No age direction"
    ),
    id_direction = case_when(
      id_log2FoldChange > 0 ~ "Higher in old ID",
      id_log2FoldChange < 0 ~ "Lower in old ID",
      TRUE ~ "No ID direction"
    ),
    same_direction = age_log2FoldChange * id_log2FoldChange > 0,
    direction_class = case_when(
      age_log2FoldChange > 0 & id_log2FoldChange > 0 ~ "Youth-like higher",
      age_log2FoldChange < 0 & id_log2FoldChange < 0 ~ "Youth-like lower",
      age_log2FoldChange > 0 & id_log2FoldChange < 0 ~ "Opposing: ID lower",
      age_log2FoldChange < 0 & id_log2FoldChange > 0 ~ "Opposing: ID higher",
      TRUE ~ "No clear direction"
    ),
    direction_class = factor(
      direction_class,
      levels = c(
        "Youth-like higher",
        "Youth-like lower",
        "Opposing: ID lower",
        "Opposing: ID higher",
        "No clear direction"
      )
    ),
    rescue_fraction = ifelse(
      abs(age_log2FoldChange) >= min_abs_age_log2fc_for_rescue_fraction,
      id_log2FoldChange / age_log2FoldChange,
      NA_real_
    ),
    rescue_interpretation = case_when(
      !same_direction ~ "Opposing direction",
      is.na(rescue_fraction) ~ "Age contrast too small",
      rescue_fraction < 0.25 ~ "Weak youth-like shift",
      rescue_fraction <= 1.50 ~ "Partial/full youth-like shift",
      rescue_fraction > 1.50 ~ "Overshoot past young CD",
      TRUE ~ "Unclassified"
    ),
    age_nominal = age_pvalue <= deg_pvalue_cutoff,
    id_nominal = id_pvalue <= deg_pvalue_cutoff,
    significance_class = case_when(
      age_nominal & id_nominal ~ "Nominal in both",
      age_nominal & !id_nominal ~ "Age contrast only",
      !age_nominal & id_nominal ~ "ID contrast only",
      TRUE ~ "Neither nominal"
    ),
    limiting_pvalue = pmax(age_pvalue, id_pvalue, .Machine$double.xmin),
    limiting_neg_log10_p = -log10(limiting_pvalue),
    shared_effect_size = pmin(abs(age_log2FoldChange), abs(id_log2FoldChange)),
    youth_shift_score = ifelse(
      same_direction,
      shared_effect_size * limiting_neg_log10_p,
      -shared_effect_size * limiting_neg_log10_p
    )
  ) %>%
  arrange(desc(youth_shift_score), gene.id)

youth_shift_summary = youth_shift_gene_scores %>%
  count(direction_class, significance_class, name = "n_genes") %>%
  arrange(direction_class, significance_class)

top_youth_shift_genes = youth_shift_gene_scores %>%
  dplyr::filter(same_direction) %>%
  arrange(desc(youth_shift_score), gene.id) %>%
  slice_head(n = top_youth_shift_genes_to_table)

write.csv(
  youth_shift_gene_scores,
  here("data/processed/youth_shift/deseq2 youth-like old ID shift gene scores.csv"),
  row.names = FALSE
)

write.csv(
  youth_shift_summary,
  here("data/processed/youth_shift/deseq2 youth-like old ID shift summary.csv"),
  row.names = FALSE
)

write.csv(
  top_youth_shift_genes,
  here("data/processed/youth_shift/deseq2 youth-like old ID shift top genes.csv"),
  row.names = FALSE
)


# ========== 6.0 - Youth-shift concordance plot ==========
youth_shift_plot_data = youth_shift_gene_scores %>%
  mutate(
    highlight_group = case_when(
      same_direction & age_nominal & id_nominal ~ "Youth-like, nominal in both",
      same_direction ~ "Youth-like direction",
      !same_direction & (age_nominal | id_nominal) ~ "Opposing, nominal in at least one",
      TRUE ~ "Other"
    ),
    highlight_group = factor(
      highlight_group,
      levels = c(
        "Youth-like, nominal in both",
        "Youth-like direction",
        "Opposing, nominal in at least one",
        "Other"
      )
    )
  )

youth_shift_label_genes = youth_shift_plot_data %>%
  dplyr::filter(same_direction) %>%
  arrange(desc(youth_shift_score), gene.id) %>%
  slice_head(n = 30)

youth_shift_colors = c(
  "Youth-like, nominal in both" = "#D95F02",
  "Youth-like direction" = "#7570B3",
  "Opposing, nominal in at least one" = "#1B9E77",
  "Other" = "grey78"
)

youth_shift_scatter = ggplot(
  youth_shift_plot_data,
  aes(age_log2FoldChange, id_log2FoldChange)
) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", color = "grey40", linewidth = 0.45) +
  geom_point(aes(color = highlight_group), alpha = 0.62, size = 1.15, stroke = 0) +
  ggrepel::geom_text_repel(
    data = youth_shift_label_genes,
    aes(label = gene.id),
    color = "black",
    size = 3,
    box.padding = 0.45,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.alpha = 0.45,
    show.legend = FALSE
  ) +
  scale_color_manual(values = youth_shift_colors, name = NULL) +
  labs(
    title = "DESeq2 contrast concordance",
    subtitle = "Genes in upper-right or lower-left shift old ID toward young CD relative to old CD",
    x = expression(log[2]~fold~change~"(young CD / old CD)"),
    y = expression(log[2]~fold~change~"(old ID / old CD)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

youth_shift_scatter

ggsave(
  here("plots/youth_shift/deseq2 youth-like old ID shift concordance.pdf"),
  youth_shift_scatter,
  width = 8,
  height = 6.5
)


# ========== 7.0 - Heatmap of top youth-like shift genes ==========
meta = read.csv(here("data/metadata/metadata.csv")) %>%
  mutate(
    sample = paste0("r", sample),
    cohort = factor(cohort),
    age.grp = factor(age.grp, levels = c("young", "old")),
    diet = factor(diet, levels = c("CD", "ID")),
    analysis_group = case_when(
      age.grp == "old" & diet == "CD" ~ "Old CD",
      age.grp == "old" & diet == "ID" ~ "Old ID",
      age.grp == "young" & diet == "CD" ~ "Young CD",
      TRUE ~ paste(age.grp, diet)
    ),
    analysis_group = factor(analysis_group, levels = c("Old CD", "Old ID", "Young CD"))
  ) %>%
  dplyr::filter(!is.na(analysis_group)) %>%
  arrange(analysis_group, cohort, sample)

counts = read.csv(here("data/processed/counts/counts matrix with mouse gene symbols.csv")) %>%
  column_to_rownames(var = "gene_id")

counts = counts[, meta$sample]
stopifnot(identical(colnames(counts), meta$sample))

dge = DGEList(counts = counts, group = meta$analysis_group)
keep_genes = filterByExpr(dge, group = meta$analysis_group)
dge = dge[keep_genes, , keep.lib.sizes = FALSE]
dge = calcNormFactors(dge, method = "TMM")
logcpm = cpm(dge, log = TRUE, prior.count = 1)

top_youth_shift_genes_for_heatmap = top_youth_shift_genes %>%
  slice_head(n = top_youth_shift_genes_to_plot)

top_youth_shift_logcpm = logcpm[top_youth_shift_genes_for_heatmap$gene.id, meta$sample]

top_youth_shift_row_z = t(scale(t(top_youth_shift_logcpm)))
top_youth_shift_row_z[is.na(top_youth_shift_row_z)] = 0

row_split = top_youth_shift_genes_for_heatmap$direction_class
row_split = factor(row_split, levels = c("Youth-like higher", "Youth-like lower"))

heatmap_column_labels = meta$sample %>%
  str_remove("^r") %>%
  str_remove("_mus$")

colnames(top_youth_shift_row_z) = heatmap_column_labels

youth_shift_heatmap = Heatmap(
  top_youth_shift_row_z,
  name = "row z-score",
  col = circlize::colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "white", "#B2182B")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_split = row_split,
  column_split = meta$analysis_group,
  column_title = "%s",
  column_gap = grid::unit(3, "mm"),
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 7.5),
  row_names_gp = grid::gpar(fontsize = 7),
  border = FALSE,
  heatmap_legend_param = list(title = "Row z-score")
)

pdf(
  here("plots/youth_shift/deseq2 youth-like old ID shift top gene heatmap.pdf"),
  width = 9.5,
  height = 9
)
draw(
  youth_shift_heatmap,
  heatmap_legend_side = "right"
)
dev.off()

youth_shift_summary
