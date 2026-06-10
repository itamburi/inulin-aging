library(here)
library(tidyverse)
library(cowplot)
library(ggforce)

here::i_am("code/03 DEG overlap Venn diagrams.R")


# ========== 0.0 - Analysis settings ==========
# This script compares significant DEGs from:
#   1. old ID vs old CD: positive log2FoldChange means higher in old ID.
#   2. young CD vs old CD: positive log2FoldChange means higher in young CD.
#
# The nominal P cutoff matches the existing old ID vs CD DEG workflow in
# code/01 CDvID deseq2.R.

deg_pvalue_cutoff = 0.05


# ========== 1.0 - Read saved DESeq2 results ==========
id_cd_results = read.csv(here("data/processed/deseq2 ID vs CD old only cohort adjusted.csv")) %>%
  as_tibble()

age_cd_results = read.csv(here("data/processed/deseq2 young CD vs old CD.csv")) %>%
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
  here("data/processed/deseq2 DEG venn overlap summary.csv"),
  row.names = FALSE
)

write.csv(
  deg_venn_gene_lists,
  here("data/processed/deseq2 DEG venn overlap gene lists.csv"),
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
  here("plots/deseq2 DEG venn overlap up down old ID vs CD and young CD vs old CD.pdf"),
  deg_venn_figure,
  width = 9.4,
  height = 4.9,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  here("plots/deseq2 DEG venn overlap up down old ID vs CD and young CD vs old CD.png"),
  deg_venn_figure,
  width = 9.4,
  height = 4.9,
  dpi = 600,
  bg = "white"
)

deg_venn_summary
