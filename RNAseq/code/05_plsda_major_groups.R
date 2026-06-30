library(here)
library(tidyverse)
library(edgeR)
library(mixOmics)
library(ComplexHeatmap)

here::i_am("code/05_plsda_major_groups.R")


# ========== 0.0 - Analysis settings ==========
# PLS-DA is run on TMM-normalized log2 CPM values from gene-level counts.
# This is preferable to TPM here because the available RNAseq inputs are count
# matrices, and TMM logCPM is appropriate for between-sample ordination.

top_vip_genes_to_plot = 100
plsda_ncomp = 2

group_levels = c("Young CD", "Old CD", "Old ID")
group_colors = c(
  "Young CD" = "#2C7FB8",
  "Old CD" = "#4D9221",
  "Old ID" = "#D95F02"
)

dir.create(here("data/processed/plsda"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("plots/plsda"), recursive = TRUE, showWarnings = FALSE)


# ========== 1.0 - Read counts and metadata ==========
meta = read.csv(here("data/metadata/metadata.csv")) %>%
  mutate(
    sample = paste0("r", sample),
    cohort = factor(cohort),
    age.grp = factor(age.grp, levels = c("young", "old")),
    diet = factor(diet, levels = c("CD", "ID")),
    analysis_group = case_when(
      age.grp == "young" & diet == "CD" ~ "Young CD",
      age.grp == "old" & diet == "CD" ~ "Old CD",
      age.grp == "old" & diet == "ID" ~ "Old ID",
      TRUE ~ paste(age.grp, diet)
    ),
    analysis_group = factor(analysis_group, levels = group_levels)
  ) %>%
  dplyr::filter(!is.na(analysis_group)) %>%
  arrange(analysis_group, cohort, sample)

counts = read.csv(here("data/processed/counts/counts matrix with mouse gene symbols.csv")) %>%
  column_to_rownames(var = "gene_id")

counts = counts[, meta$sample]
stopifnot(identical(colnames(counts), meta$sample))


# ========== 2.0 - TMM-normalized log2 CPM expression matrix ==========
dge = DGEList(counts = counts, group = meta$analysis_group)

keep_genes = filterByExpr(dge, group = meta$analysis_group)
dge = dge[keep_genes, , keep.lib.sizes = FALSE]
dge = calcNormFactors(dge, method = "TMM")

logcpm = cpm(dge, log = TRUE, prior.count = 1)

# mixOmics expects samples in rows and features in columns.
plsda_x = t(logcpm)
plsda_y = meta$analysis_group


# ========== 3.0 - PLS-DA across major study groups ==========
set.seed(20260626)
plsda_fit = plsda(
  X = plsda_x,
  Y = plsda_y,
  ncomp = plsda_ncomp,
  scale = TRUE
)

plsda_scores = as.data.frame(plsda_fit$variates$X) %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  dplyr::rename(
    PLSDA1 = comp1,
    PLSDA2 = comp2
  ) %>%
  left_join(meta, by = "sample") %>%
  arrange(analysis_group, cohort, sample)

write.csv(
  plsda_scores,
  here("data/processed/plsda/plsda major groups sample scores.csv"),
  row.names = FALSE
)

component_variance = round(100 * plsda_fit$prop_expl_var$X[seq_len(plsda_ncomp)], 1)

plsda_score_plot = ggplot(
  plsda_scores,
  aes(PLSDA1, PLSDA2, color = analysis_group)
) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  stat_ellipse(type = "norm", linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 3.4, alpha = 0.95) +
  ggrepel::geom_text_repel(
    aes(label = paste0(sample, "\nC", cohort)),
    color = "grey20",
    size = 3,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    segment.alpha = 0.45,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors, name = NULL) +
  labs(
    title = "PLS-DA: major RNA-seq groups",
    subtitle = "TMM-normalized log2 CPM expression",
    x = paste0("PLS-DA component 1 (", component_variance[1], "% X variance)"),
    y = paste0("PLS-DA component 2 (", component_variance[2], "% X variance)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

plsda_score_plot

ggsave(
  here("plots/plsda/plsda major groups ordination.pdf"),
  plsda_score_plot,
  width = 7.6,
  height = 6.2
)


# ========== 4.0 - VIP-ranked genes ==========
vip_scores = mixOmics::vip(plsda_fit) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  as_tibble() %>%
  dplyr::rename(
    VIP_component_1 = comp1,
    VIP_component_2 = comp2
  ) %>%
  mutate(
    VIP_max_component_1_2 = pmax(VIP_component_1, VIP_component_2, na.rm = TRUE),
    VIP_mean_component_1_2 = rowMeans(dplyr::select(., VIP_component_1, VIP_component_2), na.rm = TRUE)
  ) %>%
  arrange(desc(VIP_max_component_1_2), desc(VIP_mean_component_1_2), gene_id)

write.csv(
  vip_scores,
  here("data/processed/plsda/plsda major groups VIP scores.csv"),
  row.names = FALSE
)

top_vip_genes = vip_scores %>%
  slice_head(n = top_vip_genes_to_plot)

write.csv(
  top_vip_genes,
  here("data/processed/plsda/plsda major groups top 100 VIP genes.csv"),
  row.names = FALSE
)


# ========== 5.0 - Top VIP gene expression heatmap ==========
top_vip_logcpm = logcpm[top_vip_genes$gene_id, meta$sample]

top_vip_row_z = t(scale(t(top_vip_logcpm)))
top_vip_row_z[is.na(top_vip_row_z)] = 0

heatmap_column_labels = meta$sample %>%
  str_remove("^r") %>%
  str_remove("_mus$")

colnames(top_vip_row_z) = heatmap_column_labels

heatmap_column_split = factor(meta$analysis_group, levels = group_levels)

top_vip_heatmap = Heatmap(
  top_vip_row_z,
  name = "row z-score",
  col = circlize::colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "white", "#B2182B")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = heatmap_column_split,
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
  here("plots/plsda/plsda major groups top 100 VIP gene heatmap.pdf"),
  width = 9.5,
  height = 12
)
draw(
  top_vip_heatmap,
  heatmap_legend_side = "right"
)
dev.off()

top_vip_genes
