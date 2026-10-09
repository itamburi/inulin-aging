library(here)
library(tidyverse)
library(edgeR)
library(DESeq2)

here::i_am("code/08_plot_selected_gene_expression.R")

plot_objects = readRDS(here("data/processed/counts/rnaseq filtered counts qc objects.rds"))

selected_genes = c("Ttr", "Serpina1c", "Ahsg", "Lep")
group_levels = c("Old CD", "Old ID", "Young CD")
plot_orders = list(
  "old CD first" = c("Old CD", "Old ID", "Young CD"),
  "young CD first" = c("Young CD", "Old CD", "Old ID")
)
group_axis_labels = c(
  "Old CD" = "Old\nCD",
  "Old ID" = "Old\nID",
  "Young CD" = "Young\nCD"
)

counts = plot_objects$filtered_counts
meta = plot_objects$meta

stopifnot(all(selected_genes %in% rownames(counts)))
stopifnot(all(colnames(counts) %in% rownames(meta)))

meta = meta[colnames(counts), , drop = FALSE] %>%
  as.data.frame() %>%
  mutate(
    age.grp = factor(age.grp, levels = c("old", "young")),
    diet = factor(diet, levels = c("CD", "ID"))
  )

sample_groups = meta %>%
  rownames_to_column("sample") %>%
  mutate(
    group = case_when(
      age.grp == "old" & diet == "CD" ~ "Old CD",
      age.grp == "old" & diet == "ID" ~ "Old ID",
      age.grp == "young" & diet == "CD" ~ "Young CD",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(sample, cohort, age.grp, diet, group)

dge = DGEList(counts = round(counts))
dge = calcNormFactors(dge)
cpm_values = cpm(dge, normalized.lib.sizes = TRUE)

dds = DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = meta,
  design = ~ 1
)
dds = estimateSizeFactors(dds)
vsd = vst(dds, blind = TRUE)

cpm_data = cpm_values[selected_genes, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "cpm"
  ) %>%
  left_join(sample_groups, by = "sample") %>%
  dplyr::filter(!is.na(group)) %>%
  mutate(
    gene = factor(gene, levels = selected_genes),
    group = factor(group, levels = group_levels)
  )

vst_data = SummarizedExperiment::assay(vsd)[selected_genes, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "vst_expression"
  ) %>%
  left_join(sample_groups, by = "sample") %>%
  dplyr::filter(!is.na(group)) %>%
  mutate(
    gene = factor(gene, levels = selected_genes),
    group = factor(group, levels = group_levels)
  )

diet_deseq = read_csv(
  here("data/processed/deseq2/deseq2 ID vs CD old only cohort adjusted.csv"),
  show_col_types = FALSE
) %>%
  dplyr::filter(gene.id %in% selected_genes) %>%
  transmute(
    gene = gene.id,
    group = "Old ID",
    comparison,
    log2_fold_change_vs_old_cd = log2FoldChange,
    lfcSE,
    pvalue,
    padj
  )

age_deseq = read_csv(
  here("data/processed/deseq2/deseq2 young CD vs old CD.csv"),
  show_col_types = FALSE
) %>%
  dplyr::filter(gene.id %in% selected_genes) %>%
  transmute(
    gene = gene.id,
    group = "Young CD",
    comparison,
    log2_fold_change_vs_old_cd = log2FoldChange,
    lfcSE,
    pvalue,
    padj
  )

old_id_young_cd_deseq = read_csv(
  here("data/processed/deseq2/deseq2 old ID vs young CD.csv"),
  show_col_types = FALSE
) %>%
  dplyr::filter(gene.id %in% selected_genes) %>%
  transmute(
    gene = gene.id,
    group = "Old ID",
    comparison,
    log2_fold_change_vs_young_cd = log2FoldChange,
    lfcSE,
    pvalue,
    padj
  )

fold_change_data = bind_rows(
  tibble(
    gene = selected_genes,
    group = "Old CD",
    comparison = "reference",
    log2_fold_change_vs_old_cd = 0,
    lfcSE = NA_real_,
    pvalue = NA_real_,
    padj = NA_real_
  ),
  diet_deseq,
  age_deseq
) %>%
  mutate(
    gene = factor(gene, levels = selected_genes),
    group = factor(group, levels = group_levels),
    fold_change_vs_old_cd = 2^log2_fold_change_vs_old_cd,
    fold_change_lower_1se = ifelse(
      is.na(lfcSE),
      NA_real_,
      2^(log2_fold_change_vs_old_cd - lfcSE)
    ),
    fold_change_upper_1se = ifelse(
      is.na(lfcSE),
      NA_real_,
      2^(log2_fold_change_vs_old_cd + lfcSE)
    )
  )

fold_change_young_cd_data = bind_rows(
  age_deseq %>%
    transmute(
      gene,
      group = "Old CD",
      comparison = "old CD vs young CD",
      log2_fold_change_vs_young_cd = -log2_fold_change_vs_old_cd,
      lfcSE,
      pvalue,
      padj
    ),
  old_id_young_cd_deseq,
  tibble(
    gene = selected_genes,
    group = "Young CD",
    comparison = "reference",
    log2_fold_change_vs_young_cd = 0,
    lfcSE = NA_real_,
    pvalue = NA_real_,
    padj = NA_real_
  )
) %>%
  mutate(
    gene = factor(gene, levels = selected_genes),
    group = factor(group, levels = group_levels),
    fold_change_vs_young_cd = 2^log2_fold_change_vs_young_cd,
    fold_change_lower_1se = ifelse(
      is.na(lfcSE),
      NA_real_,
      2^(log2_fold_change_vs_young_cd - lfcSE)
    ),
    fold_change_upper_1se = ifelse(
      is.na(lfcSE),
      NA_real_,
      2^(log2_fold_change_vs_young_cd + lfcSE)
    )
  )

cpm_summary = cpm_data %>%
  group_by(gene, group) %>%
  summarise(
    n = n(),
    mean_cpm = mean(cpm),
    sem_cpm = sd(cpm) / sqrt(n),
    .groups = "drop"
  )

vst_summary = vst_data %>%
  group_by(gene, group) %>%
  summarise(
    n = n(),
    mean_vst_expression = mean(vst_expression),
    sem_vst_expression = sd(vst_expression) / sqrt(n),
    .groups = "drop"
  )

dir.create(here("data/processed/expression"), recursive = TRUE, showWarnings = FALSE)
write_csv(cpm_data, here("data/processed/expression/selected gene CPM Ttr Serpina1c Ahsg Lep sample values.csv"))
write_csv(cpm_summary, here("data/processed/expression/selected gene CPM Ttr Serpina1c Ahsg Lep summary.csv"))
write_csv(vst_data, here("data/processed/expression/selected gene VST Ttr Serpina1c Ahsg Lep sample values.csv"))
write_csv(vst_summary, here("data/processed/expression/selected gene VST Ttr Serpina1c Ahsg Lep summary.csv"))
write_csv(fold_change_data, here("data/processed/expression/selected gene expression Ttr Serpina1c Ahsg Lep deseq2 fold change values.csv"))
write_csv(fold_change_young_cd_data, here("data/processed/expression/selected gene expression Ttr Serpina1c Ahsg Lep deseq2 fold change values vs young CD.csv"))

group_colors = c("Old CD" = "grey80", "Old ID" = "#90EE90", "Young CD" = "#FF638E")
group_point_colors = c("Old CD" = "grey35", "Old ID" = "#2E8B57", "Young CD" = "#C2185B")

theme_selected_gene = function() {
  theme_classic(base_size = 12, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 15, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.y = element_text(color = "black", size = 14),
      panel.spacing = unit(12, "mm"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

p_to_stars = function(pvalue) {
  case_when(
    pvalue < 0.001 ~ "***",
    pvalue < 0.01 ~ "**",
    pvalue < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

make_fc_brackets = function(data, reference_group, comparison_groups, plot_levels,
                            value_col, lower_col, upper_col, log2_scale = FALSE) {
  value_col = rlang::ensym(value_col)
  lower_col = rlang::ensym(lower_col)
  upper_col = rlang::ensym(upper_col)

  bracket_data = data %>%
    dplyr::filter(as.character(group) %in% comparison_groups) %>%
    group_by(gene) %>%
    mutate(
      y_top = max(!!upper_col, !!value_col, ifelse(log2_scale, 0, 1), na.rm = TRUE),
      y_bottom = min(!!lower_col, !!value_col, ifelse(log2_scale, 0, 1), na.rm = TRUE),
      y_range = y_top - y_bottom,
      y_range = ifelse(y_range == 0, 1, y_range),
      bracket_index = match(as.character(group), comparison_groups),
      y = if (log2_scale) {
        y_top + y_range * (0.14 + 0.30 * (bracket_index - 1))
      } else {
        y_top * (1.16 + 0.36 * (bracket_index - 1))
      },
      y_tip = if (log2_scale) {
        y - y_range * 0.035
      } else {
        y - y_top * 0.035
      },
      x_start = match(reference_group, plot_levels),
      x_end = match(as.character(group), plot_levels),
      x_mid = (x_start + x_end) / 2,
      p_label = p_to_stars(pvalue)
    ) %>%
    ungroup()

  list(
    labels = bracket_data,
    segments = bind_rows(
      bracket_data %>% transmute(gene, x = x_start, xend = x_end, y = y, yend = y),
      bracket_data %>% transmute(gene, x = x_start, xend = x_start, y = y_tip, yend = y),
      bracket_data %>% transmute(gene, x = x_end, xend = x_end, y = y_tip, yend = y)
    )
  )
}

make_expression_plot = function(summary_data, sample_data, mean_col, sem_col, sample_col,
                                 y_label, plot_levels) {
  mean_col = rlang::ensym(mean_col)
  sem_col = rlang::ensym(sem_col)
  sample_col = rlang::ensym(sample_col)

  ggplot(summary_data, aes(group, !!mean_col, fill = group)) +
    geom_col(width = 0.84, color = "grey15", linewidth = 0.75) +
    geom_errorbar(aes(ymin = !!mean_col - !!sem_col, ymax = !!mean_col + !!sem_col),
                  width = 0.2, linewidth = 0.55, color = "grey15") +
    geom_point(
      data = sample_data,
      aes(x = group, y = !!sample_col, color = group),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.12, height = 0, seed = 1),
      shape = 16,
      size = 2.4,
      alpha = 0.9
    ) +
    facet_wrap(~ gene, ncol = 2, scales = "free_y", axes = "all_x") +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_color_manual(values = group_point_colors, guide = "none") +
    scale_x_discrete(
      limits = plot_levels,
      labels = group_axis_labels[plot_levels],
      expand = expansion(add = 0.65)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = y_label) +
    theme_selected_gene()
}

make_fc_plot = function(data, value_col, lower_col, upper_col, reference_group,
                         comparison_groups, y_label, baseline, plot_levels,
                         log2_scale = FALSE) {
  value_col = rlang::ensym(value_col)
  lower_col = rlang::ensym(lower_col)
  upper_col = rlang::ensym(upper_col)

  brackets = make_fc_brackets(
    data = data,
    reference_group = reference_group,
    comparison_groups = comparison_groups,
    plot_levels = plot_levels,
    value_col = !!value_col,
    lower_col = !!lower_col,
    upper_col = !!upper_col,
    log2_scale = log2_scale
  )

  ggplot(data, aes(group, !!value_col, fill = group)) +
    geom_hline(yintercept = baseline, linetype = "22", color = "grey45", linewidth = 0.4) +
    geom_col(width = 0.84, color = "grey15", linewidth = 0.75) +
    geom_errorbar(
      aes(ymin = !!lower_col, ymax = !!upper_col),
      data = data %>% dplyr::filter(!is.na(lfcSE)),
      width = 0.2,
      linewidth = 0.55,
      color = "grey15"
    ) +
    geom_segment(
      data = brackets$segments,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.45,
      color = "black"
    ) +
    geom_text(
      data = brackets$labels,
      aes(x = x_mid, y = y, label = p_label),
      inherit.aes = FALSE,
      vjust = -0.25,
      size = 4.2,
      fontface = "bold",
      color = "black"
    ) +
    facet_wrap(~ gene, ncol = 2, scales = "free_y", axes = "all_x") +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_x_discrete(
      limits = plot_levels,
      labels = group_axis_labels[plot_levels],
      expand = expansion(add = 0.65)
    ) +
    scale_y_continuous(
      expand = expansion(mult = if (log2_scale) c(0.08, 0.36) else c(0, 0.36))
    ) +
    labs(x = NULL, y = y_label) +
    theme_selected_gene()
}

selected_gene_plot_dir = here("plots/diet DEx/selected gene expression")
expression_plot_dir = file.path(selected_gene_plot_dir, "expression values")
old_cd_reference_plot_dir = file.path(selected_gene_plot_dir, "old CD reference")
young_cd_reference_plot_dir = file.path(selected_gene_plot_dir, "young CD reference")

save_pdf = function(plot, dir_path, filename) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(dir_path, filename),
    plot,
    width = 4.32,
    height = 6.2,
    device = grDevices::cairo_pdf
  )
}

for (order_name in names(plot_orders)) {
  plot_levels = plot_orders[[order_name]]

  save_pdf(
    make_expression_plot(cpm_summary, cpm_data, mean_cpm, sem_cpm, cpm, "CPM", plot_levels),
    file.path(expression_plot_dir, order_name),
    paste0("selected gene CPM Ttr Serpina1c Ahsg Lep - ", order_name, ".pdf")
  )

  save_pdf(
    make_expression_plot(vst_summary, vst_data, mean_vst_expression, sem_vst_expression,
                         vst_expression, "VST expression", plot_levels),
    file.path(expression_plot_dir, order_name),
    paste0("selected gene VST Ttr Serpina1c Ahsg Lep - ", order_name, ".pdf")
  )

  save_pdf(
    make_fc_plot(
      fold_change_data,
      fold_change_vs_old_cd,
      fold_change_lower_1se,
      fold_change_upper_1se,
      "Old CD",
      c("Old ID", "Young CD"),
      "Fold change vs Old CD",
      1,
      plot_levels
    ),
    file.path(old_cd_reference_plot_dir, order_name),
    paste0("selected gene fold change vs old CD Ttr Serpina1c Ahsg Lep - ", order_name, ".pdf")
  )

  save_pdf(
    make_fc_plot(
      fold_change_data %>%
        mutate(
          log2_fold_change_lower_1se = log2_fold_change_vs_old_cd - lfcSE,
          log2_fold_change_upper_1se = log2_fold_change_vs_old_cd + lfcSE
        ),
      log2_fold_change_vs_old_cd,
      log2_fold_change_lower_1se,
      log2_fold_change_upper_1se,
      "Old CD",
      c("Old ID", "Young CD"),
      expression(log[2]*"FC vs Old CD"),
      0,
      plot_levels,
      log2_scale = TRUE
    ),
    file.path(old_cd_reference_plot_dir, order_name),
    paste0("selected gene log2 fold change vs old CD Ttr Serpina1c Ahsg Lep - ", order_name, ".pdf")
  )

  save_pdf(
    make_fc_plot(
      fold_change_young_cd_data,
      fold_change_vs_young_cd,
      fold_change_lower_1se,
      fold_change_upper_1se,
      "Young CD",
      c("Old CD", "Old ID"),
      "Fold change vs Young CD",
      1,
      plot_levels
    ),
    file.path(young_cd_reference_plot_dir, order_name),
    paste0("selected gene fold change vs young CD Ttr Serpina1c Ahsg Lep - ", order_name, ".pdf")
  )

  save_pdf(
    make_fc_plot(
      fold_change_young_cd_data %>%
        mutate(
          log2_fold_change_lower_1se = log2_fold_change_vs_young_cd - lfcSE,
          log2_fold_change_upper_1se = log2_fold_change_vs_young_cd + lfcSE
        ),
      log2_fold_change_vs_young_cd,
      log2_fold_change_lower_1se,
      log2_fold_change_upper_1se,
      "Young CD",
      c("Old CD", "Old ID"),
      expression(log[2]*"FC vs Young CD"),
      0,
      plot_levels,
      log2_scale = TRUE
    ),
    file.path(young_cd_reference_plot_dir, order_name),
    paste0("selected gene log2 fold change vs young CD Ttr Serpina1c Ahsg Lep - ", order_name, ".pdf")
  )
}
