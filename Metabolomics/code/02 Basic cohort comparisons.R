# ****** Inulin Aging Mouse Tissue (Liver and Muscle)
# ****** Spring 2026 - Ian T, Jessie K, Randa L



# 
library(tidyverse)
library(here)
library(broom)
library(ComplexHeatmap)
library(circlize)

here::i_am("code/02 Basic cohort comparisons.R")



# ========== 0.0 - Enviornment ==========
# -- start from the cleaned, long-format, normalized metabolomics table
expdata = readr::read_csv(
  here::here("data/processed/ion counts/tidy and normalized ion count data.csv"),
  show_col_types = FALSE
) 

sample_parts = stringr::str_match(expdata$sample, "^([124]M)_?0*([0-9]+)$")

expdata = expdata %>%
  mutate(
    sample.clean = dplyr::if_else(
      !is.na(sample_parts[, 1]),
      paste0(sample_parts[, 2], sample_parts[, 3]),
      sample
    ),
    cohort = stringr::str_extract(sample.clean, "^[0-9]+"),
    cohort = factor(cohort, levels = c("1", "2", "4"))
  )

dir.create(here::here("data/libraries"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("plots"), recursive = TRUE, showWarnings = FALSE)

metadata = tibble(
  sample = c(
    paste0("1M", c(1:4, 10:11)),
    paste0("1M", 5:9),
    paste0("2M", 16:20),
    paste0("2M", 22:33),
    paste0("4M", 1:5)
  ),
  diet = c(
    rep("CD", 6),
    rep("ID", 5),
    rep("CD", 5),
    rep("ID", 12),
    rep("CD", 5)
  ),
  cohort = c(
    rep("1", 11),
    rep("2", 17),
    rep("4", 5)
  )
) %>%
  mutate(
    animal = as.integer(stringr::str_extract(sample, "[0-9]+$")),
    age = dplyr::if_else(cohort == "4", "Young", "Old"),
    group = dplyr::case_when(
      cohort == "4" ~ "Young CD",
      diet == "CD" ~ "Old CD",
      diet == "ID" ~ "Old ID"
    ),
    cohort = factor(cohort, levels = c("1", "2", "4")),
    diet = factor(diet, levels = c("CD", "ID")),
    age = factor(age, levels = c("Young", "Old")),
    group = factor(group, levels = c("Young CD", "Old CD", "Old ID"))
  ) %>%
  arrange(cohort, sample)

readr::write_csv(metadata, here::here("data/libraries/metadata.csv"))

expdata = expdata %>%
  left_join(
    metadata %>%
      select(sample, diet, animal, age, group),
    by = c("sample.clean" = "sample")
  )


# ========== 1.0 - Old mouse diet comparison ==========
# -- compare old inulin diet vs control diet using sample-level median adjusted ion counts
old_mouse_sample_abundance = expdata %>%
  filter(age == "Old", !is.na(diet), !is.na(med.adj.ic)) %>%
  group_by(tissue, name, sample.clean, diet, cohort, animal) %>%
  summarise(
    med.adj.ic = median(med.adj.ic, na.rm = TRUE),
    n_obs = dplyr::n(),
    .groups = "drop"
  )

old_mouse_diet_ttest = old_mouse_sample_abundance %>%
  group_by(tissue, name) %>%
  nest() %>%
  mutate(
    stats = purrr::map(
      data,
      \(df) {
        cd_vals = df %>%
          filter(diet == "CD") %>%
          pull(med.adj.ic)
        id_vals = df %>%
          filter(diet == "ID") %>%
          pull(med.adj.ic)

        tibble(
          n_old_cd = length(cd_vals),
          n_old_id = length(id_vals),
          median_old_cd = median(cd_vals, na.rm = TRUE),
          median_old_id = median(id_vals, na.rm = TRUE),
          mean_old_cd = mean(cd_vals, na.rm = TRUE),
          mean_old_id = mean(id_vals, na.rm = TRUE),
          log2_fc_id_vs_cd = log2(median_old_id / median_old_cd)
        ) %>%
          bind_cols(
            if (length(cd_vals) >= 2 && length(id_vals) >= 2) {
              tryCatch(
                broom::tidy(t.test(med.adj.ic ~ diet, data = df)),
                error = \(e) tibble(
                  estimate = NA_real_,
                  estimate1 = NA_real_,
                  estimate2 = NA_real_,
                  statistic = NA_real_,
                  p.value = NA_real_,
                  parameter = NA_real_,
                  conf.low = NA_real_,
                  conf.high = NA_real_,
                  method = "Welch Two Sample t-test",
                  alternative = "two.sided"
                )
              )
            } else {
              tibble(
                estimate = NA_real_,
                estimate1 = NA_real_,
                estimate2 = NA_real_,
                statistic = NA_real_,
                p.value = NA_real_,
                parameter = NA_real_,
                conf.low = NA_real_,
                conf.high = NA_real_,
                method = "Welch Two Sample t-test",
                alternative = "two.sided"
              )
            }
          )
      }
    )
  ) %>%
  select(-data) %>%
  unnest(stats) %>%
  group_by(tissue) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(tissue, p.value, desc(abs(log2_fc_id_vs_cd)))

readr::write_csv(
  old_mouse_sample_abundance,
  here::here("data/processed/old mouse sample-level median adjusted ion counts.csv")
)

readr::write_csv(
  old_mouse_diet_ttest,
  here::here("data/processed/old mice inulin vs control t test median adjusted ion count.csv")
)


# ========== 2.0 - Old mouse CD vs ID heatmaps ==========
# -- keep only metabolites with nominal old CD vs ID differences, then z-scale within tissue
heatmap_p_cutoff = 0.05

old_mouse_heatmap_hits = old_mouse_diet_ttest %>%
  filter(!is.na(p.value), p.value < heatmap_p_cutoff) %>%
  select(tissue, name, p.value)

old_mouse_heatmap_data = old_mouse_sample_abundance %>%
  inner_join(old_mouse_heatmap_hits, by = c("tissue", "name")) %>%
  group_by(tissue, name) %>%
  mutate(
    z_abundance = if (sum(!is.na(med.adj.ic)) >= 2 & sd(med.adj.ic, na.rm = TRUE) > 0) {
      as.numeric(scale(med.adj.ic))
    } else {
      rep(0, dplyr::n())
    }
  ) %>%
  ungroup()

old_mouse_heatmaps = list()

for (current_tissue in unique(old_mouse_heatmap_data$tissue)) {

  current_column_data = old_mouse_heatmap_data %>%
    filter(tissue == current_tissue) %>%
    distinct(sample.clean, diet, cohort, animal) %>%
    arrange(diet, cohort, animal, sample.clean)

  current_matrix = old_mouse_heatmap_data %>%
    filter(tissue == current_tissue) %>%
    select(name, sample.clean, z_abundance) %>%
    tidyr::pivot_wider(
      names_from = sample.clean,
      values_from = z_abundance
    ) %>%
    tibble::column_to_rownames("name") %>%
    as.matrix()

  current_matrix = current_matrix[, current_column_data$sample.clean, drop = FALSE]
  current_matrix = current_matrix[rowSums(!is.na(current_matrix)) > 0, , drop = FALSE]
  current_row_order = old_mouse_heatmap_hits %>%
    filter(tissue == current_tissue, name %in% rownames(current_matrix)) %>%
    arrange(p.value, name) %>%
    pull(name)
  current_matrix = current_matrix[current_row_order, , drop = FALSE]
  current_row_labels = current_row_order %>%
    stringr::str_replace("^\\([^)]{1,25}\\)-", "") %>%
    stringr::str_replace("^\\([^)]{1,25}\\)", "") %>%
    stringr::str_replace_all("α", "alpha") %>%
    stringr::str_replace_all("β", "beta") %>%
    stringr::str_replace_all("γ", "gamma") %>%
    stringr::str_replace_all("δ", "delta") %>%
    stringr::str_squish()
  current_row_labels = ifelse(
    nchar(current_row_labels) > 50,
    stringr::str_trunc(current_row_labels, width = 50),
    current_row_labels
  )
  current_row_labels = make.unique(current_row_labels, sep = " | ")

  current_annotation = ComplexHeatmap::HeatmapAnnotation(
    diet = current_column_data$diet,
    cohort = current_column_data$cohort,
    col = list(
      diet = c(CD = "#1b9e77", ID = "#d95f02"),
      cohort = c("1" = "#7570b3", "2" = "#e7298a")
    ),
    annotation_legend_param = list(
      diet = list(direction = "horizontal", title_position = "leftcenter"),
      cohort = list(direction = "horizontal", title_position = "leftcenter")
    ),
    annotation_name_side = "left"
  )

  current_heatmap = ComplexHeatmap::Heatmap(
    current_matrix,
    name = "z",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_split = factor(current_column_data$diet, levels = c("CD", "ID")),
    top_annotation = current_annotation,
    na_col = "grey90",
    row_labels = current_row_labels,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    row_names_max_width = grid::unit(10, "cm"),
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8),
    row_title = paste(current_tissue, "metabolites"),
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "leftcenter"
    ),
    use_raster = TRUE,
    raster_quality = 2
  )

  old_mouse_heatmaps[[current_tissue]] = current_heatmap

  pdf(
    here::here("plots", paste0("old mice ", current_tissue, " CD vs ID z score heatmap.pdf")),
    width = 9,
    height = 12
  )

  ComplexHeatmap::draw(
    current_heatmap,
    heatmap_legend_side = "top",
    annotation_legend_side = "top",
    merge_legends = TRUE,
    legend_grouping = "original",
    padding = grid::unit(c(5, 5, 22, 5), "mm")
  )

  grid::grid.text(
    label = current_tissue,
    x = grid::unit(5, "mm"),
    y = grid::unit(1, "npc") - grid::unit(4, "mm"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  )

  grid::grid.text(
    label = paste0(
      "Old mice only | CD vs ID t-test hits | raw p < ",
      heatmap_p_cutoff,
      " | tiles: row z-scored median-adjusted ion counts"
    ),
    x = grid::unit(5, "mm"),
    y = grid::unit(1, "npc") - grid::unit(10, "mm"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9)
  )

  dev.off()
}
