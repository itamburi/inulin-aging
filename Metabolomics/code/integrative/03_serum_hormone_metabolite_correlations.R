# ****** Inulin Aging Mouse Serum Hormone-Metabolite Correlations
# ****** Summer 2026 - Ian T, Jessie K, Randa L
#
# Pairings:
# - PV and systemic sacrifice serum metabolites are matched to SAC hormone levels.
# - TP4 serum metabolites are matched to GTT hormone levels.
# - Metabolite abundance uses median-adjusted ion count when available, with raw
#   ion count as fallback for samples lacking normalization, such as TP4.


here::i_am("code/integrative/03_serum_hormone_metabolite_correlations.R")


# ========== 0.0 - Environment ==========
# --
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)

output_dir = here::here("data/processed/integrative/serum_hormone_metabolite_correlations")
plot_dir = here::here("plots/integrative/serum_hormone_metabolite_correlations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

top_n_scatter_hits = 96

# Cairo PDF needs a writable font cache when this script runs in a sandboxed shell.
Sys.setenv(XDG_CACHE_HOME = tempdir())


# ========== 1.0 - Read processed serum metabolite and raw hormone data ==========
# --

official_sample_mapping = readr::read_csv(
  here::here("data/libraries/hormone to metabolomics sample mapping MLL.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  filter(
    metabolomics_tissue %in% c("Serum PV", "Serum systemic SAC", "Serum TP4")
  ) %>%
  mutate(
    sample = stringr::str_extract(hormone_sample_id, "^[^_]+"),
    metabolite_timepoint = dplyr::case_when(
      metabolomics_tissue == "Serum PV" ~ "pv",
      metabolomics_tissue == "Serum systemic SAC" ~ "sys_sac",
      metabolomics_tissue == "Serum TP4" ~ "tp4",
      TRUE ~ NA_character_
    ),
    hormone_sampling_time = dplyr::case_when(
      hormone_tissue == "Blood SAC" ~ "SAC",
      hormone_tissue == "Blood GTT" ~ "GTT",
      TRUE ~ NA_character_
    ),
    metabolite_diet = factor(diet, levels = c("CD", "ID")),
    metabolite_age = factor(age, levels = c("Young", "Old")),
    metabolite_sex = sex,
    metabolite_cohort = cohort,
    metabolite_animal = as.integer(stringr::str_extract(sample, "\\d+$"))
  ) %>%
  select(
    sample,
    hormone_sample_id,
    metabolite_sample_id = metabolomics_sample_id,
    metabolite_timepoint,
    hormone_sampling_time,
    metabolite_cohort,
    metabolite_animal,
    metabolite_diet,
    metabolite_age,
    metabolite_sex
  )

serum_metabolite_data = readr::read_csv(
  here::here("data/processed/serum_metabolomics/ion_counts/tidy and normalized serum ion count data.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  mutate(
    metabolite_abundance = dplyr::coalesce(med.adj.ic, ion.count),
    metabolite_abundance_source = dplyr::case_when(
      !is.na(med.adj.ic) ~ "med.adj.ic",
      is.na(med.adj.ic) & !is.na(ion.count) ~ "ion.count",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    metabolomics_mouse_id = sample,
    metabolite_sample_id = sample.id,
    metabolite_timepoint = timepoint,
    metabolite_source = source,
    metabolite_ion_mode = ion.mode,
    metabolite = name,
    metabolite_formula = formula,
    metabolite_mz = mz,
    metabolite_rt = rt,
    metabolite_abundance,
    metabolite_abundance_source,
    metabolite_raw_ion_count = ion.count,
    metabolite_med_adj_ic = med.adj.ic,
    parsed_metabolite_cohort = cohort
  ) %>%
  filter(
    metabolite_timepoint %in% c("pv", "sys_sac", "tp4"),
    !is.na(metabolite_abundance)
  )

hormone_data = readr::read_csv(
  here::here("data/raw/hormones/luminex_hormone_inulin_aging_raw.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "-")
) %>%
  rename(
    plate = Plate,
    assay_group = Group,
    well_id = `Well ID`,
    raw_sample_id = `Sample ID`
  ) %>%
  mutate(
    raw_sample_id = stringr::str_trim(raw_sample_id),
    sample = stringr::str_extract(raw_sample_id, "^[^_]+"),
    hormone_sampling_time = stringr::str_extract(raw_sample_id, "(?<=_).*"),
    hormone_cohort = stringr::str_extract(sample, "^\\d+"),
    hormone_animal = as.integer(stringr::str_extract(sample, "\\d+$"))
  ) %>%
  filter(hormone_sampling_time %in% c("SAC", "GTT")) %>%
  pivot_longer(
    cols = c(
      `C-Peptide 2`, Ghrelin, GIP, `GLP-1`, Glucagon, Insulin,
      Leptin, Secretin, PP, PYY, Resistin, Amylin
    ),
    names_to = "hormone",
    values_to = "hormone_abundance"
  ) %>%
  filter(
    !(sample == "2M33" & raw_sample_id == "2M33_SAC" & hormone == "Secretin")
  ) %>%
  select(
    hormone_mouse_id = sample, raw_sample_id, hormone_sampling_time,
    hormone, hormone_abundance,
    hormone_cohort, hormone_animal,
    plate, assay_group, well_id
  )


# ========== 2.0 - Join matching serum metabolite and hormone observations ==========
# --

matched_observations = official_sample_mapping %>%
  inner_join(
    serum_metabolite_data,
    by = c("metabolite_sample_id", "metabolite_timepoint"),
    relationship = "many-to-many"
  ) %>%
  inner_join(
    hormone_data,
    by = c(
      "hormone_sample_id" = "raw_sample_id",
      "hormone_sampling_time"
    ),
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(hormone_abundance))

readr::write_csv(
  matched_observations,
  file.path(output_dir, "matched serum hormone metabolite observations.csv"),
  na = "N/A"
)

matched_sample_summary = matched_observations %>%
  distinct(sample, metabolite_timepoint, hormone_sampling_time) %>%
  count(metabolite_timepoint, hormone_sampling_time, name = "n_matched_samples")

readr::write_csv(
  matched_sample_summary,
  file.path(output_dir, "matched serum hormone sample summary.csv"),
  na = "N/A"
)


# ========== 3.0 - Calculate Pearson and Spearman correlations ==========
# --

calculate_correlations = function(data, correlation_method) {
  data %>%
    group_by(
      metabolite_timepoint, hormone_sampling_time,
      metabolite_source, metabolite_ion_mode,
      metabolite, metabolite_formula, metabolite_mz, metabolite_rt,
      hormone
    ) %>%
    summarise(
      n = n(),
      n_unique_samples = n_distinct(sample),
      n_med_adj_ic = sum(metabolite_abundance_source == "med.adj.ic", na.rm = TRUE),
      n_raw_ion_count_fallback = sum(metabolite_abundance_source == "ion.count", na.rm = TRUE),
      metabolite_min = min(metabolite_abundance, na.rm = TRUE),
      metabolite_max = max(metabolite_abundance, na.rm = TRUE),
      hormone_min = min(hormone_abundance, na.rm = TRUE),
      hormone_max = max(hormone_abundance, na.rm = TRUE),
      test = list(
        if (
          n() >= 3 &&
            dplyr::n_distinct(metabolite_abundance) > 1 &&
            dplyr::n_distinct(hormone_abundance) > 1
        ) {
          cor.test(
            metabolite_abundance,
            hormone_abundance,
            method = correlation_method,
            exact = FALSE
          )
        } else {
          NULL
        }
      ),
      .groups = "drop"
    ) %>%
    mutate(
      correlation_method = correlation_method,
      correlation = purrr::map_dbl(test, ~ if (is.null(.x)) NA_real_ else unname(.x$estimate)),
      p_value = purrr::map_dbl(test, ~ if (is.null(.x)) NA_real_ else .x$p.value),
      p_adj_bh = p.adjust(p_value, method = "BH")
    ) %>%
    select(
      correlation_method,
      metabolite_timepoint, hormone_sampling_time,
      metabolite_source, metabolite_ion_mode,
      metabolite, metabolite_formula, metabolite_mz, metabolite_rt,
      hormone,
      n, n_unique_samples, n_med_adj_ic, n_raw_ion_count_fallback,
      correlation, p_value, p_adj_bh,
      metabolite_min, metabolite_max,
      hormone_min, hormone_max
    ) %>%
    arrange(p_value, metabolite_timepoint, metabolite_source, metabolite, hormone)
}

pearson_correlations = calculate_correlations(matched_observations, "pearson")

spearman_correlations = calculate_correlations(matched_observations, "spearman")

readr::write_csv(
  pearson_correlations,
  file.path(output_dir, "pearson serum hormone metabolite correlations.csv"),
  na = "N/A"
)

readr::write_csv(
  spearman_correlations,
  file.path(output_dir, "spearman serum hormone metabolite correlations.csv"),
  na = "N/A"
)


# ========== 4.0 - Export nominal Pearson hits for quick review ==========
# --

pearson_nominal_correlations = pearson_correlations %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  arrange(p_value, metabolite_timepoint, hormone_sampling_time, metabolite, hormone)

readr::write_csv(
  pearson_nominal_correlations,
  file.path(output_dir, "pearson nominal p 0.05 serum hormone metabolite correlations.csv"),
  na = "N/A"
)


# ========== 5.0 - Recalculate Pearson correlations after outlier triage ==========
# --
# Match the tissue-hormone triage: flag observations with absolute modified
# z-score > 3.5 within each serum metabolite-hormone scatter, then recompute
# Pearson after removing flagged points.
calculate_modified_z = function(value) {
  value_median = median(value, na.rm = TRUE)
  value_mad = stats::mad(value, constant = 1, na.rm = TRUE)

  if (is.na(value_mad) || value_mad == 0) {
    rep(0, length(value))
  } else {
    0.6745 * (value - value_median) / value_mad
  }
}

pearson_outlier_flagged_observations = matched_observations %>%
  group_by(
    metabolite_timepoint, hormone_sampling_time,
    metabolite_source, metabolite_ion_mode,
    metabolite, metabolite_formula, metabolite_mz, metabolite_rt,
    hormone
  ) %>%
  mutate(
    metabolite_modified_z = calculate_modified_z(metabolite_abundance),
    hormone_modified_z = calculate_modified_z(hormone_abundance),
    metabolite_outlier = abs(metabolite_modified_z) > 3.5,
    hormone_outlier = abs(hormone_modified_z) > 3.5,
    pearson_triage_outlier = metabolite_outlier | hormone_outlier
  ) %>%
  ungroup()

readr::write_csv(
  pearson_outlier_flagged_observations %>%
    filter(pearson_triage_outlier) %>%
    arrange(desc(abs(metabolite_modified_z)), desc(abs(hormone_modified_z))),
  file.path(output_dir, "pearson triage outlier flagged serum hormone metabolite observations.csv"),
  na = "N/A"
)

pearson_outlier_summary = pearson_outlier_flagged_observations %>%
  group_by(
    metabolite_timepoint, hormone_sampling_time,
    metabolite_source, metabolite_ion_mode,
    metabolite, metabolite_formula, metabolite_mz, metabolite_rt,
    hormone
  ) %>%
  summarise(
    n_before_outlier_filter = n(),
    n_outlier_removed = sum(pearson_triage_outlier),
    removed_samples = paste(sample[pearson_triage_outlier], collapse = ";"),
    removed_for_metabolite = paste(sample[metabolite_outlier], collapse = ";"),
    removed_for_hormone = paste(sample[hormone_outlier], collapse = ";"),
    .groups = "drop"
  )

readr::write_csv(
  pearson_outlier_summary,
  file.path(output_dir, "pearson outlier summary serum hormone metabolite correlations.csv"),
  na = "N/A"
)

pearson_outlier_filtered_observations = pearson_outlier_flagged_observations %>%
  filter(!pearson_triage_outlier)

calculate_outlier_filtered_pearson = function(data) {
  data %>%
    group_by(
      metabolite_timepoint, hormone_sampling_time,
      metabolite_source, metabolite_ion_mode,
      metabolite, metabolite_formula, metabolite_mz, metabolite_rt,
      hormone
    ) %>%
    summarise(
      n = n(),
      n_unique_samples = n_distinct(sample),
      n_med_adj_ic = sum(metabolite_abundance_source == "med.adj.ic", na.rm = TRUE),
      n_raw_ion_count_fallback = sum(metabolite_abundance_source == "ion.count", na.rm = TRUE),
      n_cd = n_distinct(sample[metabolite_diet == "CD"]),
      n_id = n_distinct(sample[metabolite_diet == "ID"]),
      n_young = n_distinct(sample[metabolite_age == "Young"]),
      n_old = n_distinct(sample[metabolite_age == "Old"]),
      metabolite_min = min(metabolite_abundance, na.rm = TRUE),
      metabolite_max = max(metabolite_abundance, na.rm = TRUE),
      hormone_min = min(hormone_abundance, na.rm = TRUE),
      hormone_max = max(hormone_abundance, na.rm = TRUE),
      test = list(
        if (
          n() >= 8 &&
            dplyr::n_distinct(metabolite_abundance) > 1 &&
            dplyr::n_distinct(hormone_abundance) > 1
        ) {
          cor.test(
            metabolite_abundance,
            hormone_abundance,
            method = "pearson",
            exact = FALSE
          )
        } else {
          NULL
        }
      ),
      .groups = "drop"
    ) %>%
    mutate(
      correlation_method = "pearson",
      correlation = purrr::map_dbl(test, ~ if (is.null(.x)) NA_real_ else unname(.x$estimate)),
      p_value = purrr::map_dbl(test, ~ if (is.null(.x)) NA_real_ else .x$p.value),
      p_adj_bh = p.adjust(p_value, method = "BH")
    ) %>%
    left_join(
      pearson_outlier_summary,
      by = c(
        "metabolite_timepoint", "hormone_sampling_time",
        "metabolite_source", "metabolite_ion_mode",
        "metabolite", "metabolite_formula", "metabolite_mz", "metabolite_rt",
        "hormone"
      )
    ) %>%
    select(
      correlation_method,
      metabolite_timepoint, hormone_sampling_time,
      metabolite_source, metabolite_ion_mode,
      metabolite, metabolite_formula, metabolite_mz, metabolite_rt,
      hormone,
      n, n_unique_samples,
      n_before_outlier_filter, n_outlier_removed,
      removed_samples, removed_for_metabolite, removed_for_hormone,
      n_med_adj_ic, n_raw_ion_count_fallback,
      n_cd, n_id, n_young, n_old,
      correlation, p_value, p_adj_bh,
      metabolite_min, metabolite_max,
      hormone_min, hormone_max
    ) %>%
    arrange(p_value, metabolite_timepoint, metabolite_source, metabolite, hormone)
}

pearson_outlier_filtered_correlations = calculate_outlier_filtered_pearson(
  pearson_outlier_filtered_observations
)

pearson_outlier_filtered_nominal_correlations = pearson_outlier_filtered_correlations %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  arrange(p_value, metabolite_timepoint, hormone_sampling_time, metabolite, hormone)

readr::write_csv(
  pearson_outlier_filtered_correlations,
  file.path(output_dir, "pearson outlier filtered serum hormone metabolite correlations.csv"),
  na = "N/A"
)

readr::write_csv(
  pearson_outlier_filtered_nominal_correlations,
  file.path(output_dir, "pearson outlier filtered nominal p 0.05 serum hormone metabolite correlations.csv"),
  na = "N/A"
)


# ========== 6.0 - Plot top outlier-filtered Pearson hits by serum-hormone pairing ==========
# --

make_plot_metabolite_label = function(metabolite) {
  metabolite %>%
    stringr::str_replace("^\\([^)]{1,25}\\)-", "") %>%
    stringr::str_replace("^\\([^)]{1,25}\\)", "") %>%
    stringr::str_replace_all("α", "alpha") %>%
    stringr::str_replace_all("β", "beta") %>%
    stringr::str_replace_all("γ", "gamma") %>%
    stringr::str_replace_all("δ", "delta") %>%
    stringr::str_squish() %>%
    stringr::str_trunc(width = 55)
}

plot_top_pearson_hits = function(current_metabolite_timepoint, current_hormone_sampling_time, output_file) {
  current_correlations = pearson_outlier_filtered_correlations %>%
    filter(
      metabolite_timepoint == current_metabolite_timepoint,
      hormone_sampling_time == current_hormone_sampling_time,
      !is.na(p_value)
    ) %>%
    arrange(p_value, desc(abs(correlation)), metabolite, hormone) %>%
    slice_head(n = top_n_scatter_hits) %>%
    mutate(
      plot_label = paste0(
        make_plot_metabolite_label(metabolite), " | ", hormone,
        "\nr = ", round(correlation, 2),
        ", p = ", signif(p_value, 3),
        ", p_adj = ", signif(p_adj_bh, 3),
        ", n = ", n_unique_samples,
        ", outliers removed = ", n_outlier_removed,
        "\n", metabolite_source, " ", metabolite_ion_mode,
        " | med.adj n = ", n_med_adj_ic,
        ", raw n = ", n_raw_ion_count_fallback
      ),
      plot_page = ceiling(row_number() / 12)
    )

  if (nrow(current_correlations) == 0) {
    if (file.exists(output_file)) {
      file.remove(output_file)
    }
    return(invisible(NULL))
  }

  grDevices::cairo_pdf(
    output_file,
    width = 14,
    height = 10
  )

  for (current_page in sort(unique(current_correlations$plot_page))) {
    page_correlations = current_correlations %>%
      filter(plot_page == current_page)

    page_observations = pearson_outlier_filtered_observations %>%
      inner_join(
        page_correlations %>%
          select(
            metabolite_timepoint, hormone_sampling_time,
            metabolite_source, metabolite_ion_mode,
            metabolite, metabolite_formula, metabolite_mz, metabolite_rt,
            hormone, plot_label
          ),
        by = c(
          "metabolite_timepoint", "hormone_sampling_time",
          "metabolite_source", "metabolite_ion_mode",
          "metabolite", "metabolite_formula", "metabolite_mz", "metabolite_rt",
          "hormone"
        )
      ) %>%
      mutate(
        plot_label = factor(plot_label, levels = page_correlations$plot_label)
      )

    print(
      ggplot(
        page_observations,
        aes(x = metabolite_abundance, y = hormone_abundance)
      ) +
        geom_point(aes(color = metabolite_diet, shape = metabolite_age), size = 2.1, alpha = 0.85) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.35) +
        facet_wrap(~ plot_label, scales = "free", ncol = 4) +
        scale_color_manual(
          values = c(CD = "#0057b8", ID = "#00843d"),
          na.value = "grey50"
        ) +
        scale_shape_manual(
          values = c("Young" = 17, "Old" = 16),
          na.value = 4
        ) +
        labs(
          title = paste0(
            "Top ", top_n_scatter_hits, " outlier-filtered Pearson serum hormone-metabolite hits: ",
            current_metabolite_timepoint, " metabolites vs ", current_hormone_sampling_time, " hormones"
          ),
          subtitle = paste("Page", current_page, "of", max(current_correlations$plot_page)),
          x = "Serum metabolite abundance",
          y = "Hormone abundance",
          color = "Diet",
          shape = "Age"
        ) +
        theme_bw(base_size = 9) +
        theme(
          strip.text = element_text(size = 6.5),
          axis.text = element_text(size = 6),
          legend.position = "bottom",
          plot.title = element_text(face = "bold")
        )
    )
  }

  dev.off()
}

plot_top_pearson_hits(
  "pv",
  "SAC",
  file.path(plot_dir, "pearson outlier filtered top 96 hits pv vs sac hormones.pdf")
)

plot_top_pearson_hits(
  "sys_sac",
  "SAC",
  file.path(plot_dir, "pearson outlier filtered top 96 hits sys sac vs sac hormones.pdf")
)

plot_top_pearson_hits(
  "tp4",
  "GTT",
  file.path(plot_dir, "pearson outlier filtered top 96 hits tp4 vs gtt hormones.pdf")
)
