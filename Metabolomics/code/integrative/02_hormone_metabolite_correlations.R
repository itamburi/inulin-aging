# ****** Inulin Aging Mouse Hormone-Metabolite Correlations
# ****** Summer 2026 - Ian T, Jessie K, Randa L
#
# NOTE: Secretin value for 2M33_SAC is excluded as an assay outlier.
# The value is 1747.92, while all other matched sacrifice Secretin values
# range from 11.69 to 60.33. This exclusion applies only to Secretin for
# sample 2M33; other hormones for 2M33 remain in the analysis.
#
# TRIAGED SET DEFINITION:
# The collaborator-facing triaged set starts from matched SAC observations
# only: 13 mice with hormone data and liver/muscle metabolomics. Pearson
# correlations are calculated for each tissue-metabolite-hormone combination.
# For triage, each scatter is filtered independently by removing observations
# where either metabolite abundance or hormone abundance has absolute modified
# z-score > 3.5. Pearson correlations are then recomputed after this outlier
# removal, requiring n >= 8. The final triaged set keeps correlations with
# nominal Pearson p < 0.05 and old-diet stratification, defined as absolute
# Cohen's d >= 1 for Old CD vs Old ID in either metabolite abundance or hormone
# abundance. This yielded 215 triaged correlations in the current dataset.


here::i_am("code/integrative/02_hormone_metabolite_correlations.R")


# ========== 0.0 - Environment ==========
# --
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)

dir.create(here::here("data/processed/integrative/hormone_metabolite_correlations"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("plots/integrative/hormone_metabolite_correlations"), recursive = TRUE, showWarnings = FALSE)

# Cairo PDF needs a writable font cache when this script runs in a sandboxed shell.
Sys.setenv(XDG_CACHE_HOME = tempdir())


# ========== 1.0 - Read processed metabolite and raw hormone data ==========
# --
metabolite_data = readr::read_csv(
  here::here("data/processed/tissue_metabolomics/cohort_comparisons/all group sample-level median adjusted ion counts.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  rename(
    sample = sample.clean,
    metabolite = name,
    metabolite_tissue = tissue,
    metabolite_abundance = med.adj.ic
  ) %>%
  select(
    sample, metabolite_tissue, metabolite,
    metabolite_abundance, metabolite_n_obs = n_obs,
    metabolite_group = group, metabolite_cohort = cohort,
    metabolite_age = age, metabolite_diet = diet,
    metabolite_animal = animal
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
    raw_sample_id = str_trim(raw_sample_id),
    sample = str_extract(raw_sample_id, "^[^_]+"),
    sampling_time = str_extract(raw_sample_id, "(?<=_).*"),
    hormone_cohort = str_extract(sample, "^\\d+"),
    hormone_animal = as.integer(str_extract(sample, "\\d+$"))
  ) %>%
  filter(sampling_time == "SAC") %>%
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
    sample, raw_sample_id, sampling_time,
    hormone, hormone_abundance,
    hormone_cohort, hormone_animal,
    plate, assay_group, well_id
  )


# ========== 2.0 - Join matching sacrifice mouse-level observations ==========
# --
matched_observations = metabolite_data %>%
  inner_join(hormone_data, by = "sample", relationship = "many-to-many") %>%
  filter(
    !is.na(metabolite_abundance),
    !is.na(hormone_abundance)
  )

readr::write_csv(
  matched_observations,
  here::here("data/processed/integrative/hormone_metabolite_correlations/matched sacrifice hormone metabolite observations.csv"),
  na = "N/A"
)


# ========== 3.0 - Calculate Pearson and Spearman correlations ==========
# --
calculate_correlations = function(data, correlation_method) {
  data %>%
    group_by(
      metabolite_tissue, metabolite,
      hormone
    ) %>%
    summarise(
      n = n(),
      n_unique_samples = n_distinct(sample),
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
      metabolite_tissue, metabolite,
      hormone,
      n, n_unique_samples,
      correlation, p_value, p_adj_bh,
      metabolite_min, metabolite_max,
      hormone_min, hormone_max
    ) %>%
    arrange(p_value, metabolite_tissue, metabolite, hormone)
}

pearson_correlations = calculate_correlations(matched_observations, "pearson")

spearman_correlations = calculate_correlations(matched_observations, "spearman")


# ========== 4.0 - Export correlation tables ==========
# --
readr::write_csv(
  pearson_correlations,
  here::here("data/processed/integrative/hormone_metabolite_correlations/pearson sacrifice hormone metabolite correlations.csv"),
  na = "N/A"
)

readr::write_csv(
  spearman_correlations,
  here::here("data/processed/integrative/hormone_metabolite_correlations/spearman sacrifice hormone metabolite correlations.csv"),
  na = "N/A"
)


# ========== 5.0 - Plot BH-significant Pearson correlations ==========
# --
pearson_significant_correlations = pearson_correlations %>%
  filter(!is.na(p_adj_bh), p_adj_bh < 0.05) %>%
  mutate(
    correlation_direction = if_else(correlation > 0, "positive", "negative"),
    plot_label = paste0(
      metabolite_tissue, " | ", metabolite, " | ", hormone,
      "\nr = ", round(correlation, 2),
      ", p_adj = ", signif(p_adj_bh, 3),
      ", n = ", n
    )
  )

plot_significant_correlations = function(correlation_direction) {
  plot_file = here::here(
    "plots/integrative/hormone_metabolite_correlations",
    paste0("pearson BH significant ", correlation_direction, " hormone metabolite correlations.pdf")
  )

  significant_correlations = pearson_significant_correlations %>%
    filter(correlation_direction == !!correlation_direction) %>%
    arrange(p_adj_bh, metabolite_tissue, metabolite, hormone) %>%
    mutate(plot_page = ceiling(row_number() / 12))

  if (nrow(significant_correlations) == 0) {
    if (file.exists(plot_file)) {
      file.remove(plot_file)
    }
    return(invisible(NULL))
  }

  grDevices::cairo_pdf(
    plot_file,
    width = 14,
    height = 10
  )

  for (current_page in sort(unique(significant_correlations$plot_page))) {
    page_correlations = significant_correlations %>%
      filter(plot_page == current_page)

    page_observations = matched_observations %>%
      inner_join(
        page_correlations %>%
          select(metabolite_tissue, metabolite, hormone, plot_label),
        by = c("metabolite_tissue", "metabolite", "hormone")
      ) %>%
      mutate(
        plot_label = factor(plot_label, levels = page_correlations$plot_label)
      )

    print(
      ggplot(
        page_observations,
        aes(x = metabolite_abundance, y = hormone_abundance)
      ) +
        geom_point(aes(color = metabolite_diet, shape = metabolite_age), size = 2.2, alpha = 0.85) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.35) +
        facet_wrap(~ plot_label, scales = "free", ncol = 4) +
        labs(
          title = paste(
            "Pearson BH-significant",
            correlation_direction,
            "hormone-metabolite correlations"
          ),
          subtitle = paste("Page", current_page, "of", max(significant_correlations$plot_page)),
          x = "Median-adjusted metabolite ion count",
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

plot_significant_correlations("positive")
plot_significant_correlations("negative")


# ========== 6.0 - Plot nominal p-value Pearson correlations ==========
# --
pearson_nominal_correlations = pearson_correlations %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  mutate(
    correlation_direction = if_else(correlation > 0, "positive", "negative"),
    plot_label = paste0(
      metabolite_tissue, " | ", metabolite, " | ", hormone,
      "\nr = ", round(correlation, 2),
      ", p = ", signif(p_value, 3),
      ", p_adj = ", signif(p_adj_bh, 3),
      ", n = ", n
    )
  )

plot_nominal_correlations = function(correlation_direction) {
  plot_file = here::here(
    "plots/integrative/hormone_metabolite_correlations",
    paste0("pearson nominal p 0.05 ", correlation_direction, " hormone metabolite correlations.pdf")
  )

  nominal_correlations = pearson_nominal_correlations %>%
    filter(correlation_direction == !!correlation_direction) %>%
    arrange(p_value, metabolite_tissue, metabolite, hormone) %>%
    mutate(plot_page = ceiling(row_number() / 12))

  if (nrow(nominal_correlations) == 0) {
    if (file.exists(plot_file)) {
      file.remove(plot_file)
    }
    return(invisible(NULL))
  }

  grDevices::cairo_pdf(
    plot_file,
    width = 14,
    height = 10
  )

  for (current_page in sort(unique(nominal_correlations$plot_page))) {
    page_correlations = nominal_correlations %>%
      filter(plot_page == current_page)

    page_observations = matched_observations %>%
      inner_join(
        page_correlations %>%
          select(metabolite_tissue, metabolite, hormone, plot_label),
        by = c("metabolite_tissue", "metabolite", "hormone")
      ) %>%
      mutate(
        plot_label = factor(plot_label, levels = page_correlations$plot_label)
      )

    print(
      ggplot(
        page_observations,
        aes(x = metabolite_abundance, y = hormone_abundance)
      ) +
        geom_point(aes(color = metabolite_diet, shape = metabolite_age), size = 2.2, alpha = 0.85) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.35) +
        facet_wrap(~ plot_label, scales = "free", ncol = 4) +
        labs(
          title = paste(
            "Pearson nominal p < 0.05",
            correlation_direction,
            "hormone-metabolite correlations"
          ),
          subtitle = paste("Page", current_page, "of", max(nominal_correlations$plot_page)),
          x = "Median-adjusted metabolite ion count",
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

plot_nominal_correlations("positive")
plot_nominal_correlations("negative")


# ========== 7.0 - Triage Pearson correlations with outlier and diet filters ==========
# --
# Outlier rule: remove observations with absolute modified z-score > 3.5
# within each tissue-metabolite-hormone pairing for either metabolite abundance
# or hormone abundance. Diet stratification is primarily evaluated among old
# CD vs old ID mice because young mice are CD only in the matched dataset.
calculate_modified_z = function(value) {
  value_median = median(value, na.rm = TRUE)
  value_mad = stats::mad(value, constant = 1, na.rm = TRUE)

  if (is.na(value_mad) || value_mad == 0) {
    rep(0, length(value))
  } else {
    0.6745 * (value - value_median) / value_mad
  }
}

calculate_cohens_d = function(value, group) {
  complete_observations = !is.na(value) & !is.na(group) & group %in% c("CD", "ID")
  value = value[complete_observations]
  group = group[complete_observations]

  cd_values = value[group == "CD"]
  id_values = value[group == "ID"]

  if (length(cd_values) < 2 || length(id_values) < 2 || dplyr::n_distinct(value) < 2) {
    return(NA_real_)
  }

  pooled_sd = sqrt(
    ((length(cd_values) - 1) * stats::var(cd_values) +
       (length(id_values) - 1) * stats::var(id_values)) /
      (length(cd_values) + length(id_values) - 2)
  )

  if (is.na(pooled_sd) || pooled_sd == 0) {
    NA_real_
  } else {
    (mean(id_values) - mean(cd_values)) / pooled_sd
  }
}

pearson_outlier_flagged_observations = matched_observations %>%
  group_by(metabolite_tissue, metabolite, hormone) %>%
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
  here::here("data/processed/integrative/hormone_metabolite_correlations/pearson triage outlier flagged observations.csv"),
  na = "N/A"
)

pearson_outlier_summary = pearson_outlier_flagged_observations %>%
  group_by(metabolite_tissue, metabolite, hormone) %>%
  summarise(
    n_before_outlier_filter = n(),
    n_outlier_removed = sum(pearson_triage_outlier),
    removed_samples = paste(sample[pearson_triage_outlier], collapse = ";"),
    removed_for_metabolite = paste(sample[metabolite_outlier], collapse = ";"),
    removed_for_hormone = paste(sample[hormone_outlier], collapse = ";"),
    .groups = "drop"
  )

pearson_outlier_filtered_observations = pearson_outlier_flagged_observations %>%
  filter(!pearson_triage_outlier)

calculate_outlier_filtered_pearson = function(data) {
  data %>%
    group_by(metabolite_tissue, metabolite, hormone) %>%
    summarise(
      n = n(),
      n_unique_samples = n_distinct(sample),
      n_cd = n_distinct(sample[metabolite_diet == "CD"]),
      n_id = n_distinct(sample[metabolite_diet == "ID"]),
      n_old_cd = n_distinct(sample[metabolite_age == "Old" & metabolite_diet == "CD"]),
      n_old_id = n_distinct(sample[metabolite_age == "Old" & metabolite_diet == "ID"]),
      metabolite_min = min(metabolite_abundance, na.rm = TRUE),
      metabolite_max = max(metabolite_abundance, na.rm = TRUE),
      hormone_min = min(hormone_abundance, na.rm = TRUE),
      hormone_max = max(hormone_abundance, na.rm = TRUE),
      old_diet_metabolite_cohens_d = calculate_cohens_d(
        metabolite_abundance[metabolite_age == "Old"],
        metabolite_diet[metabolite_age == "Old"]
      ),
      old_diet_hormone_cohens_d = calculate_cohens_d(
        hormone_abundance[metabolite_age == "Old"],
        metabolite_diet[metabolite_age == "Old"]
      ),
      all_diet_metabolite_cohens_d = calculate_cohens_d(metabolite_abundance, metabolite_diet),
      all_diet_hormone_cohens_d = calculate_cohens_d(hormone_abundance, metabolite_diet),
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
      p_adj_bh = p.adjust(p_value, method = "BH"),
      old_diet_stratifying = abs(old_diet_metabolite_cohens_d) >= 1 |
        abs(old_diet_hormone_cohens_d) >= 1,
      all_diet_stratifying = abs(all_diet_metabolite_cohens_d) >= 1 |
        abs(all_diet_hormone_cohens_d) >= 1
    ) %>%
    left_join(
      pearson_outlier_summary,
      by = c("metabolite_tissue", "metabolite", "hormone")
    ) %>%
    select(
      correlation_method,
      metabolite_tissue, metabolite, hormone,
      n, n_unique_samples,
      n_before_outlier_filter, n_outlier_removed,
      removed_samples, removed_for_metabolite, removed_for_hormone,
      n_cd, n_id, n_old_cd, n_old_id,
      correlation, p_value, p_adj_bh,
      old_diet_stratifying, all_diet_stratifying,
      old_diet_metabolite_cohens_d, old_diet_hormone_cohens_d,
      all_diet_metabolite_cohens_d, all_diet_hormone_cohens_d,
      metabolite_min, metabolite_max,
      hormone_min, hormone_max
    ) %>%
    arrange(p_value, metabolite_tissue, metabolite, hormone)
}

pearson_outlier_filtered_correlations = calculate_outlier_filtered_pearson(
  pearson_outlier_filtered_observations
)

pearson_outlier_filtered_diet_triage_correlations = pearson_outlier_filtered_correlations %>%
  filter(
    !is.na(p_value),
    p_value < 0.05,
    n_old_cd >= 3,
    n_old_id >= 3,
    old_diet_stratifying
  ) %>%
  arrange(p_value, desc(abs(old_diet_metabolite_cohens_d)), metabolite_tissue, metabolite, hormone)

readr::write_csv(
  pearson_outlier_filtered_correlations,
  here::here("data/processed/integrative/hormone_metabolite_correlations/pearson outlier filtered hormone metabolite correlations.csv"),
  na = "N/A"
)

readr::write_csv(
  pearson_outlier_filtered_diet_triage_correlations,
  here::here("data/processed/integrative/hormone_metabolite_correlations/pearson outlier filtered nominal p 0.05 old diet stratifying correlations.csv"),
  na = "N/A"
)


# ========== 8.0 - Plot outlier-filtered diet-triaged Pearson correlations ==========
# --
pearson_triage_plot_correlations = pearson_outlier_filtered_diet_triage_correlations %>%
  mutate(
    correlation_direction = if_else(correlation > 0, "positive", "negative"),
    plot_label = paste0(
      metabolite_tissue, " | ", metabolite, " | ", hormone,
      "\nr = ", round(correlation, 2),
      ", p = ", signif(p_value, 3),
      ", n = ", n,
      ", outliers removed = ", n_outlier_removed,
      "\nold diet d(metab) = ", round(old_diet_metabolite_cohens_d, 2),
      ", old diet d(hormone) = ", round(old_diet_hormone_cohens_d, 2)
    )
  )

plot_triaged_correlations = function(correlation_direction) {
  plot_file = here::here(
    "plots/integrative/hormone_metabolite_correlations",
    paste0(
      "pearson outlier filtered nominal p 0.05 old diet stratifying ",
      correlation_direction,
      " hormone metabolite correlations.pdf"
    )
  )

  triaged_correlations = pearson_triage_plot_correlations %>%
    filter(correlation_direction == !!correlation_direction) %>%
    arrange(p_value, metabolite_tissue, metabolite, hormone) %>%
    mutate(plot_page = ceiling(row_number() / 12))

  if (nrow(triaged_correlations) == 0) {
    if (file.exists(plot_file)) {
      file.remove(plot_file)
    }
    return(invisible(NULL))
  }

  grDevices::cairo_pdf(
    plot_file,
    width = 14,
    height = 10
  )

  for (current_page in sort(unique(triaged_correlations$plot_page))) {
    page_correlations = triaged_correlations %>%
      filter(plot_page == current_page)

    page_observations = pearson_outlier_filtered_observations %>%
      inner_join(
        page_correlations %>%
          select(metabolite_tissue, metabolite, hormone, plot_label),
        by = c("metabolite_tissue", "metabolite", "hormone")
      ) %>%
      mutate(
        plot_label = factor(plot_label, levels = page_correlations$plot_label)
      )

    print(
      ggplot(
        page_observations,
        aes(x = metabolite_abundance, y = hormone_abundance)
      ) +
        geom_point(aes(color = metabolite_diet, shape = metabolite_age), size = 2.2, alpha = 0.85) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.35) +
        facet_wrap(~ plot_label, scales = "free", ncol = 4) +
        labs(
          title = paste(
            "Outlier-filtered Pearson nominal p < 0.05 old diet-stratifying",
            correlation_direction,
            "hormone-metabolite correlations"
          ),
          subtitle = paste("Page", current_page, "of", max(triaged_correlations$plot_page)),
          x = "Median-adjusted metabolite ion count",
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

plot_triaged_correlations("positive")
plot_triaged_correlations("negative")


# ========== 9.0 - Export simplified collaborator-facing triage audit ==========
# --
scale_for_clustering = function(value) {
  value_sd = stats::sd(value, na.rm = TRUE)

  if (is.na(value_sd) || value_sd == 0) {
    rep(0, length(value))
  } else {
    as.numeric(scale(value))
  }
}

calculate_centroid_distance = function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

triage_group_centroids = pearson_outlier_filtered_observations %>%
  semi_join(
    pearson_outlier_filtered_diet_triage_correlations %>%
      select(metabolite_tissue, metabolite, hormone),
    by = c("metabolite_tissue", "metabolite", "hormone")
  ) %>%
  mutate(
    plot_group = case_when(
      metabolite_age == "Young" & metabolite_diet == "CD" ~ "Young CD",
      metabolite_age == "Old" & metabolite_diet == "CD" ~ "Old CD",
      metabolite_age == "Old" & metabolite_diet == "ID" ~ "Old ID",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(plot_group)) %>%
  group_by(metabolite_tissue, metabolite, hormone) %>%
  mutate(
    metabolite_clustering_z = scale_for_clustering(log10(metabolite_abundance + 1)),
    hormone_clustering_z = scale_for_clustering(hormone_abundance)
  ) %>%
  ungroup() %>%
  group_by(metabolite_tissue, metabolite, hormone, plot_group) %>%
  summarise(
    n_group = n_distinct(sample),
    centroid_metabolite = mean(metabolite_clustering_z, na.rm = TRUE),
    centroid_hormone = mean(hormone_clustering_z, na.rm = TRUE),
    mean_metabolite_abundance = mean(metabolite_abundance, na.rm = TRUE),
    mean_hormone_abundance = mean(hormone_abundance, na.rm = TRUE),
    .groups = "drop"
  )

triage_clustering_audit = triage_group_centroids %>%
  select(
    metabolite_tissue, metabolite, hormone, plot_group,
    n_group, centroid_metabolite, centroid_hormone,
    mean_metabolite_abundance, mean_hormone_abundance
  ) %>%
  pivot_wider(
    names_from = plot_group,
    values_from = c(
      n_group, centroid_metabolite, centroid_hormone,
      mean_metabolite_abundance, mean_hormone_abundance
    ),
    names_sep = "__"
  ) %>%
  mutate(
    distance_young_cd_to_old_id = calculate_centroid_distance(
      `centroid_metabolite__Young CD`, `centroid_hormone__Young CD`,
      `centroid_metabolite__Old ID`, `centroid_hormone__Old ID`
    ),
    distance_young_cd_to_old_cd = calculate_centroid_distance(
      `centroid_metabolite__Young CD`, `centroid_hormone__Young CD`,
      `centroid_metabolite__Old CD`, `centroid_hormone__Old CD`
    ),
    distance_old_id_to_old_cd = calculate_centroid_distance(
      `centroid_metabolite__Old ID`, `centroid_hormone__Old ID`,
      `centroid_metabolite__Old CD`, `centroid_hormone__Old CD`
    ),
    closest_group_pattern = case_when(
      distance_young_cd_to_old_id <= distance_young_cd_to_old_cd &
        distance_young_cd_to_old_id <= distance_old_id_to_old_cd ~ "Old ID resembles Young CD",
      distance_young_cd_to_old_cd <= distance_young_cd_to_old_id &
        distance_young_cd_to_old_cd <= distance_old_id_to_old_cd ~ "CD groups cluster; Old ID separate",
      distance_old_id_to_old_cd <= distance_young_cd_to_old_id &
        distance_old_id_to_old_cd <= distance_young_cd_to_old_cd ~ "Old groups cluster; Young CD separate",
      TRUE ~ "Mixed/no clear centroid clustering"
    ),
    clustering_strength = case_when(
      pmin(
        distance_young_cd_to_old_id,
        distance_young_cd_to_old_cd,
        distance_old_id_to_old_cd,
        na.rm = TRUE
      ) == 0 ~ NA_real_,
      TRUE ~ pmax(
        distance_young_cd_to_old_id,
        distance_young_cd_to_old_cd,
        distance_old_id_to_old_cd,
        na.rm = TRUE
      ) /
        pmin(
          distance_young_cd_to_old_id,
          distance_young_cd_to_old_cd,
          distance_old_id_to_old_cd,
          na.rm = TRUE
        )
    ),
    clustering_category = case_when(
      is.na(clustering_strength) | clustering_strength < 1.5 ~ "Mixed/no strong centroid clustering",
      TRUE ~ closest_group_pattern
    )
  )

collaborator_triage_audit = pearson_outlier_filtered_diet_triage_correlations %>%
  left_join(
    triage_clustering_audit,
    by = c("metabolite_tissue", "metabolite", "hormone")
  ) %>%
  mutate(
    correlation_direction = if_else(correlation > 0, "positive", "negative")
  ) %>%
  select(
    metabolite_tissue, metabolite, hormone,
    correlation_direction,
    pearson_r = correlation,
    pearson_p_value = p_value,
    pearson_p_adj_bh = p_adj_bh,
    n,
    n_outlier_removed,
    removed_samples,
    n_old_cd, n_old_id,
    clustering_category,
    old_diet_metabolite_cohens_d, old_diet_hormone_cohens_d,
    `mean_metabolite_abundance__Old CD`,
    `mean_metabolite_abundance__Old ID`,
    `mean_metabolite_abundance__Young CD`,
    `mean_hormone_abundance__Old CD`,
    `mean_hormone_abundance__Old ID`,
    `mean_hormone_abundance__Young CD`
  ) %>%
  arrange(
    clustering_category,
    pearson_p_value,
    metabolite_tissue,
    metabolite,
    hormone
  )

readr::write_csv(
  collaborator_triage_audit,
  here::here("data/processed/integrative/hormone_metabolite_correlations/collaborator pearson triage audit simplified.csv"),
  na = "N/A"
)


# ========== 10.0 - Plot collaborator-facing triage audit correlations ==========
# --
collaborator_triage_plot_correlations = collaborator_triage_audit %>%
  mutate(
    plot_label = paste0(
      metabolite_tissue, " | ", metabolite, " | ", hormone,
      "\n", clustering_category,
      "\nr = ", round(pearson_r, 2),
      ", p = ", signif(pearson_p_value, 3),
      ", n = ", n,
      ", outliers removed = ", n_outlier_removed
    ),
    plot_page = ceiling(row_number() / 12)
  )

grDevices::cairo_pdf(
  here::here("plots/integrative/hormone_metabolite_correlations/collaborator pearson triage audit scatters.pdf"),
  width = 14,
  height = 10
)

for (current_page in sort(unique(collaborator_triage_plot_correlations$plot_page))) {
  page_correlations = collaborator_triage_plot_correlations %>%
    filter(plot_page == current_page)

  page_observations = pearson_outlier_filtered_observations %>%
    inner_join(
      page_correlations %>%
        select(metabolite_tissue, metabolite, hormone, plot_label),
      by = c("metabolite_tissue", "metabolite", "hormone")
    ) %>%
    mutate(
      plot_label = factor(plot_label, levels = page_correlations$plot_label)
    )

  print(
    ggplot(
      page_observations,
      aes(x = metabolite_abundance, y = hormone_abundance)
    ) +
      geom_point(aes(color = metabolite_diet, shape = metabolite_age), size = 2.2, alpha = 0.85) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.35) +
      facet_wrap(~ plot_label, scales = "free", ncol = 4) +
      labs(
        title = "Collaborator triage audit hormone-metabolite correlations",
        subtitle = paste(
          "Outlier-filtered Pearson nominal p < 0.05 old diet-stratifying correlations; page",
          current_page,
          "of",
          max(collaborator_triage_plot_correlations$plot_page)
        ),
        x = "Median-adjusted metabolite ion count",
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


# ========== 11.0 - Plot top outlier-filtered Leptin Pearson correlations ==========
# --
top_leptin_correlations = pearson_outlier_filtered_correlations %>%
  filter(
    hormone == "Leptin",
    !is.na(p_value),
    p_value < 0.05,
    abs(correlation) >= 0.8
  ) %>%
  arrange(p_value, metabolite_tissue, metabolite) %>%
  mutate(
    correlation_direction = if_else(correlation > 0, "positive", "negative"),
    plot_label = paste0(
      metabolite_tissue, " | ", metabolite,
      "\nr = ", round(correlation, 2),
      ", p = ", signif(p_value, 3),
      ", n = ", n,
      ", outliers removed = ", n_outlier_removed
    ),
    plot_page = ceiling(row_number() / 12)
  )

grDevices::cairo_pdf(
  here::here("plots/integrative/hormone_metabolite_correlations/top leptin outlier filtered pearson correlations.pdf"),
  width = 14,
  height = 10
)

for (current_page in sort(unique(top_leptin_correlations$plot_page))) {
  page_correlations = top_leptin_correlations %>%
    filter(plot_page == current_page)

  page_observations = pearson_outlier_filtered_observations %>%
    inner_join(
      page_correlations %>%
        select(metabolite_tissue, metabolite, hormone, plot_label),
      by = c("metabolite_tissue", "metabolite", "hormone")
    ) %>%
    mutate(
      plot_label = factor(plot_label, levels = page_correlations$plot_label)
    )

  print(
    ggplot(
      page_observations,
      aes(x = metabolite_abundance, y = hormone_abundance)
    ) +
      geom_point(aes(color = metabolite_diet, shape = metabolite_age), size = 2.2, alpha = 0.85) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.35) +
      facet_wrap(~ plot_label, scales = "free", ncol = 4) +
      labs(
        title = "Top outlier-filtered Leptin-metabolite Pearson correlations",
        subtitle = paste(
          "Leptin only; nominal p < 0.05 and |r| >= 0.8; page",
          current_page,
          "of",
          max(top_leptin_correlations$plot_page)
        ),
        x = "Median-adjusted metabolite ion count",
        y = "Leptin abundance",
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
