# ****** Inulin Aging Mouse Blood Metabolite-Hormone Correlations
# ****** Summer 2026 - Ian T, Jessie K, Randa L
#
# Pairings:
# - PV serum metabolites are matched to SAC hormone levels.
# - Systemic sacrifice serum metabolites are matched to SAC hormone levels.
# - TP4 serum metabolites are matched to GTT hormone levels.
# - TP4/GTT and systemic sacrifice/SAC matched observations are combined.
#
# Correlations are run among old animals only, separately for CD, ID, and the
# combined CD/ID group. No filtering or plotting is applied in this script.


here::i_am("code/integrative/03_serum_hormone_metabolite_correlations.R")


# ========== 0.0 - Environment ==========
# --
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

output_dir = here::here("data/processed/integrative/serum_hormone_metabolite_correlations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hormone_columns = c(
  "C-Peptide 2", "Ghrelin", "GIP", "GLP-1", "Glucagon", "Insulin",
  "Leptin", "Secretin", "PP", "PYY", "Resistin", "Amylin"
)

correlation_methods = c("pearson", "spearman", "bicor")
minimum_n_for_correlation = 3


# ========== 1.0 - Read serum sample inventory, serum metabolites, and hormones ==========
# --

serum_sample_inventory = readr::read_csv(
  here::here("data/processed/sample_inventory/serum_compound_discoverer_sample_inventory.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  filter(
    dataset == "serum",
    sample_type %in% c("pv", "sys_sac", "tp4")
  ) %>%
  mutate(
    metabolite_sample_id = paste0("Serum_", sample, "_", sample_type),
    metabolite_timepoint = sample_type,
    hormone_sampling_time = case_when(
      metabolite_timepoint %in% c("pv", "sys_sac") ~ "SAC",
      metabolite_timepoint == "tp4" ~ "GTT",
      TRUE ~ NA_character_
    ),
    hormone_sample_id = paste0(sample, "_", hormone_sampling_time),
    metabolomics_tissue = case_when(
      metabolite_timepoint == "pv" ~ "Serum PV",
      metabolite_timepoint == "sys_sac" ~ "Serum systemic SAC",
      metabolite_timepoint == "tp4" ~ "Serum TP4",
      TRUE ~ NA_character_
    ),
    hormone_tissue = paste("Blood", hormone_sampling_time),
    metabolite_cohort = as.integer(cohort),
    metabolite_animal = as.integer(stringr::str_extract(sample, "\\d+$")),
    metabolite_age = factor(age, levels = c("Young", "Old")),
    metabolite_diet = factor(diet, levels = c("CD", "ID")),
    metabolite_gender = gender,
    metabolite_sex = dplyr::recode(gender, Female = "F", Male = "M", .default = gender)
  ) %>%
  group_by(
    sample, sample_type, hormone_sample_id, metabolite_sample_id,
    metabolite_timepoint, hormone_sampling_time,
    metabolomics_tissue, hormone_tissue
  ) %>%
  summarise(
    metabolite_cohort = first(metabolite_cohort),
    metabolite_animal = first(metabolite_animal),
    metabolite_diet = first(metabolite_diet),
    metabolite_age = first(metabolite_age),
    metabolite_gender = first(metabolite_gender),
    metabolite_sex = first(metabolite_sex),
    died = paste(sort(unique(died[!is.na(died)])), collapse = ";"),
    timepoint_months = first(timepoint_months),
    blood_type = first(blood_type),
    collection_protocol = first(collection_protocol),
    .groups = "drop"
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

raw_hormone_data = readr::read_csv(
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
  )

serum_metabolite_samples = serum_metabolite_data %>%
  distinct(metabolite_sample_id, metabolite_timepoint)

hormone_samples = raw_hormone_data %>%
  filter(hormone_sampling_time %in% c("SAC", "GTT")) %>%
  distinct(
    hormone_sample_id = raw_sample_id,
    hormone_sampling_time,
    hormone_sample_mouse_id = sample,
    hormone_sample_cohort = hormone_cohort,
    hormone_sample_animal = hormone_animal
  )

serum_hormone_sample_matches = serum_sample_inventory %>%
  inner_join(
    serum_metabolite_samples,
    by = c("metabolite_sample_id", "metabolite_timepoint")
  ) %>%
  inner_join(
    hormone_samples,
    by = c("hormone_sample_id", "hormone_sampling_time")
  ) %>%
  arrange(metabolite_timepoint, sample, hormone_sample_id)

readr::write_csv(
  serum_hormone_sample_matches %>%
    transmute(
      sample,
      metabolite_sample_id,
      `metabolite timepoint` = metabolite_timepoint,
      `hormone sample` = hormone_sample_id,
      `collection protocol` = collection_protocol,
      `blood type` = blood_type,
      `timepoint months` = timepoint_months
    ),
  file.path(output_dir, "serum hormone metabolomics sample matches.csv"),
  na = "N/A"
)

hormone_data = raw_hormone_data %>%
  filter(hormone_sampling_time %in% c("SAC", "GTT")) %>%
  pivot_longer(
    cols = all_of(hormone_columns),
    names_to = "hormone",
    values_to = "hormone_abundance"
  ) %>%
  filter(
    !(sample == "2M33" & raw_sample_id == "2M33_SAC" & hormone == "Secretin")
  ) %>%
  select(
    hormone_mouse_id = sample,
    hormone_sample_id = raw_sample_id,
    hormone_sampling_time,
    hormone,
    hormone_abundance,
    hormone_cohort,
    hormone_animal
  )


# ========== 2.0 - Join old-animal matched serum metabolite and hormone observations ==========
# --

matched_observations = serum_hormone_sample_matches %>%
  inner_join(
    serum_metabolite_data,
    by = c("metabolite_sample_id", "metabolite_timepoint"),
    relationship = "many-to-many"
  ) %>%
  inner_join(
    hormone_data,
    by = c("hormone_sample_id", "hormone_sampling_time"),
    relationship = "many-to-many"
  ) %>%
  filter(
    metabolite_age == "Old",
    metabolite_diet %in% c("CD", "ID"),
    !is.na(metabolite_abundance),
    !is.na(hormone_abundance)
  ) %>%
  mutate(
    matched_observation_pair = paste0(
      metabolite_timepoint, "_metabolite_vs_",
      hormone_sampling_time, "_hormone"
    )
  )

matched_sample_summary = matched_observations %>%
  distinct(sample, metabolite_diet, metabolite_timepoint, hormone_sampling_time) %>%
  count(metabolite_timepoint, hormone_sampling_time, metabolite_diet, name = "n_matched_old_samples") %>%
  arrange(metabolite_timepoint, hormone_sampling_time, metabolite_diet)

readr::write_csv(
  matched_sample_summary,
  file.path(output_dir, "old serum hormone matched sample summary.csv"),
  na = "N/A"
)


# ========== 3.0 - Define requested correlation sets ==========
# --

analysis_specs = tibble::tibble(
  analysis_id = c(
    "pv_metabolite_vs_sac_hormone",
    "sys_sac_metabolite_vs_sac_hormone",
    "tp4_metabolite_vs_gtt_hormone",
    "tp4_and_sys_sac_metabolite_vs_matched_gtt_and_sac_hormone"
  ),
  analysis_label = c(
    "metabolite PV vs hormone SAC",
    "metabolite sys_sac vs hormone SAC",
    "metabolite TP4 vs hormone GTT",
    "metabolite TP4 and sys_sac vs hormone GTT and SAC"
  ),
  output_file = c(
    "old serum pv metabolite vs sac hormone correlations.csv",
    "old serum sys_sac metabolite vs sac hormone correlations.csv",
    "old serum tp4 metabolite vs gtt hormone correlations.csv",
    "old serum tp4 and sys_sac metabolite vs matched gtt and sac hormone correlations.csv"
  ),
  allowed_pairs = list(
    tibble::tibble(metabolite_timepoint = "pv", hormone_sampling_time = "SAC"),
    tibble::tibble(metabolite_timepoint = "sys_sac", hormone_sampling_time = "SAC"),
    tibble::tibble(metabolite_timepoint = "tp4", hormone_sampling_time = "GTT"),
    tibble::tibble(
      metabolite_timepoint = c("tp4", "sys_sac"),
      hormone_sampling_time = c("GTT", "SAC")
    )
  )
)


# ========== 4.0 - Calculate correlations ==========
# --

add_diet_correlation_groups = function(data) {
  bind_rows(
    data %>%
      filter(metabolite_diet == "CD") %>%
      mutate(diet_correlation_group = "CD"),
    data %>%
      filter(metabolite_diet == "ID") %>%
      mutate(diet_correlation_group = "ID"),
    data %>%
      filter(metabolite_diet %in% c("CD", "ID")) %>%
      mutate(diet_correlation_group = "CD_ID_combined")
  )
}

calculate_correlation_test = function(metabolite_abundance, hormone_abundance, correlation_method) {
  if (
    length(metabolite_abundance) >= minimum_n_for_correlation &&
      dplyr::n_distinct(metabolite_abundance) > 1 &&
      dplyr::n_distinct(hormone_abundance) > 1
  ) {
    tryCatch(
      if (correlation_method == "bicor") {
        WGCNA::bicorAndPvalue(
          metabolite_abundance,
          hormone_abundance,
          use = "pairwise.complete.obs"
        )
      } else {
        cor.test(
          metabolite_abundance,
          hormone_abundance,
          method = correlation_method,
          exact = FALSE
        )
      },
      error = function(error) NULL
    )
  } else {
    NULL
  }
}

extract_correlation = function(test, correlation_method) {
  if (is.null(test)) {
    NA_real_
  } else if (correlation_method == "bicor") {
    unname(test$bicor[1, 1])
  } else {
    unname(test$estimate)
  }
}

extract_p_value = function(test, correlation_method) {
  if (is.null(test)) {
    NA_real_
  } else if (correlation_method == "bicor") {
    unname(test$p[1, 1])
  } else {
    test$p.value
  }
}

calculate_correlations = function(data, current_correlation_method, required_observation_pairs) {
  summarised_correlations = data %>%
    group_by(
      diet_correlation_group,
      analysis_id,
      analysis_label,
      metabolite_source,
      metabolite_ion_mode,
      metabolite,
      metabolite_formula,
      metabolite_mz,
      metabolite_rt,
      hormone
    ) %>%
    summarise(
      matched_observation_pairs = paste(sort(unique(matched_observation_pair)), collapse = ";"),
      n_matched_observation_pairs = n_distinct(matched_observation_pair),
      metabolite_timepoints = paste(sort(unique(metabolite_timepoint)), collapse = ";"),
      hormone_sampling_times = paste(sort(unique(hormone_sampling_time)), collapse = ";"),
      n = n(),
      n_unique_samples = n_distinct(sample),
      n_unique_metabolite_samples = n_distinct(metabolite_sample_id),
      n_unique_hormone_samples = n_distinct(hormone_sample_id),
      n_cd = n_distinct(sample[metabolite_diet == "CD"]),
      n_id = n_distinct(sample[metabolite_diet == "ID"]),
      n_med_adj_ic = sum(metabolite_abundance_source == "med.adj.ic", na.rm = TRUE),
      n_raw_ion_count_fallback = sum(metabolite_abundance_source == "ion.count", na.rm = TRUE),
      metabolite_min = min(metabolite_abundance, na.rm = TRUE),
      metabolite_max = max(metabolite_abundance, na.rm = TRUE),
      hormone_min = min(hormone_abundance, na.rm = TRUE),
      hormone_max = max(hormone_abundance, na.rm = TRUE),
      test = list(calculate_correlation_test(
        metabolite_abundance,
        hormone_abundance,
        current_correlation_method
      )),
      .groups = "drop"
    )

  summarised_correlations %>%
    filter(n_matched_observation_pairs == required_observation_pairs) %>%
    mutate(
      correlation_method = current_correlation_method,
      correlation = purrr::map_dbl(test, ~ extract_correlation(.x, current_correlation_method)),
      p_value = purrr::map_dbl(test, ~ extract_p_value(.x, current_correlation_method))
    ) %>%
    group_by(analysis_id, diet_correlation_group, correlation_method) %>%
    mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    select(
      diet_correlation_group,
      analysis_id,
      analysis_label,
      correlation_method,
      matched_observation_pairs,
      n_matched_observation_pairs,
      metabolite_timepoints,
      hormone_sampling_times,
      metabolite_source,
      metabolite_ion_mode,
      metabolite,
      metabolite_formula,
      metabolite_mz,
      metabolite_rt,
      hormone,
      n,
      n_unique_samples,
      n_unique_metabolite_samples,
      n_unique_hormone_samples,
      n_cd,
      n_id,
      n_med_adj_ic,
      n_raw_ion_count_fallback,
      correlation,
      p_value,
      p_adj_bh,
      metabolite_min,
      metabolite_max,
      hormone_min,
      hormone_max
    ) %>%
    arrange(
      diet_correlation_group,
      correlation_method,
      p_value,
      metabolite_source,
      metabolite_ion_mode,
      metabolite,
      hormone
    )
}

calculate_analysis_correlations = function(spec_row) {
  required_observation_pairs = nrow(spec_row$allowed_pairs[[1]])

  analysis_observations = matched_observations %>%
    semi_join(
      spec_row$allowed_pairs[[1]],
      by = c("metabolite_timepoint", "hormone_sampling_time")
    ) %>%
    mutate(
      analysis_id = spec_row$analysis_id,
      analysis_label = spec_row$analysis_label
    ) %>%
    add_diet_correlation_groups()

  correlations = purrr::map_dfr(
    correlation_methods,
    ~ calculate_correlations(analysis_observations, .x, required_observation_pairs)
  )

  readr::write_csv(
    correlations,
    file.path(output_dir, spec_row$output_file),
    na = "N/A"
  )

  correlations %>%
    mutate(output_file = spec_row$output_file, .after = analysis_label)
}

correlation_results = purrr::map_dfr(
  seq_len(nrow(analysis_specs)),
  ~ calculate_analysis_correlations(analysis_specs[.x, ])
)


# ========== 5.0 - Summarise correlation output sizes and significant results ==========
# --

correlation_summary = correlation_results %>%
  group_by(
    analysis_id,
    analysis_label,
    output_file,
    diet_correlation_group,
    correlation_method
  ) %>%
  summarise(
    n_combinations = n(),
    n_tested = sum(!is.na(p_value)),
    n_significant = sum(!is.na(p_adj_bh) & p_adj_bh < 0.05),
    n_nominal_p_lt_0_05 = sum(!is.na(p_value) & p_value < 0.05),
    .groups = "drop"
  ) %>%
  mutate(significance_definition = "BH adjusted p < 0.05", .after = n_significant) %>%
  arrange(analysis_id, diet_correlation_group, correlation_method)

readr::write_csv(
  correlation_summary,
  file.path(output_dir, "old serum hormone metabolite correlation summary.csv"),
  na = "N/A"
)


# ========== 6.0 - Export compact CD+ID bicor table for TP4/GTT and sys_sac/SAC analysis ==========
# --

tp4_sys_sac_cd_id_bicor_results = correlation_results %>%
  filter(
    analysis_id == "tp4_and_sys_sac_metabolite_vs_matched_gtt_and_sac_hormone",
    diet_correlation_group == "CD_ID_combined",
    correlation_method == "bicor"
  ) %>%
  mutate(
    abs_correlation = abs(correlation),
    correlation_direction = case_when(
      correlation > 0 ~ "positive",
      correlation < 0 ~ "negative",
      TRUE ~ "zero"
    )
  ) %>%
  arrange(p_adj_bh, p_value, desc(abs_correlation), metabolite, hormone) %>%
  mutate(result_rank = row_number(), .before = 1) %>%
  select(
    result_rank,
    analysis_label,
    diet_correlation_group,
    correlation_method,
    metabolite_timepoints,
    hormone_sampling_times,
    metabolite_source,
    metabolite_ion_mode,
    metabolite,
    metabolite_formula,
    metabolite_mz,
    metabolite_rt,
    hormone,
    n,
    n_unique_samples,
    n_unique_metabolite_samples,
    n_unique_hormone_samples,
    n_cd,
    n_id,
    correlation,
    abs_correlation,
    correlation_direction,
    p_value,
    p_adj_bh,
    n_raw_ion_count_fallback,
    metabolite_min,
    metabolite_max,
    hormone_min,
    hormone_max
  ) %>%
  mutate(
    across(
      where(is.character),
      ~ trimws(gsub("[[:space:]]+", " ", .x))
    )
  )

readr::write_csv(
  tp4_sys_sac_cd_id_bicor_results,
  file.path(output_dir, "old serum tp4 sys_sac matched gtt sac CD ID combined bicor essential results.csv"),
  na = "N/A"
)


# ========== 7.0 - Export nominal bicor hits with metabolite diet statistics ==========
# --

safe_mean = function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    NA_real_
  } else {
    mean(x, na.rm = TRUE)
  }
}

calculate_welch_t_test_p_value = function(metabolite_abundance, metabolite_diet) {
  test_data = tibble::tibble(
    metabolite_abundance = metabolite_abundance,
    metabolite_diet = as.character(metabolite_diet)
  ) %>%
    filter(
      !is.na(metabolite_abundance),
      metabolite_diet %in% c("CD", "ID")
    )

  if (
    n_distinct(test_data$metabolite_diet) < 2 ||
      sum(test_data$metabolite_diet == "CD") < 2 ||
      sum(test_data$metabolite_diet == "ID") < 2
  ) {
    return(NA_real_)
  }

  tryCatch(
    t.test(
      metabolite_abundance ~ metabolite_diet,
      data = test_data
    )$p.value,
    error = function(error) NA_real_
  )
}

metabolite_diet_statistics = serum_sample_inventory %>%
  inner_join(
    serum_metabolite_data,
    by = c("metabolite_sample_id", "metabolite_timepoint")
  ) %>%
  filter(
    metabolite_age == "Old",
    metabolite_diet %in% c("CD", "ID"),
    metabolite_timepoint %in% c("tp4", "sys_sac"),
    !is.na(metabolite_abundance)
  ) %>%
  mutate(
    across(
      c(metabolite_source, metabolite_ion_mode, metabolite, metabolite_formula),
      ~ trimws(gsub("[[:space:]]+", " ", .x))
    )
  ) %>%
  group_by(
    metabolite_timepoint,
    metabolite_source,
    metabolite_ion_mode,
    metabolite,
    metabolite_formula,
    metabolite_mz,
    metabolite_rt
  ) %>%
  summarise(
    n_cd = n_distinct(sample[metabolite_diet == "CD"]),
    n_id = n_distinct(sample[metabolite_diet == "ID"]),
    mean_cd = safe_mean(metabolite_abundance[metabolite_diet == "CD"]),
    mean_id = safe_mean(metabolite_abundance[metabolite_diet == "ID"]),
    fc_id_vs_cd = if_else(!is.na(mean_cd) & mean_cd > 0, mean_id / mean_cd, NA_real_),
    log2_fc_id_vs_cd = if_else(
      !is.na(mean_cd) & mean_cd > 0 & !is.na(mean_id) & mean_id > 0,
      log2(mean_id / mean_cd),
      NA_real_
    ),
    ttest_p_value = calculate_welch_t_test_p_value(
      metabolite_abundance,
      metabolite_diet
    ),
    .groups = "drop"
  ) %>%
  mutate(
    metabolite_timepoint = recode(
      metabolite_timepoint,
      tp4 = "gtt",
      sys_sac = "sac"
    )
  ) %>%
  pivot_wider(
    id_cols = c(
      metabolite_source,
      metabolite_ion_mode,
      metabolite,
      metabolite_formula,
      metabolite_mz,
      metabolite_rt
    ),
    names_from = metabolite_timepoint,
    values_from = c(
      n_cd,
      n_id,
      mean_cd,
      mean_id,
      fc_id_vs_cd,
      log2_fc_id_vs_cd,
      ttest_p_value
    ),
    names_glue = "metabolite_{metabolite_timepoint}_{.value}"
  )

tp4_sys_sac_cd_id_bicor_nominal_results = tp4_sys_sac_cd_id_bicor_results %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  left_join(
    metabolite_diet_statistics,
    by = c(
      "metabolite_source",
      "metabolite_ion_mode",
      "metabolite",
      "metabolite_formula",
      "metabolite_mz",
      "metabolite_rt"
    )
  ) %>%
  arrange(p_value, p_adj_bh, desc(abs_correlation), metabolite, hormone) %>%
  mutate(nominal_result_rank = row_number(), .before = 1)

readr::write_csv(
  tp4_sys_sac_cd_id_bicor_nominal_results,
  file.path(
    output_dir,
    "old serum tp4 sys_sac matched gtt sac CD ID combined bicor nominal p significant with metabolite diet stats.csv"
  ),
  na = "N/A"
)
