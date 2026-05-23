# ****** Inulin Aging Mouse Cytokine-Metabolite Correlations
# ****** Spring 2026 - Ian T, Jessie K, Randa L


here::i_am("code/04 Cytokine metabolite correlations.R")


# ========== 0.0 - Environment ==========
# --
library(dplyr)
library(tidyr)
library(readr)
library(purrr)

dir.create(here::here("data/processed/cytokines"), recursive = TRUE, showWarnings = FALSE)


# ========== 1.0 - Read processed cytokine and metabolite data ==========
# --
metabolite_data = readr::read_csv(
  here::here("data/processed/all group sample-level median adjusted ion counts.csv"),
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

cytokine_data = readr::read_csv(
  here::here("data/processed/cytokines/tidy cytokine abundances.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  filter(sampling_time == "sacrificed sample") %>%
  select(
    sample, cytokine_sampling_tissue = sampling_tissue,
    cytokine, concentration_pg_ml, mfi,
    cytokine_group = group, cytokine_cohort = cohort,
    cytokine_age = age, cytokine_diet = diet,
    cytokine_animal = animal
  ) %>%
  pivot_longer(
    cols = c(concentration_pg_ml, mfi),
    names_to = "cytokine_measurement",
    values_to = "cytokine_abundance"
  )


# ========== 2.0 - Join matching mouse-level observations ==========
# --
matched_observations = metabolite_data %>%
  inner_join(cytokine_data, by = "sample", relationship = "many-to-many") %>%
  filter(
    !is.na(metabolite_abundance),
    !is.na(cytokine_abundance)
  )

readr::write_csv(
  matched_observations,
  here::here("data/processed/cytokines/matched sacrifice cytokine metabolite observations.csv"),
  na = "N/A"
)


# ========== 3.0 - Calculate Pearson and Spearman correlations ==========
# --
calculate_correlations = function(data, correlation_method) {
  data %>%
    group_by(
      metabolite_tissue, metabolite,
      cytokine_sampling_tissue, cytokine,
      cytokine_measurement
    ) %>%
    summarise(
      n = n(),
      n_unique_samples = n_distinct(sample),
      metabolite_min = min(metabolite_abundance, na.rm = TRUE),
      metabolite_max = max(metabolite_abundance, na.rm = TRUE),
      cytokine_min = min(cytokine_abundance, na.rm = TRUE),
      cytokine_max = max(cytokine_abundance, na.rm = TRUE),
      test = list(
        if (
          n() >= 3 &&
            dplyr::n_distinct(metabolite_abundance) > 1 &&
            dplyr::n_distinct(cytokine_abundance) > 1
        ) {
          cor.test(
            metabolite_abundance,
            cytokine_abundance,
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
      cytokine_sampling_tissue, cytokine, cytokine_measurement,
      n, n_unique_samples,
      correlation, p_value, p_adj_bh,
      metabolite_min, metabolite_max,
      cytokine_min, cytokine_max
    ) %>%
    arrange(p_value, metabolite_tissue, metabolite, cytokine_sampling_tissue, cytokine)
}

pearson_correlations = calculate_correlations(matched_observations, "pearson")

spearman_correlations = calculate_correlations(matched_observations, "spearman")


# ========== 4.0 - Export correlation tables ==========
# --
readr::write_csv(
  pearson_correlations,
  here::here("data/processed/cytokines/pearson sacrifice cytokine metabolite correlations.csv"),
  na = "N/A"
)

readr::write_csv(
  spearman_correlations,
  here::here("data/processed/cytokines/spearman sacrifice cytokine metabolite correlations.csv"),
  na = "N/A"
)
