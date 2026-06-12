# ****** Inulin Aging Mouse Cytokines
# ****** Spring 2026 - Ian T, Jessie K, Randa L


here::i_am("code/cytokines/01_process_cytokines.R")


# ========== 0.0 - Environment ==========
# --
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

dir.create(here::here("data/processed/cytokines"), recursive = TRUE, showWarnings = FALSE)


# ========== 1.0 - Read cytokine summary and metadata ==========
# --
cytokine_file = here::here(
  "data/raw/cytokines",
  "RANDA_20260102_Inulin_Aging_inflammation__longitudinal_Custom_panel_Luminex.xlsx"
)

cytokine_raw = readxl::read_xlsx(cytokine_file, sheet = 1)

metadata = readr::read_csv(
  here::here("data/libraries/metadata.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  mutate(
    sample = str_trim(sample),
    cohort = as.character(cohort)
  )


# ========== 2.0 - Parse sample identifiers ==========
# --
cytokine_samples = cytokine_raw %>%
  rename(
    plate = Plate,
    assay_group = Group,
    well_id = `Well ID`,
    raw_sample_id = `Sample ID`
  ) %>%
  mutate(
    raw_sample_id = str_trim(raw_sample_id),
    sample = str_extract(raw_sample_id, "^[^_]+"),
    cohort = str_extract(sample, "^\\d+"),
    animal = as.integer(str_extract(sample, "\\d+$")),
    sample_suffix = str_remove(raw_sample_id, "^[^_]+_?"),
    sampling_tissue = case_when(
      str_detect(sample_suffix, "(^|_)sys($|_)") ~ "systemic blood",
      str_detect(sample_suffix, "(^|_)pv($|_)") ~ "portal blood",
      str_detect(sample_suffix, "^\\d+mon$") ~ "systemic blood",
      TRUE ~ NA_character_
    ),
    sampling_time = case_when(
      str_detect(sample_suffix, "(^|_)sac($|_)") ~ "sacrificed sample",
      str_detect(sample_suffix, "^\\d+mon$") ~ "intermediary month sample",
      TRUE ~ NA_character_
    ),
    sampling_month = as.integer(str_match(sample_suffix, "^(\\d+)mon$")[, 2])
  )


# ========== 3.0 - Tidy cytokine abundance values ==========
# --
cytokine_tidy = cytokine_samples %>%
  pivot_longer(
    cols = matches("\\.(con|MFI)$"),
    names_to = c("cytokine", "measurement"),
    names_pattern = "^(.*)\\.(con|MFI)$",
    values_to = "value_raw",
    values_transform = list(value_raw = as.character)
  ) %>%
  mutate(
    measurement = recode(
      measurement,
      con = "concentration_pg_ml",
      MFI = "mfi"
    ),
    value_clean = if_else(value_raw %in% c("-", "", "NA", "N/A"), NA_character_, value_raw),
    value = as.numeric(value_clean)
  ) %>%
  select(
    plate, assay_group, well_id,
    raw_sample_id, sample, sample_suffix,
    sampling_tissue, sampling_time, sampling_month,
    cohort, animal,
    cytokine, measurement, value
  ) %>%
  pivot_wider(
    names_from = measurement,
    values_from = value
  ) %>%
  # Match diet, age, and group using the metadata mouse ID in the sample column.
  left_join(
    metadata %>%
      select(
        sample,
        diet,
        cohort_metadata = cohort,
        animal_metadata = animal,
        age,
        group
      ),
    by = "sample"
  ) %>%
  mutate(
    cohort = coalesce(cohort, cohort_metadata),
    animal = coalesce(animal, animal_metadata)
  ) %>%
  select(
    plate, assay_group, well_id,
    raw_sample_id, sample, sample_suffix,
    sampling_tissue, sampling_time, sampling_month,
    diet, cohort, animal, age, group,
    cytokine, concentration_pg_ml, mfi
  ) %>%
  arrange(sample, sampling_time, sampling_month, sampling_tissue, cytokine)


# ========== 4.0 - Export tidy cytokine table ==========
# --
readr::write_csv(
  cytokine_tidy,
  here::here("data/processed/cytokines/tidy cytokine abundances.csv"),
  na = "N/A"
)
