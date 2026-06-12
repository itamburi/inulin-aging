# ****** Inulin Aging Hormone-Metabolomics Sample Mapping Audit
# ****** Summer 2026 - Ian T, Jessie K, Randa L
#
# Builds a collaborator-facing crosswalk between hormone sample IDs and the
# tissue/serum metabolomics sample IDs used in correlation analyses.

here::i_am("code/integrative/00_hormone_metabolomics_sample_mapping_audit.R")


# ========== 0.0 - Environment ==========
# --
library(dplyr)
library(readr)
library(stringr)

output_dir = here::here("data/processed/integrative/sample_mapping_audits")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ========== 1.0 - Helpers and metadata ==========
# --
normalize_mouse_id = function(sample_id) {
  parsed = stringr::str_match(sample_id, "^(\\d+)M_?0*(\\d+)$")

  dplyr::if_else(
    is.na(parsed[, 1]),
    sample_id,
    paste0(parsed[, 2], "M", as.integer(parsed[, 3]))
  )
}

canonical_metadata = readr::read_csv(
  here::here("data/libraries/metadata.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    sample = as.character(sample),
    metadata_source = "Metabolomics/data/libraries/metadata.csv",
    metadata_review_status = "present in current metabolomics metadata",
    metadata_note = NA_character_
  )

supplemental_metadata = tibble::tibble(
  sample = c(paste0("3M", 1:40), paste0("5M", 1:10)),
  diet = c(rep("ID", 40), rep("CD", 10)),
  cohort = c(rep(3, 40), rep(5, 10)),
  animal = c(1:40, 1:10),
  age = c(rep("Old", 40), rep("Young", 10))
) %>%
  mutate(
    group = paste(age, diet),
    metadata_source = "supplemental proposed mapping for serum/hormone audit",
    metadata_review_status = "needs collaborator confirmation",
    metadata_note = dplyr::case_when(
      cohort == 3 ~ "Cohort 3 assigned as Old ID in current serum correlation script; not present in metabolomics metadata.csv.",
      cohort == 5 ~ "Cohort 5 assigned as Young CD in current serum correlation script; 5M1-5M4 are also Young CD in RNAseq metadata.",
      TRUE ~ NA_character_
    )
  ) %>%
  anti_join(canonical_metadata, by = "sample")

sample_metadata = bind_rows(canonical_metadata, supplemental_metadata)


# ========== 2.0 - Input sample manifests ==========
# --
hormone_samples = readr::read_csv(
  here::here("data/raw/hormones/luminex_hormone_inulin_aging_raw.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "-")
) %>%
  transmute(
    hormone_sample_id = stringr::str_trim(`Sample ID`),
    mouse_id = stringr::str_extract(hormone_sample_id, "^[^_]+"),
    hormone_sampling_time = stringr::str_extract(hormone_sample_id, "(?<=_).*"),
    hormone_plate = Plate,
    hormone_assay_group = Group,
    hormone_well_id = `Well ID`,
    hormone_location = Location
  ) %>%
  distinct()

serum_metabolomics_samples = readr::read_csv(
  here::here("data/processed/serum_metabolomics/ion_counts/tidy and normalized serum ion count data.csv"),
  show_col_types = FALSE
) %>%
  distinct(
    metabolomics_domain = tissue,
    metabolomics_tissue = tissue,
    metabolomics_original_sample = sample,
    metabolomics_timepoint = timepoint,
    metabolomics_sample_id = sample.id,
    metabolomics_sequence_id = sequence.id,
    metabolomics_source = source,
    metabolomics_ion_mode = ion.mode,
    metabolomics_raw_file = raw.file
  ) %>%
  mutate(mouse_id = normalize_mouse_id(metabolomics_original_sample)) %>%
  group_by(
    metabolomics_domain, metabolomics_tissue, metabolomics_original_sample,
    metabolomics_timepoint, metabolomics_sample_id, mouse_id
  ) %>%
  summarise(
    metabolomics_sequence_ids = paste(sort(unique(metabolomics_sequence_id)), collapse = ";"),
    metabolomics_sources = paste(sort(unique(metabolomics_source)), collapse = ";"),
    metabolomics_ion_modes = paste(sort(unique(metabolomics_ion_mode)), collapse = ";"),
    metabolomics_raw_files = paste(sort(unique(metabolomics_raw_file)), collapse = ";"),
    .groups = "drop"
  )

tissue_metabolomics_samples = readr::read_csv(
  here::here("data/processed/tissue_metabolomics/ion_counts/tidy and normalized tissue ion count data.csv"),
  show_col_types = FALSE
) %>%
  distinct(
    metabolomics_domain = tissue,
    metabolomics_tissue = tissue,
    metabolomics_original_sample = sample,
    metabolomics_sample_id = sample.id,
    metabolomics_sequence_id = sequence.id,
    metabolomics_source = source,
    metabolomics_ion_mode = ion.mode,
    metabolomics_raw_file = raw.file
  ) %>%
  mutate(
    metabolomics_timepoint = "SAC",
    mouse_id = normalize_mouse_id(metabolomics_original_sample)
  ) %>%
  group_by(
    metabolomics_domain, metabolomics_tissue, metabolomics_original_sample,
    metabolomics_timepoint, metabolomics_sample_id, mouse_id
  ) %>%
  summarise(
    metabolomics_sequence_ids = paste(sort(unique(metabolomics_sequence_id)), collapse = ";"),
    metabolomics_sources = paste(sort(unique(metabolomics_source)), collapse = ";"),
    metabolomics_ion_modes = paste(sort(unique(metabolomics_ion_mode)), collapse = ";"),
    metabolomics_raw_files = paste(sort(unique(metabolomics_raw_file)), collapse = ";"),
    .groups = "drop"
  )

metabolomics_samples = bind_rows(serum_metabolomics_samples, tissue_metabolomics_samples)


# ========== 3.0 - Build expected hormone-metabolomics pairings ==========
# --
make_pairing_map = function(metabolomics_domain_filter,
                            metabolomics_tissue_filter,
                            metabolomics_timepoint_filter,
                            hormone_sampling_time_filter,
                            pairing_label) {
  metabolomics_for_pairing = metabolomics_samples %>%
    filter(
      metabolomics_domain == metabolomics_domain_filter,
      metabolomics_tissue == metabolomics_tissue_filter,
      metabolomics_timepoint == metabolomics_timepoint_filter
    ) %>%
    mutate(
      pairing = pairing_label,
      expected_hormone_sampling_time = hormone_sampling_time_filter
    )

  hormone_for_pairing = hormone_samples %>%
    filter(hormone_sampling_time == hormone_sampling_time_filter) %>%
    mutate(pairing = pairing_label)

  full_join(
    metabolomics_for_pairing,
    hormone_for_pairing,
    by = c("pairing", "mouse_id")
  ) %>%
    mutate(
      metabolomics_domain = coalesce(metabolomics_domain, metabolomics_domain_filter),
      metabolomics_tissue = coalesce(metabolomics_tissue, metabolomics_tissue_filter),
      metabolomics_timepoint = coalesce(metabolomics_timepoint, metabolomics_timepoint_filter),
      expected_hormone_sampling_time = coalesce(
        expected_hormone_sampling_time,
        hormone_sampling_time_filter
      ),
      mapping_status = dplyr::case_when(
        !is.na(metabolomics_sample_id) & !is.na(hormone_sample_id) ~ "matched",
        !is.na(metabolomics_sample_id) & is.na(hormone_sample_id) ~ "metabolomics_only_no_hormone",
        is.na(metabolomics_sample_id) & !is.na(hormone_sample_id) ~ "hormone_only_no_metabolomics",
        TRUE ~ "unresolved"
      )
    )
}

sample_mapping = bind_rows(
  make_pairing_map("Serum", "Serum", "pv", "SAC", "serum pv metabolites vs SAC hormones"),
  make_pairing_map("Serum", "Serum", "sys_sac", "SAC", "serum systemic SAC metabolites vs SAC hormones"),
  make_pairing_map("Serum", "Serum", "tp4", "GTT", "serum TP4 metabolites vs GTT hormones"),
  make_pairing_map("Liver", "Liver", "SAC", "SAC", "liver tissue metabolites vs SAC hormones"),
  make_pairing_map("Muscle", "Muscle", "SAC", "SAC", "muscle tissue metabolites vs SAC hormones")
) %>%
  left_join(
    sample_metadata %>%
      select(
        mouse_id = sample,
        diet, age, cohort, animal, group,
        metadata_source, metadata_review_status, metadata_note
      ),
    by = "mouse_id"
  ) %>%
  mutate(
    needs_metadata_review = metadata_review_status == "needs collaborator confirmation" |
      is.na(metadata_review_status),
    metadata_review_status = coalesce(metadata_review_status, "missing from current metadata mapping"),
    metadata_note = dplyr::case_when(
      is.na(metadata_note) & is.na(diet) ~ "No diet/age metadata found for this mouse ID.",
      TRUE ~ metadata_note
    )
  ) %>%
  select(
    pairing,
    mapping_status,
    mouse_id,
    hormone_sample_id,
    hormone_sampling_time,
    expected_hormone_sampling_time,
    metabolomics_sample_id,
    metabolomics_domain,
    metabolomics_tissue,
    metabolomics_timepoint,
    metabolomics_original_sample,
    metabolomics_sequence_ids,
    metabolomics_sources,
    metabolomics_ion_modes,
    metabolomics_raw_files,
    diet, age, cohort, animal, group,
    metadata_review_status,
    needs_metadata_review,
    metadata_source,
    metadata_note,
    hormone_plate,
    hormone_assay_group,
    hormone_well_id,
    hormone_location
  ) %>%
  arrange(
    pairing,
    factor(mapping_status, levels = c("matched", "metabolomics_only_no_hormone", "hormone_only_no_metabolomics", "unresolved")),
    cohort,
    animal,
    mouse_id
  )

sample_mapping_summary = sample_mapping %>%
  count(pairing, mapping_status, metadata_review_status, name = "n_rows") %>%
  arrange(pairing, mapping_status, metadata_review_status)

sample_mapping_for_collaborator_review = sample_mapping %>%
  filter(mapping_status == "matched") %>%
  mutate(
    metabolomics_tissue = dplyr::case_when(
      metabolomics_domain == "Serum" & metabolomics_timepoint == "pv" ~ "Serum PV",
      metabolomics_domain == "Serum" & metabolomics_timepoint == "sys_sac" ~ "Serum systemic SAC",
      metabolomics_domain == "Serum" & metabolomics_timepoint == "tp4" ~ "Serum TP4",
      TRUE ~ metabolomics_tissue
    ),
    hormone_tissue = paste("Blood", hormone_sampling_time)
  ) %>%
  select(
    cohort,
    diet,
    age,
    hormone_sample_id,
    metabolomics_sample_id,
    metabolomics_tissue,
    hormone_tissue
  ) %>%
  arrange(metabolomics_tissue, hormone_tissue, hormone_sample_id, metabolomics_sample_id)


# ========== 4.0 - Export ==========
# --
readr::write_csv(
  sample_mapping_for_collaborator_review,
  file.path(output_dir, "hormone to metabolomics sample mapping.csv"),
  na = "N/A"
)

readr::write_csv(
  sample_mapping_summary,
  file.path(output_dir, "hormone to metabolomics serum and tissue sample mapping summary.csv"),
  na = "N/A"
)
