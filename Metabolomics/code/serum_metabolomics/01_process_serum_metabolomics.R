# ****** Inulin Aging Mouse Serum Metabolomics
# ****** Summer 2026 - Ian T, Jessie K, Randa L


here::i_am("code/serum_metabolomics/01_process_serum_metabolomics.R")
source(here::here("code/z-source.R"))


# ========== 0.0 - Environment ==========
# --

raw_cd_dir = here::here("data/raw/serum/compound_discoverer")
raw_maven_dir = here::here("data/raw/serum/maven")
processed_dir = here::here("data/processed/serum_metabolomics/ion_counts")
maven_merger_dir = here::here("data/processed/serum_metabolomics/ion_counts/maven_merger")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(maven_merger_dir, recursive = TRUE, showWarnings = FALSE)

# Parse positive/negative ion mode from CD filenames or serum sequence IDs.
get_serum_ion_mode_from_text <- function(x) {
  dplyr::case_when(
    grepl("positive|(?:^|_)pos(?:_|$)", x, ignore.case = TRUE, perl = TRUE) ~ "posi",
    grepl("negative|(?:^|_)neg(?:2)?(?:_|$)", x, ignore.case = TRUE, perl = TRUE) ~ "nega",
    TRUE ~ NA_character_
  )
}

# Remove CD/Maven suffixes so reruns and timestamped injections share one ID.
normalize_serum_sequence_id <- function(sequence_id) {
  sequence_id %>%
    stringr::str_remove("\\.raw.*$") %>%
    stringr::str_replace("_+$", "") %>%
    stringr::str_replace(
      "(_(?:pos|neg))(?:_(?:rerun|resume)[0-9]*|_[0-9]{8,})+$",
      "\\1"
    ) %>%
    stringr::str_replace("_+$", "")
}

# Extract the mouse sample ID, such as 3M21, from a serum sequence ID.
get_serum_sample <- function(sequence_id) {
  stringr::str_extract(sequence_id, "^[1-5]M\\d+")
}

# Extract the serum collection timepoint encoded in the sequence ID.
get_serum_timepoint <- function(sequence_id) {
  dplyr::case_when(
    grepl("_pv_", sequence_id, ignore.case = TRUE) ~ "pv",
    grepl("_sys_sac_", sequence_id, ignore.case = TRUE) ~ "sys_sac",
    grepl("_tp4_", sequence_id, ignore.case = TRUE) ~ "tp4",
    TRUE ~ NA_character_
  )
}

# Extract numeric cohort from the mouse sample ID.
get_serum_cohort <- function(sample) {
  stringr::str_extract(sample, "^\\d+")
}

# Prefer CD Group Area columns for sample ion counts; fall back to Area columns only if needed.
get_serum_cd_sample_cols <- function(col_names) {
  group_area_cols = col_names[
    grepl("^Group Area: ", col_names) &
      !grepl("blank|STD", col_names, ignore.case = TRUE)
  ]

  if (length(group_area_cols) > 0) {
    return(group_area_cols)
  }

  col_names[
    grepl("^Area: ", col_names) &
      col_names != "Area (Max.)" &
      !grepl("blank|STD", col_names, ignore.case = TRUE)
  ]
}

# Convert a CD measurement column name into the normalized serum sequence ID.
clean_serum_cd_sequence_id <- function(col_name) {
  col_name %>%
    stringr::str_remove("^Group Area: ") %>%
    stringr::str_remove("^Area: ") %>%
    normalize_serum_sequence_id()
}

# Select Maven sample columns while excluding blanks and standards.
get_serum_maven_sample_cols <- function(col_names) {
  col_names[
    grepl(
      "^[1-5]M\\d+_(?:pv|sys_sac|tp4)_(?:pos|neg)(?:_(?:rerun|resume)[0-9]*|_[0-9]{8,})?$",
      col_names,
      ignore.case = TRUE,
      perl = TRUE
    ) &
      !grepl("blank|STD", col_names, ignore.case = TRUE)
  ]
}

# Standardize Maven compound names, especially labeled valine used for normalization.
normalize_maven_compound_name <- function(compound, isotope_label) {
  dplyr::case_when(
    tolower(compound) == "valine" & isotope_label == "C13N15-label-5-1" ~ "Valine-C13N15",
    compound %in% c("L-Valine", "m37_Valine") & isotope_label == "C13N15-label-5-1" ~ "Valine-C13N15",
    compound %in% c("L-Valine", "m37_Valine") ~ "Valine",
    TRUE ~ compound
  )
}


# ========== 1.0 - Convert xlsx to csv ==========
# --

CD_files = list.files(raw_cd_dir, pattern = "\\.xlsx$", full.names = TRUE)

for (CD_file in CD_files) {
  file_name = tools::file_path_sans_ext(basename(CD_file))
  data = readxl::read_xlsx(CD_file)
  csv_file = file.path(raw_cd_dir, paste0(file_name, ".csv"))

  write.csv(data, csv_file, row.names = FALSE)
  cat("Converted", CD_file, "to", csv_file, "\n")
}


# ========== 2.0 - Format Compound Discoverer Data ==========
# --

csv_files = list.files(raw_cd_dir, pattern = "\\.csv$", full.names = TRUE)

csv_data_list = list()
for (csv_file in csv_files) {
  csv_data = read.csv(csv_file, check.names = FALSE)
  csv_data_list[[csv_file]] = csv_data
}

tidy_df_list <- function(filename) {
  cd_df = csv_data_list[[filename]]
  file_label = basename(filename)

  measurement_cols = get_serum_cd_sample_cols(names(cd_df))

  if (length(measurement_cols) == 0) {
    stop(paste0("No serum ion-count columns found for file: ", file_label), call. = FALSE)
  }

  cleaned_sequence_ids = clean_serum_cd_sequence_id(measurement_cols)

  id_cols = intersect(
    c(
      "Name", "Formula", "Calc. MW", "m/z", "RT [min]",
      "Annot. Source: mzCloud Search", "Annot. Source: MassList Search",
      "mzCloud Best Match Confidence", "Annot. DeltaMass [ppm]"
    ),
    names(cd_df)
  )

  tidy = cd_df %>%
    select(all_of(c(id_cols, measurement_cols)))

  names(tidy) = dplyr::recode(
    names(tidy),
    "Calc. MW" = "calc.mw",
    "m/z" = "mz",
    "RT [min]" = "rt.min",
    "Annot. Source: mzCloud Search" = "annot.source.mzcloud",
    "Annot. Source: MassList Search" = "annot.source.masslist",
    "mzCloud Best Match Confidence" = "mzcloud.conf",
    "Annot. DeltaMass [ppm]" = "delta.mass.ppm",
    .default = names(tidy)
  )

  tidy = tidy %>%
    pivot_longer(
      cols = all_of(measurement_cols),
      names_to = "sequence.id",
      values_to = "ion.count"
    ) %>%
    mutate(
      sequence.id = cleaned_sequence_ids[match(sequence.id, measurement_cols)],
      Formula = gsub(" ", "", Formula),
      raw.file = file_label
    ) %>%
    rename_with(tolower) %>%
    mutate(
      ion.count = as.numeric(ion.count),
      mz = as.numeric(mz),
      rt.min = as.numeric(rt.min),
      delta.mass.ppm = as.numeric(delta.mass.ppm),
      mzcloud.conf = as.numeric(mzcloud.conf)
    )

  print(paste0("Completed ", filename))

  return(tidy)
}

tidy_list = lapply(names(csv_data_list), tidy_df_list)
tidy_df = bind_rows(tidy_list)

write.csv(tidy_df, here("data/processed/serum_metabolomics/ion_counts/raw CD observations combined and melted.csv"), row.names = FALSE)

tidy_df2 = tidy_df %>%
  filter(
    is.na(name) == FALSE,
    grepl("\\[Similar to:", name) == FALSE,
    grepl("blank|STD", sequence.id, ignore.case = TRUE) == FALSE,
    abs(delta.mass.ppm) <= 10
  )

write.csv(tidy_df2, here("data/processed/serum_metabolomics/ion_counts/filtered CD observations combined and melted.csv"), row.names = FALSE)


# ========== 3.0 - Check Annotations Against Reference List, Apply mz Confidence Score Filtering ==========
# --

expdata = read.csv(
  here("data/processed/serum_metabolomics/ion_counts/filtered CD observations combined and melted.csv"),
  fileEncoding = "UTF-8-BOM"
) %>%
  mutate(name.exp = tolower(name))

cpds = expdata %>%
  distinct(raw.file, name, name.exp, formula, annot.source.mzcloud, annot.source.masslist, delta.mass.ppm, calc.mw, mz, rt.min, mzcloud.conf) %>%
  select(name, name.exp, everything())

reference_list = read.csv(here("data/libraries/metabolites list_230424_ForCDprocessing.csv"), fileEncoding = "UTF-8-BOM") %>%
  rename_with(~ tolower(.)) %>%
  select("hmdb.name", "rt") %>%
  rename(c("rt.ref" = "rt", "name.ref" = "hmdb.name")) %>%
  mutate(name.ref = tolower(trimws(name.ref)))

cpd_key = c(
  "raw.file", "name", "name.exp", "formula", "annot.source.mzcloud",
  "annot.source.masslist", "delta.mass.ppm", "calc.mw", "mz", "rt.min",
  "mzcloud.conf"
)

# Reference names can have multiple acceptable RTs, so summarize to one row per observed feature.
reference_match_summary = cpds %>%
  left_join(
    reference_list,
    by = c("name.exp" = "name.ref"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    rt.diff = abs(rt.min - rt.ref),
    rt.match = !is.na(rt.diff) & rt.diff < 3
  ) %>%
  group_by(across(all_of(cpd_key))) %>%
  summarise(
    matched.reference.name = any(!is.na(rt.ref)),
    rt.match = any(rt.match),
    min.rt.diff = ifelse(all(is.na(rt.diff)), NA_real_, min(rt.diff, na.rm = TRUE)),
    .groups = "drop"
  )

cpds2 = cpds %>%
  left_join(reference_match_summary, by = cpd_key) %>%
  filter(
    !matched.reference.name | rt.match
  )

mz_filt = cpds2 %>%
  filter(
    !matched.reference.name,
    mzcloud.conf < 40
  )

cpds3 = anti_join(cpds2, mz_filt, by = names(cpds2))

expdata2 = inner_join(
  expdata,
  cpds3,
  by = cpd_key,
  relationship = "many-to-one"
) %>%
  distinct(raw.file, sequence.id, name, formula, mz, rt.min, ion.count, .keep_all = TRUE)

write.csv(
  expdata2,
  here("data/processed/serum_metabolomics/ion_counts/CD data filtered by mismatches to ref and mz cloud conf.csv"),
  row.names = FALSE
)


# ========== 4.0 - Resolve redundant compound entries / select top isomer ==========
# --

expdata = read.csv(
  here("data/processed/serum_metabolomics/ion_counts/CD data filtered by mismatches to ref and mz cloud conf.csv")
) %>%
  mutate(
    sequence.id = normalize_serum_sequence_id(sequence.id),
    sample = get_serum_sample(sequence.id),
    cohort = get_serum_cohort(sample),
    timepoint = get_serum_timepoint(sequence.id),
    tissue = "Serum",
    ion.mode = get_serum_ion_mode_from_text(sequence.id)
  )

batch_list = expdata %>%
  group_by(raw.file) %>%
  group_split(.)

expdata2 = bind_rows(lapply(batch_list, select_top_isomer3))

expdata3 = expdata2 %>%
  filter(grepl("(:| )isomer", name) == FALSE)

expdata4 = expdata3 %>%
  group_by(raw.file, tissue, sequence.id, sample, cohort, timepoint, ion.mode, name, formula) %>%
  arrange(rt.min, .by_group = TRUE) %>%
  summarise(
    rt = first(rt.min),
    mz = first(mz),
    ion.count = max(ion.count, na.rm = TRUE),
    .groups = "drop"
  )

expdata5 = expdata4 %>%
  select(raw.file, tissue, ion.mode, sequence.id, sample, cohort, timepoint, name, formula, rt, mz, ion.count)

write.csv(
  expdata5,
  here("data/processed/serum_metabolomics/ion_counts/annotated and dereplicated CD obs without repeat isomers.csv"),
  row.names = FALSE
)


# ========== 5.0 - Consolidate Negative and Positive Ion Modes ==========
# --

expdata = read.csv(
  here("data/processed/serum_metabolomics/ion_counts/annotated and dereplicated CD obs without repeat isomers.csv"),
  fileEncoding = "UTF-8-BOM"
)

ic_medians = expdata %>%
  group_by(name, ion.mode, formula) %>%
  summarise(median.ic = median(ion.count), .groups = "drop") %>%
  group_by(name) %>%
  mutate(n = n())

top_mode = ic_medians %>%
  group_by(name) %>%
  top_n(1, median.ic)

expdata2 = inner_join(expdata, top_mode %>% select(-median.ic, -n))

expdata3 = expdata2 %>%
  select(raw.file, tissue, ion.mode, sequence.id, sample, cohort, timepoint, name, formula, rt, mz, ion.count)

write.csv(
  expdata3,
  here("data/processed/serum_metabolomics/ion_counts/annotated and processed CD obs consolidated by top ion mode.csv"),
  row.names = FALSE
)


# ========== 6.0 - Merge and Consolidate CD and Maven Observations ==========
# --

mav_files = list.files(raw_maven_dir, pattern = "\\.csv$", full.names = TRUE)

mav_obs_list = lapply(mav_files, function(csv_file) {
  csv_data = read.csv(csv_file, check.names = FALSE)
  sample_cols = get_serum_maven_sample_cols(names(csv_data))

  csv_data %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sequence.id",
      values_to = "ion.count"
    ) %>%
    transmute(
      raw.file = basename(csv_file),
      tissue = "Serum",
      sequence.id = normalize_serum_sequence_id(sequence.id),
      ion.mode = get_serum_ion_mode_from_text(sequence.id),
      sample = get_serum_sample(sequence.id),
      cohort = get_serum_cohort(sample),
      timepoint = get_serum_timepoint(sequence.id),
      isotope.label = isotopeLabel,
      name = normalize_maven_compound_name(compound, isotopeLabel),
      formula = gsub(" ", "", formula),
      medrt = as.numeric(medRt),
      medmz = as.numeric(medMz),
      ion.count = as.numeric(ion.count)
    ) %>%
    filter(
      is.na(sequence.id) == FALSE,
      is.na(ion.mode) == FALSE,
      is.na(sample) == FALSE,
      is.na(timepoint) == FALSE,
      is.na(ion.count) == FALSE
    ) %>%
    distinct()
})

mav_obs = bind_rows(mav_obs_list) %>%
  group_by(tissue, raw.file, sequence.id, ion.mode, sample, cohort, timepoint, name, formula) %>%
  summarise(
    medrt = median(medrt, na.rm = TRUE),
    medmz = median(medmz, na.rm = TRUE),
    ion.count = max(ion.count, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(mav_obs, here("data/processed/serum_metabolomics/ion_counts/maven_merger/all maven observations.csv"), row.names = FALSE)

mav_compounds = mav_obs %>%
  group_by(ion.mode, name, formula) %>%
  summarise(
    source = "Maven",
    mz = median(medmz, na.rm = TRUE),
    rt = median(medrt, na.rm = TRUE),
    med.ic = median(ion.count, na.rm = TRUE),
    n.samples = n(),
    .groups = "drop"
  ) %>%
  select(source, ion.mode, name, formula, mz, rt, med.ic, n.samples)

write.csv(
  mav_compounds,
  here("data/processed/serum_metabolomics/ion_counts/maven_merger/unique compounds and formulas from maven.csv"),
  row.names = FALSE
)

cd_data = read.csv(
  here("data/processed/serum_metabolomics/ion_counts/annotated and processed CD obs consolidated by top ion mode.csv"),
  fileEncoding = "UTF-8-BOM"
) %>%
  mutate(
    source = "CD",
    sequence.id = normalize_serum_sequence_id(sequence.id),
    sample = get_serum_sample(sequence.id),
    cohort = get_serum_cohort(sample),
    timepoint = get_serum_timepoint(sequence.id),
    tissue = "Serum"
  ) %>%
  group_by(tissue, raw.file, sequence.id, ion.mode, sample, cohort, timepoint, name, formula, mz, rt, source) %>%
  summarise(
    ion.count = max(ion.count, na.rm = TRUE),
    .groups = "drop"
  )

cd_compounds = cd_data %>%
  group_by(ion.mode, name, formula, mz, rt) %>%
  summarise(
    source = "CD",
    med.ic = median(ion.count, na.rm = TRUE),
    n.samples = n(),
    .groups = "drop"
  ) %>%
  select(source, ion.mode, name, formula, mz, rt, med.ic, n.samples)

write.csv(
  cd_compounds,
  here("data/processed/serum_metabolomics/ion_counts/maven_merger/unique compounds and formulas from compound discoverer.csv"),
  row.names = FALSE
)

all_cpds = bind_rows(mav_compounds, cd_compounds) %>%
  arrange(ion.mode, name, formula, source, rt, mz)

write.csv(
  all_cpds,
  here("data/processed/serum_metabolomics/ion_counts/maven_merger/unique maven and cd compounds before consolidation.csv"),
  row.names = FALSE
)

overlap_audit = inner_join(
  cd_compounds %>%
    select(ion.mode, name, formula, mz.cd = mz, rt.cd = rt, med.ic.cd = med.ic, n.samples.cd = n.samples),
  mav_compounds %>%
    select(ion.mode, name, formula, mz.mav = mz, rt.mav = rt, med.ic.mav = med.ic, n.samples.mav = n.samples),
  by = c("ion.mode", "name", "formula")
) %>%
  mutate(
    merge.action = "Maven kept; overlapping CD observations removed"
  ) %>%
  arrange(ion.mode, name)

write.csv(
  overlap_audit,
  here("data/processed/serum_metabolomics/ion_counts/maven_merger/automatic overlap audit between CD and maven.csv"),
  row.names = FALSE
)

mav_data = mav_obs %>%
  transmute(
    tissue,
    raw.file,
    source = "Maven",
    sequence.id,
    sample,
    cohort,
    timepoint,
    ion.mode,
    name,
    formula,
    mz = medmz,
    rt = medrt,
      ion.count
  )

overlap_keys = overlap_audit %>%
  distinct(ion.mode, name, formula)

cd_data_for_merge = cd_data %>%
  anti_join(overlap_keys, by = c("ion.mode", "name", "formula"))

expdata = bind_rows(
  cd_data_for_merge %>% select(tissue, raw.file, source, sequence.id, sample, cohort, timepoint, ion.mode, name, formula, mz, rt, ion.count),
  mav_data
) %>%
  arrange(timepoint, sample, ion.mode, name, source)

write.csv(
  expdata,
  here("data/processed/serum_metabolomics/ion_counts/annotated and processed CD plus MAV obs.csv"),
  row.names = FALSE
)


# ========== 7.0 - Valine and Median Normalization ==========
# --

expdata = read.csv(here("data/processed/serum_metabolomics/ion_counts/annotated and processed CD plus MAV obs.csv"))

expdata2 = expdata %>%
  mutate(sample.id = paste(tissue, sample, timepoint, sep = "_"))

valine = expdata2 %>%
  filter(name == "Valine-C13N15") %>%
  group_by(tissue, sequence.id, sample, sample.id, timepoint, ion.mode) %>%
  summarise(
    val.ic = max(ion.count, na.rm = TRUE),
    .groups = "drop"
  )

expdata3 = expdata2 %>%
  left_join(valine, by = c("tissue", "sequence.id", "sample", "sample.id", "timepoint", "ion.mode"))

missing_valine = expdata3 %>%
  filter(is.na(val.ic)) %>%
  distinct(tissue, sequence.id, sample, sample.id, timepoint, ion.mode)

write.csv(
  missing_valine,
  here("data/processed/serum_metabolomics/ion_counts/maven_merger/samples missing valine after merge.csv"),
  row.names = FALSE
)

expdata4 = expdata3 %>%
  group_by(ion.mode) %>%
  mutate(
    med.val = median(val.ic, na.rm = TRUE),
    val.adj.ic = dplyr::if_else(
      is.na(val.ic) | val.ic <= 0,
      NA_real_,
      (ion.count / val.ic) * med.val
    )
  ) %>%
  ungroup()

expdata5 = expdata4 %>%
  group_by(tissue, name) %>%
  mutate(med.metab = median(val.adj.ic, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    log10.ic.over.med.metab = dplyr::if_else(
      is.na(val.adj.ic) | is.na(med.metab) | val.adj.ic <= 0 | med.metab <= 0,
      NA_real_,
      log10(val.adj.ic / med.metab)
    )
  ) %>%
  group_by(tissue, sample.id) %>%
  mutate(
    med.sample = median(log10.ic.over.med.metab, na.rm = TRUE),
    med.adj.ic = dplyr::if_else(
      is.na(log10.ic.over.med.metab) | is.na(med.sample),
      val.adj.ic,
      (10^(log10.ic.over.med.metab - med.sample)) * med.metab
    )
  ) %>%
  ungroup() %>%
  select(
    tissue, raw.file, source, sequence.id, sample, cohort, timepoint, sample.id, ion.mode,
    name, formula, mz, rt, ion.count, val.ic, med.val, val.adj.ic, med.adj.ic
  )

write.csv(
  expdata5,
  here("data/processed/serum_metabolomics/ion_counts/tidy and normalized serum ion count data.csv"),
  row.names = FALSE
)

expdata_wide = expdata5 %>%
  mutate(
    wide.ic = dplyr::coalesce(med.adj.ic, ion.count)
  ) %>%
  group_by(name, formula, ion.mode, source, mz, rt, sample.id) %>%
  summarise(
    wide.ic = ifelse(
      all(is.na(wide.ic)),
      NA_real_,
      max(wide.ic, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  select(name, formula, ion.mode, source, mz, rt, sample.id, wide.ic) %>%
  pivot_wider(names_from = sample.id, values_from = wide.ic)

write.csv(
  expdata_wide,
  here("data/processed/serum_metabolomics/ion_counts/tidy and normalized serum ion count data wide.csv"),
  row.names = FALSE
)
