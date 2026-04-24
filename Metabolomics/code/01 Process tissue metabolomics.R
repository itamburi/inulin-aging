# ****** Inulin Aging Mouse Tissue (Liver and Muscle)
# ****** Spring 2026 - Ian T, Jessie K, Randa L



# 
here::i_am("code/01 Process tissue metabolomics.R")
source(here::here("code/z-source.R"))
#

# ========== 0.0 - Enviornment ==========
# --

dir.create(here::here("data/processed/ion counts"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("data/processed/ion counts/maven merger"), recursive = TRUE, showWarnings = FALSE)

get_tissue_from_raw_file <- function(raw_file) {
  dplyr::case_when(
    grepl("muscle", raw_file, ignore.case = TRUE) ~ "Muscle",
    grepl("liver", raw_file, ignore.case = TRUE) ~ "Liver",
    TRUE ~ NA_character_
  )
}

get_ion_mode_from_text <- function(x) {
  dplyr::case_when(
    grepl("discover_positive\\.csv$", x, ignore.case = TRUE) ~ "posi",
    grepl("discover_negative\\.csv$", x, ignore.case = TRUE) ~ "nega",
    grepl("_(?:pos)(?:_|$)", x, ignore.case = TRUE, perl = TRUE) ~ "posi",
    grepl("_(?:neg)(?:_|$)", x, ignore.case = TRUE, perl = TRUE) ~ "nega",
    TRUE ~ NA_character_
  )
}

get_tissue_cohort <- function(sequence_id) {
  stringr::str_extract(sequence_id, "^\\d+")
}

get_tissue_sample <- function(sequence_id) {
  sub("_(neg|pos)", "", sequence_id)
}

normalize_tissue_sequence_id <- function(sequence_id) {
  sequence_id %>%
    stringr::str_replace("(_(?:neg|pos))(?:_rerun|_[0-9]{8,})$", "\\1")
}

normalize_tissue_sample_id <- function(sample_id) {
  sample_id %>%
    stringr::str_replace("(?:_rerun|_[0-9]{8,})$", "")
}

get_maven_sample_cols <- function(col_names) {
  col_names[
    grepl(
      "^(?:[124]M\\d+|[124]M_\\d+)_(?:pos|neg)(?:_rerun|_[0-9]{8,})?$",
      col_names,
      ignore.case = TRUE
    ) &
      !grepl("DDA|Blank|STD", col_names, ignore.case = TRUE)
  ]
}

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
library(readxl)

# Loop through each .xlsx file and convert to .csv
CD_files <- list.files(here("data/raw tissue"), pattern = "\\.xlsx$", full.names = TRUE)

for (CD_file in CD_files) {

  file_name = tools::file_path_sans_ext(basename(CD_file))   # Extract the file name without extension

  data = readxl::read_xlsx(CD_file)   # Read .xlsx file
  
  csv_file = file.path(here("data/raw tissue"), paste0(file_name, ".csv"))   # Set the output .csv file path

  write.csv(data, csv_file, row.names = FALSE)  # Write data to .csv file
  cat("Converted", CD_file, "to", csv_file, "\n")
}

# ========== 2.0 - Format Compound Disvoverer Data ==========
# --

# List all CD .csv files in the raw directory
csv_files = list.files(here("data/raw tissue"), pattern = "\\.csv$", full.names = TRUE)

# Loop through each .csv file and read into a list of data frames
csv_data_list = list()
for (csv_file in csv_files) {
  csv_data = read.csv(csv_file)
  csv_data_list[[csv_file]] = csv_data
}


# Function to tidy each data frame. Remove/rename columns, and reshape the structure
tidy_df_list <- function(filename) {

  cd_df = csv_data_list[[filename]]
  file_label = basename(filename)

  measurement_cols = character()
  if (grepl("muscle", file_label, ignore.case = TRUE)) {
    measurement_cols = names(cd_df)[
      grepl("^Area\\.\\.", names(cd_df)) &
        names(cd_df) != "Area..Max.."
    ]
  } else if (grepl("liver", file_label, ignore.case = TRUE)) {
    measurement_cols = names(cd_df)[
      grepl("^Group\\.Area\\.\\.", names(cd_df))
    ]
  }

  if (length(measurement_cols) == 0) {
    stop(paste0("No ion-count columns found for file: ", file_label), call. = FALSE)
  }

  cleaned_sequence_ids = measurement_cols
  if (grepl("muscle", file_label, ignore.case = TRUE)) {
    cleaned_sequence_ids = cleaned_sequence_ids %>%
      str_remove("^Area\\.\\.") %>%
      str_remove("\\.raw\\.\\..*$")
  } else if (grepl("liver", file_label, ignore.case = TRUE)) {
    cleaned_sequence_ids = cleaned_sequence_ids %>%
      str_remove("^Group\\.Area\\.\\.")
  }

  id_cols = intersect(
    c(
      "Name", "Formula", "Calc..MW", "m.z", "RT..min.",
      "Annot..Source..mzCloud.Search", "Annot..Source..MassList.Search",
      "mzCloud.Best.Match.Confidence", "Annot..DeltaMass..ppm."
    ),
    names(cd_df)
  )

  tidy = cd_df %>%
    select(all_of(c(id_cols, measurement_cols))) %>%
    rename(
      mz = m.z,
      calc.mw = Calc..MW,
      rt.min = RT..min.,
      annot.source.mzcloud = Annot..Source..mzCloud.Search,
      annot.source.masslist = Annot..Source..MassList.Search,
      mzcloud.conf = mzCloud.Best.Match.Confidence,
      delta.mass.ppm = Annot..DeltaMass..ppm.
    ) %>%
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
    rename_with(tolower)
  
  print(paste0("Completed ", filename))
  
  return(tidy)
  
}

# for every dataframe within the list of dfs, lapply() the tidy function above
tidy_list = lapply(names(csv_data_list), tidy_df_list)
tidy_df = bind_rows(tidy_list)

# save all features:
write.csv(tidy_df, here("data/processed/ion counts/raw CD observations combined and melted.csv"), row.names = FALSE)

# Apply initial filtering
tidy_df2 = tidy_df %>%
  filter( 
    # (1) remove unannotated names
    is.na(name) == FALSE,
    # (2) remove similar to compounds
    grepl("\\[Similar to:", name) == FALSE,
    # (3) remove blanks
    grepl("blank", sequence.id, ignore.case = TRUE) == FALSE,
    # (4) Remove delta mass outside of +/- 10 ppm
    abs(delta.mass.ppm) <= 10
  )

# save annotated compounds:
write.csv(tidy_df2, here("data/processed/ion counts/filtered CD observations combined and melted.csv"), row.names = FALSE)


# --- decide if we want to take out the chemspider annotations

# ========== 3.0 - Check Annotations Against Reference List, Apply mz Confidence Score Filtering ==========
# -- NOTE FOR FUTURE CLEANUP:
# -- The current reference-list matching can create duplicated annotations downstream for a small number of metabolites
# -- when the reference library contains the same compound name at multiple acceptable RTs.
# -- In that case, `keep` can retain multiple rows for the same observed feature (same name.exp + mz + rt.min),
# -- which then duplicates rows in `cpds3` and inflates the later `inner_join(expdata, cpds3)`.
# -- Current saved data show this for 9-Decenoic acid and L-Lactic acid.

# -- load annotated ion count data and make a df of unique compound annotations
expdata  = read.csv(here("data/processed/ion counts/filtered CD observations combined and melted.csv"), fileEncoding = "UTF-8-BOM") %>%
  mutate(name.exp = tolower(name))

cpds = expdata %>%
  distinct(raw.file, name, name.exp, formula, annot.source.mzcloud, annot.source.masslist, delta.mass.ppm, calc.mw, mz, rt.min, mzcloud.conf) %>%
  select(name, name.exp, everything())
length(unique(cpds$name))

 # -- 1) Check against reference list.
# -- > If a name matches our reference list but the RT is > 3 min away from the reference list, we will flag it for removal
reference_list = read.csv(here("data/libraries/metabolites list_230424_ForCDprocessing.csv"), fileEncoding = "UTF-8-BOM") %>%
  #filter(is.na(RT) == FALSE) %>%
  rename_with(~ tolower(.)) %>%
  select( "hmdb.name", "rt", "formula") %>%
  rename( c("rt.ref" = "rt", "name.ref" = "hmdb.name")) %>%
  mutate( name.ref = trimws(name.ref) )

intersect = inner_join( cpds, reference_list, by = c("name.exp" = "name.ref"), keep = TRUE) %>%
  select("name.exp", "name.ref", "mz", "rt.min", "rt.ref") %>%
  mutate(
    rt.diff = abs(rt.min - rt.ref),
    rt.match = ifelse( rt.diff < 3, TRUE, FALSE)
  )

length(unique(subset(intersect, rt.match == TRUE)$name.exp)) # 792*
length(unique(subset(intersect, rt.match == FALSE)$name.exp)) # 52*

# -- These are the CD metabolties that matched our reference list but the RT was > 3 minutes different. We will remove these! 
discard = intersect %>%
  filter( rt.match == FALSE ) %>%
  group_by(name.exp, name.ref, mz, rt.min, rt.ref, rt.match) %>%
  summarise(
    count = n()
  )
length(unique(discard$name.exp))
# write.csv(discard, here("data/processed/discarded/summary of annotations removed from expdata due to non-matching RTs.csv"), row.names = FALSE )

# -- this is a summary of the CD metabolite-RT combinations that DID match our reference list. We will keep these!
keep = intersect %>%
  filter( rt.match == TRUE ) %>%
  group_by(name.exp, name.ref, mz, rt.min, rt.ref, rt.match) %>%
  summarise(
    count = n()
  )
length(unique(keep$name.exp))

# -- anti_join() will remove the observations in 'discard' from 'expcpds'
cpds2 = anti_join(cpds, discard) %>%
  left_join(., keep[,c("name.exp","mz","rt.min","rt.match")])

# -- mz cloud filtering among the cases that didnt match our ref list
mz_filt = cpds2 %>%
  filter(
    is.na(rt.match), # if rt.match == TRUE we keep regardless of mzcloud.conf. FALSE cases were already removed
    mzcloud.conf < 40
  ) 

# -- anti_join() will remove the observations from 'mz_filt'
cpds3 = anti_join(cpds2, mz_filt)


# -- Finally, Keep only the full expdata with annotations in cpds3
expdata2 = inner_join(expdata, cpds3)
length(unique(expdata2$name))
setdiff(expdata$name, expdata2$name)

write.csv(expdata2, here("data/processed/ion counts/CD data filtered by mismatches to ref and mz cloud conf.csv"),row.names = FALSE)



# ========== 4.0 - Resolve redundant compound entries / select top isomer ==========
# -- Manage repeated names with different RTs (but not yet consolidating by ion mode)

expdata = read.csv(here("data/processed/ion counts/CD data filtered by mismatches to ref and mz cloud conf.csv")) %>%
  mutate(
    cohort = get_tissue_cohort(sequence.id),
    tissue = get_tissue_from_raw_file(raw.file),
    ion.mode = get_ion_mode_from_text(raw.file)
  )

unique(expdata$cohort)
expdata %>% filter(is.na(cohort)) %>% distinct(sequence.id) %>% pull(sequence.id)

# --- How many replicate annotations for same name? ---
n = expdata %>%
  distinct(raw.file, name, mz, formula, rt.min) %>%
  group_by(raw.file, name) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  filter(n > 1)
# .. many metabolites still have multiple annotations per ion mode

# --- (A) Apply a function select_top_isomer() ---
# When multiple features share a compound name, select_top_isomer() selects the preferred peak and labels others as isomers.


# -- Make a list of data frames for each tissue acquisition batch.
# -- Run isomer selection separately within each cohort and raw file to preserve within-batch consistency.
batch_list = expdata %>%
  group_by(raw.file) %>%
  group_split(.)

expdata2 = bind_rows(lapply(batch_list, select_top_isomer3))

# -- data frame with secondary isomers removed
expdata3 = expdata2 %>% filter(grepl("(:| )isomer", name) == FALSE)

# -- consolidate the RT/MZ for the primary cluster and select the best ion count in the cluster for the sample
expdata4 = expdata3 %>%
  group_by(raw.file, sequence.id, name, formula) %>%
  arrange(rt.min, .by_group = TRUE) %>%   # ensures rows are ordered by rt.min
  summarise(
    rt = first(rt.min),           # == representative_rt for the primary cluster
    mz = first(mz),                   # corresponding MZ for the first RT
    ion.count = max(ion.count, na.rm = TRUE)  # Largest annotation PER SAMPLE
  ) %>%
  ungroup()



# --- check the number of occurences of each metabolite Name ---
t = expdata4 %>%
  group_by( name, raw.file, formula, , mz, rt ) %>%
  summarise(median.ic = median(ion.count))
n = t %>% group_by(name) %>% summarise(n=n()) %>% arrange(desc(n))
# ...Should Max 4 for everything, depending on n ion modes for the metabolite, and n tissues, if we did everything correctly

# --- Evaluate RT drift between the tissues (per ion mode) ---
expdata5 = expdata4 %>%
  mutate(
    tissue = case_when(
      grepl("muscle", raw.file, ignore.case = TRUE) ~ "Muscle",
      grepl("liver", raw.file, ignore.case = TRUE) ~ "Liver",
      TRUE ~ NA_character_
    ),
    ion.mode = get_ion_mode_from_text(raw.file),
  )


rt_drift_tissues = expdata5 %>%
  group_by( name, formula, ion.mode) %>%
  summarise( rt_drift = diff(range(rt, na.rm = TRUE)) )

remove = rt_drift_tissues %>% filter(rt_drift > 1)
# these metabolites (per ion mode) have RT difference > 1 minute between the tissues and are probably unreliable
# any rt drift < 1 min bw tissues is acceptable


expdata6 = expdata5 %>%
  anti_join(., remove) %>%
  # let's consolidate the rt drift < 1 min here
  group_by(ion.mode, name, formula) %>%
  arrange(rt, .by_group = TRUE) %>%
  mutate(
    rt = first(rt),
    mz = first(mz)
  ) %>%
  ungroup() %>%
  select(raw.file, tissue, ion.mode, sequence.id, name, formula, rt, mz, ion.count)


n = expdata6 %>%
  distinct(name, tissue, ion.mode, formula, mz, rt) %>%
  count(name, sort = TRUE)

dup_check = expdata6 %>%
  distinct(name, tissue, ion.mode, formula, mz, rt) %>%
  count(name, tissue, ion.mode, sort = TRUE) %>%
  filter(n > 1)



# save dereplicated data
write.csv(expdata6, here("data/processed/ion counts/annotated and dereplicated CD obs without repeat isomers.csv"), row.names = FALSE)
write.csv(remove, here("data/processed/ion counts/names removed due to tissue RT drift or nonsense name.csv"), row.names = FALSE)

# ========== 5.0 - Consolidate Negative and Positive Ion Modes ==========
# -- Among metabolites detected in +&- modes, select the annotation that has the highest median ion count

# -- read in dereplicated data
expdata = read.csv(here("data/processed/ion counts/annotated and dereplicated CD obs without repeat isomers.csv"), fileEncoding = "UTF-8-BOM") %>%
  mutate(
    sample = get_tissue_sample(sequence.id)
  ) 

# -- confirm the n annotations per mode as in previous step
n_modes = expdata %>%
  distinct( name, ion.mode, formula, mz, rt ) %>%
  group_by(name, ion.mode) %>%
  mutate(n.per.mode = n()) %>%
  ungroup() %>%
  arrange(desc(n.per.mode), name, ion.mode, rt)
# anything with n.per.mode == 2 had an acceptable RT drift bw the tissues.
# anything greater than 2 would be an error that we would go back and address in prior step

# -- Determine the median of every metabolite 'Name' per ion mode
ic_medians = expdata %>%
  group_by(name, ion.mode, formula) %>%
  summarise(median.ic = median(ion.count)) %>%
  ungroup() %>%
  group_by(name) %>%
  mutate(n = n())

# -- Extract the annotation with the larger median across all the samples
# ** we do not group by tissue here; the preferred ion mode should match across tissues
top_mode = ic_medians %>% group_by(name) %>% top_n(1, median.ic) 
length(unique(top_mode$name))
nrow(top_mode)
# rows == # names when a single top mode is selected for each metabolite

# -- Select from expdata only the top median annotations
expdata2 = inner_join(expdata, top_mode %>% select(-median.ic, -n) )

# -- check that we only have 1 mode per metabolite
n_modes2 = expdata2 %>%
  distinct( name, ion.mode, formula) %>%
  group_by(name) %>%
  summarise(n.obs = n())


names(expdata2)
expdata3 = expdata2 %>%
  select(raw.file, tissue, ion.mode, sequence.id, sample, name, formula, rt, mz, ion.count)

# -- save the data consolidated by ion mode
write.csv(expdata3, here("data/processed/ion counts/annotated and processed CD obs consolidated by top ion mode.csv"), row.names = FALSE)
#

# ========== 6.0 - Merge and Consolidate CD and Maven Observations ==========
# --

mav_dir = "data/raw tissue/maven export"
mav_files = list.files(mav_dir, pattern = "\\.csv$", full.names = TRUE)

mav_obs_list = lapply(mav_files, function(csv_file) {
  csv_data = read.csv(csv_file, check.names = FALSE)
  sample_cols = get_maven_sample_cols(names(csv_data))

  csv_data %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sequence.id",
      values_to = "ion.count"
    ) %>%
    transmute(
      raw.file = basename(csv_file),
      tissue = get_tissue_from_raw_file(raw.file),
      sequence.id = normalize_tissue_sequence_id(sequence.id),
      ion.mode = get_ion_mode_from_text(sequence.id),
      sample = normalize_tissue_sample_id(get_tissue_sample(sequence.id)),
      isotope.label = isotopeLabel,
      name = normalize_maven_compound_name(compound, isotopeLabel),
      formula = formula,
      medrt = as.numeric(medRt),
      medmz = as.numeric(medMz),
      ion.count = as.numeric(ion.count)
    ) %>%
    filter(
      is.na(sequence.id) == FALSE,
      is.na(ion.mode) == FALSE,
      is.na(ion.count) == FALSE
    ) %>%
    distinct()
})

mav_obs = bind_rows(mav_obs_list) %>%
  group_by(tissue, raw.file, sequence.id, ion.mode, sample, name, formula) %>%
  summarise(
    medrt = median(medrt, na.rm = TRUE),
    medmz = median(medmz, na.rm = TRUE),
    ion.count = max(ion.count, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(mav_obs, here("data/processed/ion counts/maven merger/all maven observations.csv"), row.names = FALSE)

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

write.csv(mav_compounds, here("data/processed/ion counts/maven merger/unique compounds and formulas from maven.csv"), row.names = FALSE)

cd_data = read.csv(
  here("data/processed/ion counts/annotated and processed CD obs consolidated by top ion mode.csv"),
  fileEncoding = "UTF-8-BOM"
) %>%
  mutate(
    source = "CD",
    sequence.id = normalize_tissue_sequence_id(sequence.id),
    sample = normalize_tissue_sample_id(sample),
    tissue = get_tissue_from_raw_file(raw.file)
  ) %>%
  group_by(tissue, raw.file, sequence.id, ion.mode, sample, name, formula, mz, rt, source) %>%
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

write.csv(cd_compounds, here("data/processed/ion counts/maven merger/unique compounds and formulas from compound discoverer.csv"), row.names = FALSE)

all_cpds = bind_rows(mav_compounds, cd_compounds) %>%
  arrange(ion.mode, name, formula, source, rt, mz)

write.csv(all_cpds, here("data/processed/ion counts/maven merger/unique maven and cd compounds before consolidation.csv"), row.names = FALSE)

overlap_audit = inner_join(
  cd_compounds %>%
    select(ion.mode, name, formula, mz.cd = mz, rt.cd = rt, med.ic.cd = med.ic, n.samples.cd = n.samples),
  mav_compounds %>%
    select(ion.mode, name, formula, mz.mav = mz, rt.mav = rt, med.ic.mav = med.ic, n.samples.mav = n.samples),
  by = c("ion.mode", "name", "formula")
) %>%
  arrange(ion.mode, name)

write.csv(overlap_audit, here("data/processed/ion counts/maven merger/automatic overlap audit between CD and maven.csv"), row.names = FALSE)

mav_data = mav_obs %>%
  transmute(
    tissue,
    raw.file,
    source = "Maven",
    sequence.id,
    sample,
    ion.mode,
    name,
    formula,
    mz = medmz,
    rt = medrt,
    ion.count
  )

expdata = bind_rows(
  cd_data %>% select(tissue, raw.file, source, sequence.id, sample, ion.mode, name, formula, mz, rt, ion.count),
  mav_data
) %>%
  arrange(tissue, sample, ion.mode, name, source)

write.csv(expdata, here("data/processed/ion counts/annotated and processed CD plus MAV obs.csv"), row.names = FALSE)

# ========== 7.0 - Valine and Median Normalization ==========
# --
expdata = read.csv(here("data/processed/ion counts/annotated and processed CD plus MAV obs.csv"))

expdata2 = expdata %>%
  mutate(
    sample.id = paste(tissue, sample, sep = "_")
  )

valine = expdata2 %>%
  filter(name == "Valine-C13N15") %>%
  group_by(tissue, sequence.id, sample, sample.id, ion.mode) %>%
  summarise(
    val.ic = max(ion.count, na.rm = TRUE),
    .groups = "drop"
  )

expdata3 = expdata2 %>%
  left_join(valine, by = c("tissue", "sequence.id", "sample", "sample.id", "ion.mode"))

missing_valine = expdata3 %>%
  filter(is.na(val.ic)) %>%
  distinct(tissue, sequence.id, sample, ion.mode)

# --- Valine Normalization
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

# --- Median Normalization
expdata5 = expdata4 %>%
  group_by(tissue, name) %>%
  mutate(
    med.metab = median(val.adj.ic, na.rm = TRUE)
  ) %>%
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
    tissue, raw.file, source, sequence.id, sample, sample.id, ion.mode,
    name, formula, mz, rt, ion.count, val.ic, med.val, val.adj.ic, med.adj.ic
  )

write.csv(expdata5, here("data/processed/ion counts/tidy and normalized ion count data.csv"), row.names = FALSE)

expdata_wide = expdata5 %>%
  group_by(name, formula, ion.mode, source, mz, rt, sample.id) %>%
  summarise(
    med.adj.ic = ifelse(
      all(is.na(med.adj.ic)),
      NA_real_,
      max(med.adj.ic, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  select(name, formula, ion.mode, source, mz, rt, sample.id, med.adj.ic) %>%
  pivot_wider(names_from = sample.id, values_from = med.adj.ic)

write.csv(expdata_wide, here("data/processed/ion counts/tidy and normalized ion count data wide.csv"), row.names = FALSE)

