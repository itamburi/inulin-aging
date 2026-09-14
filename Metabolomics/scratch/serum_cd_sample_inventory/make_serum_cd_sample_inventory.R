# Build a sample/run inventory from serum Compound Discoverer exports.

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(script_file_arg) > 0) {
  normalizePath(sub("^--file=", "", script_file_arg[[1]]))
} else {
  normalizePath("metabolomics/scratch/serum_cd_sample_inventory/make_serum_cd_sample_inventory.R")
}

scratch_dir <- dirname(script_path)
metabolomics_dir <- normalizePath(file.path(scratch_dir, "..", ".."))

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("The readxl package is required to read Compound Discoverer xlsx files.", call. = FALSE)
}

inputs <- data.frame(
  dataset = c("serum", "serum", "serum_longitudinal", "serum_longitudinal"),
  compound_discoverer_file = c(
    "data/raw/serum/compound_discoverer/0100_Inulin Aging_serum_pos_original.xlsx",
    "data/raw/serum/compound_discoverer/0100_Inulin_Aging_serum_neg_original.xlsx",
    "data/raw/serum_longitudinal/compound discoverer/0300 Longitudinal Inulin Aging Area Pos.xlsx",
    "data/raw/serum_longitudinal/compound discoverer/0300 Longitudinal Inulin Aging Area Neg.xlsx"
  ),
  stringsAsFactors = FALSE
)

died_samples <- c(
  "3M3", "3M5", "3M9", "3M10", "3M12", "3M19", "3M20", "3M22",
  "3M23", "3M24", "3M25", "3M29", "3M30", "3M32", "3M35", "3M36"
)

normalize_sequence_id <- function(sequence_id) {
  sequence_id <- sub("\\.raw.*$", "", sequence_id, ignore.case = TRUE)
  sequence_id <- sub("_+$", "", sequence_id)
  sequence_id <- sub(
    "(_(?:pos|neg))(?:_(?:rerun|resume)[0-9]*|_[0-9]{8,})+$",
    "\\1",
    sequence_id,
    perl = TRUE,
    ignore.case = TRUE
  )
  sub("_+$", "", sequence_id)
}

get_sample <- function(sequence_id) {
  sample <- regmatches(sequence_id, regexpr("^[1-5]M[0-9]+", sequence_id, perl = TRUE))
  ifelse(length(sample) == 0, NA_character_, sample)
}

get_timepoint <- function(sequence_id, dataset) {
  if (dataset == "serum") {
    if (grepl("_pv_", sequence_id, ignore.case = TRUE)) return("pv")
    if (grepl("_sys_sac_", sequence_id, ignore.case = TRUE)) return("sys_sac")
    if (grepl("_tp4_", sequence_id, ignore.case = TRUE)) return("tp4")
    return(NA_character_)
  }

  sample <- get_sample(sequence_id)
  if (grepl("_mid_", sequence_id, ignore.case = TRUE)) return("mid")
  if (grepl("_old_", sequence_id, ignore.case = TRUE)) return("old")
  if (grepl("^5M[0-9]+_(?:pos|neg)$", sequence_id, ignore.case = TRUE, perl = TRUE)) return("young")
  if (!is.na(sample) && grepl("^5M", sample, ignore.case = TRUE)) return("young")
  NA_character_
}

get_ion_mode <- function(text) {
  if (grepl("positive|\\bpos\\b|(?:^|_)pos(?:_|$)", text, ignore.case = TRUE, perl = TRUE)) return("posi")
  if (grepl("negative|\\bneg\\b|(?:^|_)neg(?:2)?(?:_|$)", text, ignore.case = TRUE, perl = TRUE)) return("nega")
  NA_character_
}

get_collection_protocol <- function(dataset, timepoint) {
  protocol <- rep(NA_character_, length(timepoint))

  protocol[timepoint %in% c("pv", "sys_sac")] <- "sacrifice"
  protocol[dataset == "serum" & timepoint == "tp4"] <- "GTT"
  protocol[dataset == "serum_longitudinal" & timepoint %in% c("mid", "old")] <- "DNL"
  protocol[dataset == "serum_longitudinal" & timepoint == "young"] <- "sacrifice"

  protocol
}

get_timepoint_months <- function(dataset, timepoint) {
  months <- rep(NA_character_, length(timepoint))

  months[dataset == "serum" & timepoint %in% c("pv", "sys_sac")] <- "24-26 months"
  months[dataset == "serum" & timepoint == "tp4"] <- "13 months"
  months[dataset == "serum_longitudinal" & timepoint == "young"] <- "7-14 months"
  months[dataset == "serum_longitudinal" & timepoint == "mid"] <- "14 months"
  months[dataset == "serum_longitudinal" & timepoint == "old"] <- "22 months"

  months
}

get_diet <- function(cohort, animal_id) {
  diet <- rep(NA_character_, length(cohort))
  cohort <- as.character(cohort)

  diet[cohort == "1" & animal_id %in% c(1, 2, 4, 10, 11)] <- "CD"
  diet[cohort == "1" & animal_id %in% c(5, 6, 7, 9)] <- "ID"

  diet[cohort == "2" & animal_id %in% c(16, 17, 20)] <- "CD"
  diet[cohort == "2" & animal_id %in% c(22, 23, 26, 31, 33)] <- "ID"

  diet[cohort == "3" & (
    animal_id %in% 1:10 | animal_id %in% 21:30
  )] <- "CD"
  diet[cohort == "3" & (
    animal_id %in% 11:20 | animal_id %in% 31:40
  )] <- "ID"

  diet[cohort %in% c("4", "5")] <- "CD"

  diet
}

get_gender <- function(cohort, animal_id) {
  gender <- rep(NA_character_, length(cohort))
  cohort <- as.character(cohort)

  gender[cohort %in% c("1", "2", "4")] <- "Female"
  gender[cohort == "3" & animal_id %in% 1:20] <- "Female"
  gender[cohort == "3" & animal_id %in% 21:40] <- "Male"
  gender[cohort == "5"] <- "Male"

  gender
}

parse_cd_header <- function(header) {
  if (grepl("^Area: ", header)) {
    label <- sub("^Area: ", "", header)
    injection_position <- if (grepl("\\([A-Za-z]+[0-9]+\\)$", label)) {
      sub("^.*\\(([A-Za-z]+[0-9]+)\\)$", "\\1", label)
    } else {
      NA_character_
    }
    raw_file <- trimws(sub("\\s*\\([A-Za-z]+[0-9]+\\)$", "", label))

    return(data.frame(
      cd_column_type = "Area",
      cd_column = header,
      raw_file = raw_file,
      raw_file_source = "area_header",
      injection_position = injection_position,
      stringsAsFactors = FALSE
    ))
  }

  label <- sub("^Group Area: ", "", header)
  data.frame(
    cd_column_type = "Group Area",
    cd_column = header,
    raw_file = paste0(label, ".raw"),
    raw_file_source = "group_area_inferred",
    injection_position = NA_character_,
    stringsAsFactors = FALSE
  )
}

extract_sample_inventory <- function(dataset, compound_discoverer_file) {
  input_path <- file.path(metabolomics_dir, compound_discoverer_file)
  if (!file.exists(input_path)) {
    stop("Missing input file: ", input_path, call. = FALSE)
  }

  headers <- names(readxl::read_xlsx(input_path, n_max = 0, .name_repair = "minimal"))
  area_cols <- headers[grepl("^Area: ", headers) & headers != "Area (Max.)"]
  group_area_cols <- headers[grepl("^Group Area: ", headers)]
  sample_cols <- if (length(area_cols) > 0) area_cols else group_area_cols
  sample_cols <- sample_cols[!grepl("blank|std", sample_cols, ignore.case = TRUE)]

  if (length(sample_cols) == 0) {
    stop("No sample Area or Group Area columns found in: ", input_path, call. = FALSE)
  }

  parsed <- do.call(rbind, lapply(sample_cols, parse_cd_header))
  parsed$dataset <- dataset
  parsed$compound_discoverer_file <- basename(compound_discoverer_file)
  parsed$compound_discoverer_path <- compound_discoverer_file
  parsed$source_column_index <- match(parsed$cd_column, headers)
  parsed$sequence_id <- sub("\\.raw.*$", "", parsed$raw_file, ignore.case = TRUE)
  parsed$normalized_sequence_id <- normalize_sequence_id(parsed$sequence_id)
  parsed$sample <- vapply(parsed$normalized_sequence_id, get_sample, character(1))
  parsed$cohort <- sub("^([1-5]).*$", "\\1", parsed$sample)
  parsed$cohort[is.na(parsed$sample)] <- NA_character_
  parsed$age <- ifelse(parsed$cohort %in% c("1", "2", "3", "4"), "Old", "Young")
  parsed$age[is.na(parsed$cohort)] <- NA_character_
  animal_id <- as.integer(sub("^[1-5]M", "", parsed$sample))
  parsed$diet <- get_diet(parsed$cohort, animal_id)
  parsed$gender <- get_gender(parsed$cohort, animal_id)
  parsed$died <- ifelse(parsed$sample %in% died_samples, "died", "")
  parsed$timepoint <- vapply(
    parsed$normalized_sequence_id,
    get_timepoint,
    character(1),
    dataset = dataset
  )
  parsed$sample_type <- parsed$timepoint
  parsed$blood_type <- ifelse(parsed$timepoint == "pv", "portal blood", "tail blood")
  parsed$collection_protocol <- get_collection_protocol(parsed$dataset, parsed$timepoint)
  parsed$timepoint_months <- get_timepoint_months(parsed$dataset, parsed$timepoint)
  parsed$ion_mode <- vapply(parsed$normalized_sequence_id, get_ion_mode, character(1))

  parsed[
    c(
      "dataset", "compound_discoverer_file", "compound_discoverer_path",
      "ion_mode", "source_column_index", "cd_column_type", "cd_column", "raw_file",
      "raw_file_source", "injection_position", "sequence_id",
      "normalized_sequence_id", "sample", "cohort", "age", "diet", "gender",
      "died", "sample_type", "timepoint_months", "blood_type",
      "collection_protocol"
    )
  ]
}

inventory <- do.call(
  rbind,
  Map(extract_sample_inventory, inputs$dataset, inputs$compound_discoverer_file)
)

inventory <- inventory[order(
  inventory$dataset,
  inventory$compound_discoverer_file,
  inventory$source_column_index
), ]
rownames(inventory) <- NULL

output_file <- file.path(scratch_dir, "serum_and_serum_longitudinal_cd_sample_inventory.csv")
write.csv(inventory, output_file, row.names = FALSE, na = "")

dataset_order <- c("serum", "serum_longitudinal")
timepoint_order <- c("young", "tp4", "mid", "old", "pv", "sys_sac")

collapse_ordered_unique <- function(x, order = NULL) {
  x <- unique(x[!is.na(x) & x != ""])
  if (!is.null(order)) {
    x <- c(order[order %in% x], sort(setdiff(x, order)))
  } else {
    x <- sort(x)
  }
  paste(x, collapse = ";")
}

sample_summary <- data.frame(
  sample = sort(unique(inventory$sample)),
  stringsAsFactors = FALSE
)
sample_summary$cohort <- sub("^([1-5]).*$", "\\1", sample_summary$sample)
sample_summary$animal_id <- as.integer(sub("^[1-5]M", "", sample_summary$sample))

sample_summary$available_timepoints <- vapply(
  sample_summary$sample,
  function(current_sample) {
    collapse_ordered_unique(
      inventory$sample_type[inventory$sample == current_sample],
      timepoint_order
    )
  },
  character(1)
)

sample_summary$dataset <- vapply(
  sample_summary$sample,
  function(current_sample) {
    collapse_ordered_unique(
      inventory$dataset[inventory$sample == current_sample],
      dataset_order
    )
  },
  character(1)
)

sample_summary <- sample_summary[order(sample_summary$cohort, sample_summary$animal_id), ]
sample_summary <- sample_summary[c("sample", "cohort", "available_timepoints", "dataset")]
rownames(sample_summary) <- NULL

sample_summary_file <- file.path(scratch_dir, "serum_and_serum_longitudinal_cd_samples_by_sample.csv")
write.csv(sample_summary, sample_summary_file, row.names = FALSE, na = "")

summary_table <- aggregate(
  raw_file ~ dataset + compound_discoverer_file + ion_mode + raw_file_source,
  inventory,
  length
)
names(summary_table)[names(summary_table) == "raw_file"] <- "n_rows"

message("Wrote ", nrow(inventory), " sample rows to: ", output_file)
message("Wrote ", nrow(sample_summary), " consolidated sample rows to: ", sample_summary_file)
print(summary_table)
