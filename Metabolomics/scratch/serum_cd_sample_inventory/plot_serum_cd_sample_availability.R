# Plot sample availability across serum Compound Discoverer datasets and sample types.

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(script_file_arg) > 0) {
  normalizePath(sub("^--file=", "", script_file_arg[[1]]))
} else {
  normalizePath("metabolomics/scratch/serum_cd_sample_inventory/plot_serum_cd_sample_availability.R")
}

scratch_dir <- dirname(script_path)
inventory_file <- file.path(scratch_dir, "serum_and_serum_longitudinal_cd_sample_inventory.csv")

required_packages <- c("dplyr", "tidyr", "ggplot2", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)

dataset_order <- c("sample_metadata", "serum_longitudinal", "serum")
sample_type_order <- c("age", "diet", "gender", "died", "young", "tp4", "mid", "old", "pv", "sys_sac")
metadata_column_order <- c("age", "diet", "gender", "died")

dataset_sample_types <- tibble::tribble(
  ~dataset, ~sample_type, ~timepoint_months,
  "serum_longitudinal", "young", "7-14 months",
  "serum", "tp4", "13 months",
  "serum_longitudinal", "mid", "14 months",
  "serum_longitudinal", "old", "22 months",
  "serum", "pv", "24-26 months",
  "serum", "sys_sac", "24-26 months"
)

sample_type_labels <- c(
  age = "Age",
  diet = "Diet",
  gender = "Gender",
  died = "Died",
  young = "young\nSAC\n(7-14 months)",
  tp4 = "tp4\nGTT\n(13 months)",
  mid = "mid\nDNL\n(14 months)",
  old = "old\nDNL\n(22 months)",
  pv = "pv\nSAC\n(24-26 months)",
  sys_sac = "sys_sac\nSAC\n(24-26 months)"
)

format_ion_modes <- function(ion_modes) {
  ordered_modes <- c("posi", "nega")
  present_modes <- ordered_modes[ordered_modes %in% ion_modes]
  paste(dplyr::recode(present_modes, posi = "P", nega = "N"), collapse = "/")
}

inventory <- read.csv(inventory_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample_type" %in% names(inventory) && "timepoint" %in% names(inventory)) {
  inventory$sample_type <- inventory$timepoint
}

cohort_gap <- 1.2

sample_positions <- inventory %>%
  distinct(sample, cohort) %>%
  mutate(
    cohort_numeric = as.integer(cohort),
    animal_id = as.integer(sub("^[1-5]M", "", sample))
  ) %>%
  arrange(cohort_numeric, animal_id) %>%
  mutate(
    cohort_rank = dense_rank(cohort_numeric),
    y_pos = row_number() + (cohort_rank - 1) * cohort_gap
  )

sample_order <- sample_positions$sample

sample_lookup <- inventory %>%
  distinct(sample, cohort, age, diet, gender, died) %>%
  left_join(sample_positions %>% select(sample, y_pos), by = "sample") %>%
  mutate(
    age = ifelse(is.na(age) | age == "", "Unknown", age),
    diet = ifelse(is.na(diet) | diet == "", "Unknown", diet),
    gender = ifelse(is.na(gender) | gender == "", "Unknown", gender),
    died = ifelse(is.na(died), "", died)
  )

availability_observed <- inventory %>%
  group_by(sample, dataset, sample_type, timepoint_months) %>%
  summarise(
    ion_modes = format_ion_modes(unique(ion_mode)),
    n_raw_files = n_distinct(raw_file),
    raw_files = paste(sort(unique(raw_file)), collapse = ";"),
    .groups = "drop"
  )

availability_grid <- tidyr::crossing(
  sample = sample_order,
  dataset_sample_types
) %>%
  left_join(sample_lookup, by = "sample") %>%
  left_join(
    availability_observed,
    by = c("sample", "dataset", "sample_type", "timepoint_months")
  ) %>%
  mutate(
    available = !is.na(ion_modes),
    ion_modes = ifelse(is.na(ion_modes), "", ion_modes),
    n_raw_files = ifelse(is.na(n_raw_files), 0L, n_raw_files),
    raw_files = ifelse(is.na(raw_files), "", raw_files),
    dataset = factor(dataset, levels = dataset_order),
    sample_type = factor(sample_type, levels = sample_type_order),
    availability_fill = case_when(
      available & dataset == "serum" ~ "serum",
      available & dataset == "serum_longitudinal" ~ "serum_longitudinal",
      TRUE ~ "not_available"
    )
  ) %>%
  arrange(
    match(as.character(sample), sample_order),
    match(as.character(sample_type), sample_type_order),
    match(as.character(dataset), dataset_order)
  )

metadata_grid <- sample_lookup %>%
  select(sample, cohort, y_pos, age, diet, gender, died) %>%
  pivot_longer(
    cols = all_of(metadata_column_order),
    names_to = "sample_type",
    values_to = "metadata_value"
  ) %>%
  mutate(
    dataset = factor("sample_metadata", levels = dataset_order),
    sample_type = factor(sample_type, levels = sample_type_order),
    timepoint_months = "",
    available = TRUE,
    ion_modes = metadata_value,
    n_raw_files = NA_integer_,
    raw_files = "",
    availability_fill = case_when(
      sample_type == "age" & metadata_value == "Old" ~ "age_old",
      sample_type == "age" & metadata_value == "Young" ~ "age_young",
      sample_type == "diet" & metadata_value == "CD" ~ "diet_cd",
      sample_type == "diet" & metadata_value == "ID" ~ "diet_id",
      sample_type == "gender" & metadata_value == "Female" ~ "gender_female",
      sample_type == "gender" & metadata_value == "Male" ~ "gender_male",
      sample_type == "died" & metadata_value == "died" ~ "died",
      sample_type == "died" ~ "not_died",
      TRUE ~ "metadata_unknown"
    )
  ) %>%
  arrange(
    match(as.character(sample), sample_order),
    match(as.character(sample_type), sample_type_order)
  )

plot_grid <- bind_rows(
  metadata_grid %>%
    select(sample, cohort, y_pos, dataset, sample_type, timepoint_months, available, ion_modes, n_raw_files, raw_files, availability_fill),
  availability_grid %>%
    select(sample, cohort, y_pos, dataset, sample_type, timepoint_months, available, ion_modes, n_raw_files, raw_files, availability_fill)
)

availability_export <- availability_grid %>%
  select(
    sample, cohort, age, diet, gender, died, dataset, sample_type, timepoint_months,
    available, ion_modes, n_raw_files, raw_files
  )

plot_data_file <- file.path(scratch_dir, "serum_cd_sample_dataset_timepoint_availability.csv")
write.csv(availability_export, plot_data_file, row.names = FALSE, na = "")

availability_plot <- ggplot(
  plot_grid,
  aes(x = sample_type, y = y_pos, fill = availability_fill)
) +
  geom_tile(color = "white", linewidth = 0.2, width = 1, height = 0.9) +
  geom_text(
    aes(label = ion_modes),
    color = "#111827",
    size = 2.15,
    fontface = "bold"
  ) +
  scale_x_discrete(
    labels = sample_type_labels,
    expand = expansion(add = 0.5),
    drop = TRUE
  ) +
  scale_y_reverse(
    breaks = sample_positions$y_pos,
    labels = sample_positions$sample,
    expand = expansion(add = 0.4)
  ) +
  scale_fill_manual(
    values = c(
      serum = "#5B8FA8",
      serum_longitudinal = "#6FA86F",
      not_available = "#ECEFF1",
      not_died = "#ECEFF1",
      diet_cd = "#7E6B4D",
      diet_id = "#D08A3C",
      age_old = "#8875B5",
      age_young = "#D7BF46",
      gender_female = "#C76E8C",
      gender_male = "#4A6FA5",
      died = "#C14953",
      metadata_unknown = "#AEB6C2"
    ),
    breaks = c(
      "age_old", "age_young", "diet_cd", "diet_id", "gender_female",
      "gender_male", "died", "serum_longitudinal", "serum", "not_available"
    ),
    labels = c(
      "Age: Old", "Age: Young", "Diet: CD", "Diet: ID", "Gender: Female",
      "Gender: Male", "Died", "Serum longitudinal present", "Serum original present",
      "Not present"
    ),
    drop = FALSE
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  labs(
    title = "Serum Compound Discoverer Sample Availability",
    subtitle = "Rows are samples sorted by cohort and animal ID; blood sample columns are ordered by collection age.",
    caption = "Metadata comes from the curated sample inventory. Blood cells show ion-mode raw files: P/N, P, or N.",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, color = "#111827"),
    axis.text.y = element_text(color = "#111827", size = 6.5),
    strip.text = element_text(face = "bold", color = "#111827", size = 10),
    plot.title = element_text(face = "bold", color = "#111827", size = 14),
    plot.subtitle = element_text(color = "#374151", size = 9),
    plot.caption = element_text(color = "#4B5563", size = 8, hjust = 0),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    panel.spacing.x = unit(0, "lines")
  )

png_file <- file.path(scratch_dir, "serum_cd_sample_availability_heatmap.png")
pdf_file <- file.path(scratch_dir, "serum_cd_sample_availability_heatmap.pdf")

ggsave(png_file, availability_plot, width = 6.8, height = 14.6, dpi = 300)
ggsave(pdf_file, availability_plot, width = 6.8, height = 14.6)

message("Wrote plot data to: ", plot_data_file)
message("Wrote PNG to: ", png_file)
message("Wrote PDF to: ", pdf_file)
