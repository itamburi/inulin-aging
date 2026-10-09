# ****** Inulin Aging Mouse Serum Hormone-Metabolite Bicor Figure
# ****** Summer 2026 - Ian T, Jessie K, Randa L
#
# Figure focus:
# - Old mice only.
# - Serum TP4 metabolites matched to GTT hormones.
# - Serum systemic sacrifice metabolites matched to SAC hormones.
# - Bicor results are summarized for CD, ID, and combined CD/ID groups.
#
# The heatmap uses a conservative candidate set: combined CD/ID FDR < 0.05,
# same CD and ID direction, nominal p < 0.05 in both CD and ID, and nominal
# p < 0.05 in both single-state checks (TP4/GTT and sys_sac/SAC), all in the
# same direction as the combined estimate.


here::i_am("code/integrative/04_plot_serum_hormone_metabolite_bicor_figure.R")


# ========== 0.0 - Environment ==========
# --
required_packages = c(
  "dplyr", "tidyr", "readr", "stringr", "ggplot2", "patchwork", "scales"
)
missing_packages = required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(patchwork)
library(scales)

output_dir = here::here("data/processed/integrative/serum_hormone_metabolite_correlations")
plot_dir = here::here("plots/integrative/serum_hormone_metabolite_correlations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Cairo PDF needs a writable font cache when this script runs in a sandboxed shell.
Sys.setenv(XDG_CACHE_HOME = tempdir())

hormone_columns = c(
  "C-Peptide 2", "Ghrelin", "GIP", "GLP-1", "Glucagon", "Insulin",
  "Leptin", "Secretin", "PP", "PYY", "Resistin", "Amylin"
)

diet_group_labels = c(
  CD = "CD",
  ID = "ID",
  CD_ID_combined = "CD+ID"
)

correlation_column_order = c(
  "CD pooled states",
  "ID pooled states",
  "CD+ID pooled states",
  "TP4/GTT CD+ID",
  "sys_sac/SAC CD+ID"
)

correlation_fill_limits = c(-1, 1)


# ========== 1.0 - Read bicor result tables ==========
# --
add_pair_id = function(data) {
  data %>%
    mutate(
      metabolite_mz_key = round(metabolite_mz, 5),
      metabolite_rt_key = round(metabolite_rt, 3),
      pair_id = paste(
        metabolite_source,
        metabolite_ion_mode,
        metabolite,
        metabolite_formula,
        metabolite_mz_key,
        metabolite_rt_key,
        hormone,
        sep = "||"
      )
    )
}

combined_correlation_results = readr::read_csv(
  file.path(
    output_dir,
    "old serum tp4 and sys_sac metabolite vs matched gtt and sac hormone correlations.csv"
  ),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  filter(correlation_method == "bicor") %>%
  add_pair_id()

tp4_gtt_correlation_results = readr::read_csv(
  file.path(output_dir, "old serum tp4 metabolite vs gtt hormone correlations.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  filter(
    correlation_method == "bicor",
    diet_correlation_group == "CD_ID_combined"
  ) %>%
  add_pair_id() %>%
  select(
    pair_id,
    tp4_gtt_correlation = correlation,
    tp4_gtt_p_value = p_value,
    tp4_gtt_p_adj_bh = p_adj_bh,
    tp4_gtt_n = n
  )

sys_sac_correlation_results = readr::read_csv(
  file.path(output_dir, "old serum sys_sac metabolite vs sac hormone correlations.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  filter(
    correlation_method == "bicor",
    diet_correlation_group == "CD_ID_combined"
  ) %>%
  add_pair_id() %>%
  select(
    pair_id,
    sys_sac_correlation = correlation,
    sys_sac_p_value = p_value,
    sys_sac_p_adj_bh = p_adj_bh,
    sys_sac_n = n
  )


# ========== 2.0 - Summarise full FDR landscape ==========
# --
significant_count_summary = combined_correlation_results %>%
  filter(!is.na(p_adj_bh), p_adj_bh < 0.05) %>%
  mutate(
    diet_correlation_group = factor(
      diet_correlation_group,
      levels = names(diet_group_labels),
      labels = diet_group_labels
    ),
    correlation_direction = if_else(correlation > 0, "positive", "negative")
  ) %>%
  count(diet_correlation_group, hormone, correlation_direction, name = "n_hits")

hormone_order = significant_count_summary %>%
  group_by(hormone) %>%
  summarise(total_hits = sum(n_hits), .groups = "drop") %>%
  arrange(total_hits, hormone) %>%
  pull(hormone)

count_plot_data = significant_count_summary %>%
  mutate(
    hormone = factor(hormone, levels = hormone_order),
    signed_hits = if_else(correlation_direction == "positive", n_hits, -n_hits)
  )

readr::write_csv(
  significant_count_summary %>%
    arrange(diet_correlation_group, hormone, correlation_direction),
  file.path(output_dir, "old serum tp4 sys_sac bicor figure significant count summary.csv"),
  na = "N/A"
)


# ========== 3.0 - Select conservative heatmap pairs ==========
# --
pair_level_correlation_summary = combined_correlation_results %>%
  select(
    pair_id,
    metabolite_source,
    metabolite_ion_mode,
    metabolite,
    metabolite_formula,
    metabolite_mz,
    metabolite_rt,
    hormone,
    diet_correlation_group,
    correlation,
    p_value,
    p_adj_bh,
    n
  ) %>%
  pivot_wider(
    id_cols = c(
      pair_id,
      metabolite_source,
      metabolite_ion_mode,
      metabolite,
      metabolite_formula,
      metabolite_mz,
      metabolite_rt,
      hormone
    ),
    names_from = diet_correlation_group,
    values_from = c(correlation, p_value, p_adj_bh, n)
  ) %>%
  rename(
    cd_correlation = correlation_CD,
    cd_p_value = p_value_CD,
    cd_p_adj_bh = p_adj_bh_CD,
    cd_n = n_CD,
    id_correlation = correlation_ID,
    id_p_value = p_value_ID,
    id_p_adj_bh = p_adj_bh_ID,
    id_n = n_ID,
    combined_correlation = correlation_CD_ID_combined,
    combined_p_value = p_value_CD_ID_combined,
    combined_p_adj_bh = p_adj_bh_CD_ID_combined,
    combined_n = n_CD_ID_combined
  ) %>%
  left_join(tp4_gtt_correlation_results, by = "pair_id") %>%
  left_join(sys_sac_correlation_results, by = "pair_id") %>%
  mutate(
    cd_id_same_direction = sign(cd_correlation) == sign(id_correlation) &
      sign(cd_correlation) == sign(combined_correlation),
    state_same_direction = sign(tp4_gtt_correlation) == sign(sys_sac_correlation) &
      sign(tp4_gtt_correlation) == sign(combined_correlation)
  )

conservative_candidate_pairs = pair_level_correlation_summary %>%
  filter(
    !is.na(combined_p_adj_bh),
    combined_p_adj_bh < 0.05,
    !is.na(cd_p_value),
    !is.na(id_p_value),
    cd_p_value < 0.05,
    id_p_value < 0.05,
    cd_id_same_direction,
    !is.na(tp4_gtt_p_value),
    !is.na(sys_sac_p_value),
    tp4_gtt_p_value < 0.05,
    sys_sac_p_value < 0.05,
    state_same_direction
  ) %>%
  arrange(combined_p_adj_bh, metabolite, hormone) %>%
  mutate(
    selection_rank = row_number(),
    selection_reason = paste(
      "Combined FDR < 0.05; CD and ID nominal p < 0.05 with same direction;",
      "TP4/GTT and sys_sac/SAC nominal p < 0.05 with same direction"
    )
  )

heatmap_pairs = conservative_candidate_pairs %>%
  mutate(
    metabolite_label = stringr::str_wrap(metabolite, width = 34),
    pair_label = paste0(metabolite_label, "\n", hormone)
  )

readr::write_csv(
  conservative_candidate_pairs,
  file.path(output_dir, "old serum tp4 sys_sac bicor conservative candidate pairs.csv"),
  na = "N/A"
)

readr::write_csv(
  heatmap_pairs,
  file.path(output_dir, "old serum tp4 sys_sac bicor figure heatmap pairs.csv"),
  na = "N/A"
)


# ========== 4.0 - Build heatmap plot data ==========
# --
heatmap_pair_order = heatmap_pairs %>%
  arrange(desc(combined_correlation), combined_p_adj_bh) %>%
  arrange(hormone, combined_p_adj_bh) %>%
  pull(pair_label)

heatmap_plot_data = bind_rows(
  heatmap_pairs %>%
    transmute(
      pair_id,
      pair_label,
      hormone,
      correlation_column = "CD pooled states",
      correlation = cd_correlation,
      p_value = cd_p_value,
      p_adj_bh = cd_p_adj_bh,
      n = cd_n
    ),
  heatmap_pairs %>%
    transmute(
      pair_id,
      pair_label,
      hormone,
      correlation_column = "ID pooled states",
      correlation = id_correlation,
      p_value = id_p_value,
      p_adj_bh = id_p_adj_bh,
      n = id_n
    ),
  heatmap_pairs %>%
    transmute(
      pair_id,
      pair_label,
      hormone,
      correlation_column = "CD+ID pooled states",
      correlation = combined_correlation,
      p_value = combined_p_value,
      p_adj_bh = combined_p_adj_bh,
      n = combined_n
    ),
  heatmap_pairs %>%
    transmute(
      pair_id,
      pair_label,
      hormone,
      correlation_column = "TP4/GTT CD+ID",
      correlation = tp4_gtt_correlation,
      p_value = tp4_gtt_p_value,
      p_adj_bh = tp4_gtt_p_adj_bh,
      n = tp4_gtt_n
    ),
  heatmap_pairs %>%
    transmute(
      pair_id,
      pair_label,
      hormone,
      correlation_column = "sys_sac/SAC CD+ID",
      correlation = sys_sac_correlation,
      p_value = sys_sac_p_value,
      p_adj_bh = sys_sac_p_adj_bh,
      n = sys_sac_n
    )
) %>%
  mutate(
    pair_label = factor(pair_label, levels = rev(heatmap_pair_order)),
    correlation_column = factor(correlation_column, levels = correlation_column_order),
    significance_level = case_when(
      !is.na(p_adj_bh) & p_adj_bh < 0.05 ~ "BH < 0.05",
      !is.na(p_value) & p_value < 0.05 ~ "p < 0.05",
      TRUE ~ "not significant"
    ),
    correlation_label = sprintf("%.2f", correlation)
  )

readr::write_csv(
  heatmap_plot_data,
  file.path(output_dir, "old serum tp4 sys_sac bicor figure heatmap plot data.csv"),
  na = "N/A"
)


# ========== 5.0 - Reconstruct matched observations for representative scatters ==========
# --
preferred_scatter_pairs = tibble::tribble(
  ~metabolite, ~hormone,
  "Dimethylguanidino valeric acid", "GLP-1",
  "Anserine", "Insulin",
  "Tetrahydrocortisol", "C-Peptide 2",
  "Uracil", "Leptin",
  "L-Valine", "Resistin",
  "5-Hydroxy-L-tryptophan", "C-Peptide 2"
)

scatter_pairs = conservative_candidate_pairs %>%
  inner_join(preferred_scatter_pairs, by = c("metabolite", "hormone")) %>%
  group_by(metabolite, hormone) %>%
  slice_min(combined_p_adj_bh, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(match(metabolite, preferred_scatter_pairs$metabolite))

if (nrow(scatter_pairs) < 6) {
  scatter_pairs = bind_rows(
    scatter_pairs,
    conservative_candidate_pairs %>%
      anti_join(scatter_pairs %>% select(pair_id), by = "pair_id") %>%
      group_by(hormone) %>%
      slice_min(combined_p_adj_bh, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(combined_p_adj_bh)
  ) %>%
    slice_head(n = 6)
}

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
    metabolite_age = factor(age, levels = c("Young", "Old")),
    metabolite_diet = factor(diet, levels = c("CD", "ID"))
  ) %>%
  group_by(
    sample,
    sample_type,
    hormone_sample_id,
    metabolite_sample_id,
    metabolite_timepoint,
    hormone_sampling_time
  ) %>%
  summarise(
    metabolite_diet = first(metabolite_diet),
    metabolite_age = first(metabolite_age),
    .groups = "drop"
  )

serum_metabolite_data = readr::read_csv(
  here::here("data/processed/serum_metabolomics/ion_counts/tidy and normalized serum ion count data.csv"),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  mutate(
    metabolite_abundance = coalesce(med.adj.ic, ion.count),
    metabolite_abundance_source = case_when(
      !is.na(med.adj.ic) ~ "med.adj.ic",
      is.na(med.adj.ic) & !is.na(ion.count) ~ "ion.count",
      TRUE ~ NA_character_
    )
  ) %>%
  transmute(
    metabolite_sample_id = sample.id,
    metabolite_timepoint = timepoint,
    metabolite_source = source,
    metabolite_ion_mode = ion.mode,
    metabolite = name,
    metabolite_formula = formula,
    metabolite_mz = mz,
    metabolite_rt = rt,
    metabolite_abundance,
    metabolite_abundance_source
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
    raw_sample_id = str_trim(raw_sample_id),
    sample = str_extract(raw_sample_id, "^[^_]+"),
    hormone_sampling_time = str_extract(raw_sample_id, "(?<=_).*")
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
  transmute(
    hormone_mouse_id = sample,
    hormone_sample_id = raw_sample_id,
    hormone_sampling_time,
    hormone,
    hormone_abundance
  )

serum_metabolite_samples = serum_metabolite_data %>%
  distinct(metabolite_sample_id, metabolite_timepoint)

hormone_samples = hormone_data %>%
  distinct(
    hormone_sample_id,
    hormone_sampling_time,
    hormone_mouse_id
  )

serum_hormone_sample_matches = serum_sample_inventory %>%
  inner_join(
    serum_metabolite_samples,
    by = c("metabolite_sample_id", "metabolite_timepoint")
  ) %>%
  inner_join(
    hormone_samples,
    by = c("hormone_sample_id", "hormone_sampling_time")
  )

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
    (metabolite_timepoint == "tp4" & hormone_sampling_time == "GTT") |
      (metabolite_timepoint == "sys_sac" & hormone_sampling_time == "SAC"),
    !is.na(metabolite_abundance),
    !is.na(hormone_abundance)
  ) %>%
  add_pair_id() %>%
  mutate(
    matched_state = case_when(
      metabolite_timepoint == "tp4" & hormone_sampling_time == "GTT" ~ "TP4/GTT",
      metabolite_timepoint == "sys_sac" & hormone_sampling_time == "SAC" ~ "sys_sac/SAC",
      TRUE ~ paste(metabolite_timepoint, hormone_sampling_time, sep = "/")
    )
  )

scatter_plot_data = matched_observations %>%
  semi_join(scatter_pairs %>% select(pair_id), by = "pair_id") %>%
  inner_join(
    scatter_pairs %>%
      select(
        pair_id,
        combined_correlation,
        combined_p_adj_bh,
        cd_correlation,
        id_correlation,
        tp4_gtt_correlation,
        sys_sac_correlation
      ),
    by = "pair_id"
  ) %>%
  filter(
    metabolite_abundance > 0,
    hormone_abundance > 0
  ) %>%
  mutate(
    log10_metabolite_abundance = log10(metabolite_abundance),
    log10_hormone_abundance = log10(hormone_abundance),
    scatter_label = paste0(
      stringr::str_wrap(metabolite, width = 28),
      "\n",
      hormone,
      "\nbicor = ",
      sprintf("%.2f", combined_correlation),
      ", q = ",
      scales::scientific(combined_p_adj_bh, digits = 2)
    ),
    scatter_label = factor(
      scatter_label,
      levels = scatter_pairs %>%
        mutate(
          scatter_label = paste0(
            stringr::str_wrap(metabolite, width = 28),
            "\n",
            hormone,
            "\nbicor = ",
            sprintf("%.2f", combined_correlation),
            ", q = ",
            scales::scientific(combined_p_adj_bh, digits = 2)
          )
        ) %>%
        pull(scatter_label)
    )
  )

readr::write_csv(
  scatter_pairs,
  file.path(output_dir, "old serum tp4 sys_sac bicor figure scatter pairs.csv"),
  na = "N/A"
)

readr::write_csv(
  scatter_plot_data,
  file.path(output_dir, "old serum tp4 sys_sac bicor figure scatter plot data.csv"),
  na = "N/A"
)


# ========== 6.0 - Make figure panels ==========
# --
count_plot = ggplot(
  count_plot_data,
  aes(x = signed_hits, y = hormone, fill = correlation_direction)
) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#374151") +
  facet_wrap(~ diet_correlation_group, ncol = 1, scales = "free_x") +
  scale_x_continuous(labels = abs) +
  scale_fill_manual(
    values = c(positive = "#B3483C", negative = "#3B6EA8"),
    breaks = c("positive", "negative"),
    labels = c("Positive", "Negative")
  ) +
  labs(
    title = "FDR-significant bicor hits",
    subtitle = "TP4/GTT plus sys_sac/SAC matched observations",
    x = "BH-significant pairs",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F3F4F6", color = "#D1D5DB"),
    strip.text = element_text(face = "bold", color = "#111827"),
    axis.text.y = element_text(size = 6.8, color = "#111827"),
    axis.text.x = element_text(size = 6.5, color = "#111827"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", color = "#111827"),
    plot.subtitle = element_text(size = 7.5, color = "#4B5563")
  )

heatmap_plot = ggplot(
  heatmap_plot_data,
  aes(x = correlation_column, y = pair_label)
) +
  geom_point(
    aes(fill = correlation, size = significance_level),
    shape = 21,
    color = "#1F2937",
    stroke = 0.25,
    alpha = 0.95
  ) +
  geom_text(
    aes(label = correlation_label),
    size = 2.1,
    color = "#111827"
  ) +
  scale_fill_gradient2(
    low = "#3B6EA8",
    mid = "white",
    high = "#B3483C",
    midpoint = 0,
    limits = correlation_fill_limits,
    name = "bicor"
  ) +
  scale_size_manual(
    values = c(
      "BH < 0.05" = 5.2,
      "p < 0.05" = 4.1,
      "not significant" = 2.4
    ),
    breaks = c("BH < 0.05", "p < 0.05", "not significant"),
    name = NULL
  ) +
  labs(
    title = "Conservative reportable candidates",
    subtitle = paste(
      "Combined FDR < 0.05; same direction and nominal support in CD, ID,",
      "TP4/GTT, and sys_sac/SAC"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_line(color = "#E5E7EB", linewidth = 0.25),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "#111827", size = 7),
    axis.text.y = element_text(color = "#111827", size = 5.7),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", color = "#111827"),
    plot.subtitle = element_text(size = 7.5, color = "#4B5563")
  )

scatter_plot = ggplot(
  scatter_plot_data,
  aes(x = log10_metabolite_abundance, y = log10_hormone_abundance)
) +
  geom_point(
    aes(color = metabolite_diet, shape = matched_state),
    size = 2,
    alpha = 0.82
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "#111827",
    linewidth = 0.35
  ) +
  facet_wrap(~ scatter_label, scales = "free", ncol = 3) +
  scale_color_manual(values = c(CD = "#4E79A7", ID = "#F28E2B")) +
  scale_shape_manual(values = c("TP4/GTT" = 16, "sys_sac/SAC" = 17)) +
  labs(
    title = "Representative raw matched observations",
    x = "log10 serum metabolite abundance",
    y = "log10 hormone abundance",
    color = "Diet",
    shape = "Matched state"
  ) +
  theme_bw(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "#F3F4F6", color = "#D1D5DB"),
    strip.text = element_text(face = "bold", color = "#111827", size = 6.2),
    axis.text = element_text(size = 6, color = "#111827"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", color = "#111827")
  )

combined_figure = (
  count_plot + heatmap_plot +
    plot_layout(widths = c(0.72, 1.45), guides = "collect")
) /
  scatter_plot +
  plot_layout(heights = c(1.18, 1), guides = "collect") +
  plot_annotation(
    title = "Old Serum Hormone-Metabolite Bicor Correlations",
    subtitle = "Metabolite TP4 and systemic sacrifice timepoints matched to GTT and sacrifice hormone measurements",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, color = "#111827"),
      plot.subtitle = element_text(size = 9, color = "#4B5563")
    )
  )


# ========== 7.0 - Save figures ==========
# --
combined_pdf_file = file.path(
  plot_dir,
  "bicor tp4 sys_sac hormone metabolite reportable correlations figure.pdf"
)
combined_png_file = file.path(
  plot_dir,
  "bicor tp4 sys_sac hormone metabolite reportable correlations figure.png"
)
summary_pdf_file = file.path(
  plot_dir,
  "bicor tp4 sys_sac hormone metabolite count and heatmap.pdf"
)
scatter_pdf_file = file.path(
  plot_dir,
  "bicor tp4 sys_sac hormone metabolite representative scatters.pdf"
)

ggsave(
  combined_pdf_file,
  combined_figure,
  width = 16,
  height = 14,
  device = grDevices::cairo_pdf
)
ggsave(
  combined_png_file,
  combined_figure,
  width = 16,
  height = 14,
  dpi = 300
)
ggsave(
  summary_pdf_file,
  count_plot + heatmap_plot + plot_layout(widths = c(0.72, 1.45), guides = "collect"),
  width = 16,
  height = 9.5,
  device = grDevices::cairo_pdf
)
ggsave(
  scatter_pdf_file,
  scatter_plot,
  width = 11,
  height = 8.5,
  device = grDevices::cairo_pdf
)

message("Wrote candidate table to: ", file.path(output_dir, "old serum tp4 sys_sac bicor conservative candidate pairs.csv"))
message("Wrote combined figure PDF to: ", combined_pdf_file)
message("Wrote combined figure PNG to: ", combined_png_file)
message("Wrote summary PDF to: ", summary_pdf_file)
message("Wrote scatter PDF to: ", scatter_pdf_file)
