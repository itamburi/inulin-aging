# ****** Inulin Aging Mouse Serum Hormone-Metabolite Bicor U Plots
# ****** Summer 2026 - Ian T, Jessie K, Randa L
#
# Figure focus:
# - Old mice only.
# - CD+ID combined bicor correlations.
# - Serum TP4 metabolites matched to GTT hormones.
# - Serum systemic sacrifice metabolites matched to SAC hormones.
# - Nominal p-values are used for the U plot threshold.
# - Carnitine-named metabolites are highlighted.


here::i_am("code/integrative/05_plot_serum_hormone_metabolite_bicor_uplots.R")


# ========== 0.0 - Environment ==========
# --
required_packages = c(
  "dplyr", "readr", "stringr", "ggplot2", "ggrepel", "scales"
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
library(readr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(scales)

output_dir = here::here("data/processed/integrative/serum_hormone_metabolite_correlations")
plot_dir = here::here("plots/integrative/serum_hormone_metabolite_correlations")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Cairo PDF needs a writable font cache when this script runs in a sandboxed shell.
Sys.setenv(XDG_CACHE_HOME = tempdir())

hormone_columns = c(
  "C-Peptide 2", "Ghrelin", "GIP", "GLP-1", "Glucagon", "Insulin",
  "Leptin", "Secretin", "PP", "PYY", "Resistin", "Amylin"
)

nominal_p_cutoff = 0.05
nominal_p_line = -log10(nominal_p_cutoff)


# ========== 1.0 - Read full combined bicor results ==========
# --
essential_results = readr::read_csv(
  file.path(
    output_dir,
    "old serum tp4 sys_sac matched gtt sac CD ID combined bicor essential results.csv"
  ),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
)

u_plot_data = essential_results %>%
  filter(
    analysis_label == "metabolite TP4 and sys_sac vs hormone GTT and SAC",
    diet_correlation_group == "CD_ID_combined",
    correlation_method == "bicor",
    !is.na(correlation),
    !is.na(p_value)
  ) %>%
  mutate(
    hormone = factor(hormone, levels = hormone_columns),
    nominal_p_significant = p_value < nominal_p_cutoff,
    carnitine_metabolite = stringr::str_detect(
      metabolite,
      stringr::regex("carnitine", ignore_case = TRUE)
    ),
    neg_log10_p_value = -log10(pmax(p_value, .Machine$double.xmin)),
    u_plot_group = case_when(
      carnitine_metabolite & nominal_p_significant ~ "Carnitine, p < 0.05",
      carnitine_metabolite ~ "Carnitine, p >= 0.05",
      TRUE ~ "Other metabolites"
    ),
    u_plot_group = factor(
      u_plot_group,
      levels = c(
        "Other metabolites",
        "Carnitine, p >= 0.05",
        "Carnitine, p < 0.05"
      )
    )
  ) %>%
  arrange(u_plot_group, p_value)

if (nrow(u_plot_data) == 0) {
  stop("No bicor rows available for U plot generation.", call. = FALSE)
}

carnitine_summary = u_plot_data %>%
  summarise(
    n_total_rows = n(),
    n_carnitine_rows = sum(carnitine_metabolite),
    n_unique_carnitines = n_distinct(metabolite[carnitine_metabolite]),
    n_nominal_carnitine_rows = sum(carnitine_metabolite & nominal_p_significant)
  )

carnitine_label_data = bind_rows(
  u_plot_data %>%
    filter(carnitine_metabolite, correlation > 0) %>%
    group_by(hormone) %>%
    slice_max(correlation, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(carnitine_label_direction = "Top positive carnitines"),
  u_plot_data %>%
    filter(carnitine_metabolite, correlation < 0) %>%
    group_by(hormone) %>%
    slice_min(correlation, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(carnitine_label_direction = "Top negative carnitines")
) %>%
  mutate(
    metabolite_label = stringr::str_wrap(metabolite, width = 18)
  )

carnitine_facet_counts = u_plot_data %>%
  filter(carnitine_metabolite) %>%
  group_by(hormone) %>%
  summarise(
    n_nominal_significant_carnitines = sum(nominal_p_significant),
    n_carnitines = n(),
    .groups = "drop"
  ) %>%
  mutate(
    carnitine_count_label = paste0(
      n_nominal_significant_carnitines,
      "/",
      n_carnitines,
      " carnitines"
    )
  )


# ========== 2.0 - Make multi-panel U plot ==========
# --
u_plot = ggplot(
  u_plot_data,
  aes(x = correlation, y = neg_log10_p_value)
) +
  geom_hline(
    yintercept = nominal_p_line,
    linewidth = 0.25,
    linetype = "dashed",
    color = "#6B7280"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.2,
    color = "#9CA3AF"
  ) +
  geom_point(
    aes(color = u_plot_group, size = u_plot_group, alpha = u_plot_group)
  ) +
  geom_text(
    data = carnitine_facet_counts,
    aes(
      x = -0.96,
      y = Inf,
      label = carnitine_count_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.35,
    size = 2.4,
    color = "#374151",
    fontface = "bold"
  ) +
  ggrepel::geom_text_repel(
    data = carnitine_label_data,
    aes(label = metabolite_label),
    size = 1.8,
    color = "#111827",
    segment.color = "#6B7280",
    segment.size = 0.18,
    min.segment.length = 0,
    box.padding = 0.18,
    point.padding = 0.08,
    force = 1.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_wrap(~ hormone, ncol = 4) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08)),
    labels = scales::label_number(accuracy = 1)
  ) +
  scale_color_manual(
    values = c(
      "Other metabolites" = "#9CA3AF",
      "Carnitine, p >= 0.05" = "#F59E0B",
      "Carnitine, p < 0.05" = "#B45309"
    ),
    name = NULL
  ) +
  scale_size_manual(
    values = c(
      "Other metabolites" = 0.45,
      "Carnitine, p >= 0.05" = 1.15,
      "Carnitine, p < 0.05" = 1.75
    ),
    name = NULL
  ) +
  scale_alpha_manual(
    values = c(
      "Other metabolites" = 0.26,
      "Carnitine, p >= 0.05" = 0.85,
      "Carnitine, p < 0.05" = 0.95
    ),
    name = NULL
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1)),
    size = "none",
    alpha = "none"
  ) +
  labs(
    title = "Old Serum Hormone-Metabolite Bicor U Plots",
    subtitle = "CD+ID combined TP4/GTT and sys_sac/SAC correlations; dashed line marks nominal p = 0.05",
    x = "bicor correlation",
    y = "-log10(nominal p-value)"
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.25),
    strip.background = element_rect(fill = "#F3F4F6", color = "#D1D5DB"),
    strip.text = element_text(face = "bold", color = "#111827"),
    axis.text = element_text(size = 6.5, color = "#111827"),
    axis.title = element_text(size = 8, color = "#111827"),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 12, color = "#111827"),
    plot.subtitle = element_text(size = 8, color = "#4B5563")
  )


# ========== 3.0 - Save PDF ==========
# --
u_plot_pdf_file = file.path(
  plot_dir,
  "bicor tp4 sys_sac hormone metabolite nominal p u plots carnitines.pdf"
)

ggsave(
  u_plot_pdf_file,
  u_plot,
  width = 15,
  height = 11,
  device = grDevices::cairo_pdf
)

message("Wrote U plot PDF to: ", u_plot_pdf_file)
message(
  "Carnitine rows highlighted: ",
  carnitine_summary$n_carnitine_rows,
  " rows across ",
  carnitine_summary$n_unique_carnitines,
  " unique metabolites; ",
  carnitine_summary$n_nominal_carnitine_rows,
  " rows have nominal p < 0.05."
)
message(
  "Carnitine labels added: ",
  nrow(carnitine_label_data),
  " rows, up to 5 positive and 5 negative bicor tails per hormone."
)
