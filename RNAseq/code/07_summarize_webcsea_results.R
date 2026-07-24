library(here)
library(tidyverse)

here::i_am("code/07_summarize_webcsea_results.R")


# ========== 0.0 - Analysis settings ==========
webcsea_dir = here("data/processed/webCSEA")
output_dir = here("data/processed/webCSEA/summary")
plot_dir = here("plots/webCSEA")
top_n_per_direction = 25
top_primary_tissues_to_plot = 15
top_primary_cell_types_to_plot = 15
developmental_stages_to_plot = c("Adult", "Fetal")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# ========== 1.0 - Find webCSEA result folders ==========
webcsea_jobs = list.dirs(webcsea_dir, recursive = FALSE, full.names = TRUE)

webcsea_jobs = webcsea_jobs[file.exists(file.path(webcsea_jobs, "WebCSEA_All_tissue_cell_type.txt"))]

if (length(webcsea_jobs) == 0) {
  stop("No webCSEA job folders found under: ", webcsea_dir)
}


# ========== 2.0 - Read and combine tissue/cell-type enrichment results ==========
read_webcsea_job = function(job_dir) {
  job_name = basename(job_dir)
  
  direction = case_when(
    str_detect(str_to_lower(job_name), "up") ~ "Upregulated",
    str_detect(str_to_lower(job_name), "down") ~ "Downregulated",
    TRUE ~ "Unknown"
  )
  
  read_tsv(
    file.path(job_dir, "WebCSEA_All_tissue_cell_type.txt"),
    show_col_types = FALSE
  ) %>%
    mutate(
      webcsea_job = job_name,
      direction = direction,
      .before = 1
    )
}

read_webcsea_gene_hits = function(job_dir) {
  job_name = basename(job_dir)
  
  direction = case_when(
    str_detect(str_to_lower(job_name), "up") ~ "Upregulated",
    str_detect(str_to_lower(job_name), "down") ~ "Downregulated",
    TRUE ~ "Unknown"
  )
  
  read_tsv(
    file.path(job_dir, "WebCSEA_all_gene_list_result.tsv"),
    skip = 1,
    col_names = c("Tissue_cell_type_name", "human_symbol"),
    show_col_types = FALSE
  ) %>%
    separate_rows(human_symbol, sep = ";") %>%
    mutate(
      human_symbol = str_trim(human_symbol),
      webcsea_job = job_name,
      direction = direction,
      .before = 1
    ) %>%
    dplyr::filter(human_symbol != "")
}

webcsea_results = map_dfr(webcsea_jobs, read_webcsea_job) %>%
  arrange(direction, input_list_combined_p, desc(input_list_combined_log10p))

write.csv(
  webcsea_results,
  file.path(output_dir, "combined webCSEA tissue cell type results.csv"),
  row.names = FALSE
)


# ========== 3.0 - Save top results per direction ==========
top_webcsea_results = webcsea_results %>%
  group_by(direction) %>%
  slice_min(input_list_combined_p, n = top_n_per_direction, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    direction,
    Tissue_cell_type_name,
    Developmental_stage,
    Organ_system,
    General_tissue_name,
    Specific_cell_type,
    General_cell_type,
    input_list_combined_p,
    input_list_combined_log10p,
    input_list_raw_p,
    input_list_permutated_p,
    input_list_rank_p,
    webcsea_job
  )

write.csv(
  top_webcsea_results,
  file.path(output_dir, "top webCSEA tissue cell type results by direction.csv"),
  row.names = FALSE
)


# ========== 4.0 - Save simple summaries by organ system and general cell type ==========
organ_system_summary = webcsea_results %>%
  group_by(direction, Organ_system) %>%
  summarise(
    best_combined_p = min(input_list_combined_p, na.rm = TRUE),
    best_combined_log10p = max(input_list_combined_log10p, na.rm = TRUE),
    significant_cell_type_n = sum(input_list_combined_p <= 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(direction, best_combined_p)

general_cell_type_summary = webcsea_results %>%
  group_by(direction, General_cell_type) %>%
  summarise(
    best_combined_p = min(input_list_combined_p, na.rm = TRUE),
    best_combined_log10p = max(input_list_combined_log10p, na.rm = TRUE),
    significant_tissue_cell_type_n = sum(input_list_combined_p <= 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(direction, best_combined_p)

write.csv(
  organ_system_summary,
  file.path(output_dir, "webCSEA organ system summary.csv"),
  row.names = FALSE
)

write.csv(
  general_cell_type_summary,
  file.path(output_dir, "webCSEA general cell type summary.csv"),
  row.names = FALSE
)


# ========== 5.0 - Assign each gene to its strongest tissue/cell-type hit ==========
# This is an exploratory "primary webCSEA hit" assignment. It does not prove the
# source tissue of an individual gene; it identifies the strongest webCSEA
# tissue/cell-type signature containing that gene within each DEG direction.

webcsea_gene_hits = map_dfr(webcsea_jobs, read_webcsea_gene_hits)

gene_primary_webcsea_hit = webcsea_gene_hits %>%
  left_join(
    webcsea_results,
    by = c("webcsea_job", "direction", "Tissue_cell_type_name")
  ) %>%
  arrange(
    direction,
    human_symbol,
    input_list_combined_p,
    desc(input_list_combined_log10p),
    input_list_raw_p
  ) %>%
  group_by(direction, human_symbol) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(
    direction,
    Developmental_stage,
    human_symbol,
    primary_tissue = General_tissue_name,
    primary_organ_system = Organ_system,
    primary_specific_cell_type = Specific_cell_type,
    primary_general_cell_type = General_cell_type,
    primary_tissue_cell_type_name = Tissue_cell_type_name,
    input_list_combined_p,
    input_list_combined_log10p,
    input_list_raw_p,
    input_list_permutated_p,
    input_list_rank_p,
    webcsea_job
  )

write.csv(
  gene_primary_webcsea_hit,
  file.path(output_dir, "gene primary webCSEA tissue hit assignments.csv"),
  row.names = FALSE
)

gene_primary_webcsea_hit_by_stage = webcsea_gene_hits %>%
  left_join(
    webcsea_results,
    by = c("webcsea_job", "direction", "Tissue_cell_type_name")
  ) %>%
  dplyr::filter(Developmental_stage %in% developmental_stages_to_plot) %>%
  arrange(
    direction,
    Developmental_stage,
    human_symbol,
    input_list_combined_p,
    desc(input_list_combined_log10p),
    input_list_raw_p
  ) %>%
  group_by(direction, Developmental_stage, human_symbol) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(
    direction,
    Developmental_stage,
    human_symbol,
    primary_tissue = General_tissue_name,
    primary_organ_system = Organ_system,
    primary_specific_cell_type = Specific_cell_type,
    primary_general_cell_type = General_cell_type,
    primary_tissue_cell_type_name = Tissue_cell_type_name,
    input_list_combined_p,
    input_list_combined_log10p,
    input_list_raw_p,
    input_list_permutated_p,
    input_list_rank_p,
    webcsea_job
  )

write.csv(
  gene_primary_webcsea_hit_by_stage,
  file.path(output_dir, "gene primary webCSEA tissue hit assignments adult fetal.csv"),
  row.names = FALSE
)

gene_muscle_webcsea_hits = webcsea_gene_hits %>%
  left_join(
    webcsea_results,
    by = c("webcsea_job", "direction", "Tissue_cell_type_name")
  ) %>%
  dplyr::filter(str_detect(General_tissue_name, regex("muscle", ignore_case = TRUE))) %>%
  arrange(direction, human_symbol, input_list_combined_p) %>%
  group_by(direction, human_symbol) %>%
  summarise(
    best_muscle_combined_p = min(input_list_combined_p, na.rm = TRUE),
    best_muscle_combined_log10p = max(input_list_combined_log10p, na.rm = TRUE),
    muscle_tissue_cell_type_hits = paste(unique(Tissue_cell_type_name), collapse = "; "),
    muscle_specific_cell_type_hits = paste(unique(Specific_cell_type), collapse = "; "),
    muscle_general_cell_type_hits = paste(unique(General_cell_type), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(direction, best_muscle_combined_p, human_symbol)

write.csv(
  gene_muscle_webcsea_hits,
  file.path(output_dir, "genes with any webCSEA muscle tissue hit.csv"),
  row.names = FALSE
)

gene_significant_muscle_webcsea_hits = gene_muscle_webcsea_hits %>%
  dplyr::filter(best_muscle_combined_p <= 0.05)

write.csv(
  gene_significant_muscle_webcsea_hits,
  file.path(output_dir, "genes with significant webCSEA muscle tissue hit.csv"),
  row.names = FALSE
)

gene_significant_muscle_webcsea_hits_by_stage = webcsea_gene_hits %>%
  left_join(
    webcsea_results,
    by = c("webcsea_job", "direction", "Tissue_cell_type_name")
  ) %>%
  dplyr::filter(
    Developmental_stage %in% developmental_stages_to_plot,
    str_detect(General_tissue_name, regex("muscle", ignore_case = TRUE)),
    input_list_combined_p <= 0.05
  ) %>%
  arrange(direction, Developmental_stage, human_symbol, input_list_combined_p) %>%
  group_by(direction, Developmental_stage, human_symbol) %>%
  summarise(
    best_muscle_combined_p = min(input_list_combined_p, na.rm = TRUE),
    best_muscle_combined_log10p = max(input_list_combined_log10p, na.rm = TRUE),
    muscle_tissue_cell_type_hits = paste(unique(Tissue_cell_type_name), collapse = "; "),
    muscle_specific_cell_type_hits = paste(unique(Specific_cell_type), collapse = "; "),
    muscle_general_cell_type_hits = paste(unique(General_cell_type), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(direction, Developmental_stage, best_muscle_combined_p, human_symbol)

write.csv(
  gene_significant_muscle_webcsea_hits_by_stage,
  file.path(output_dir, "genes with significant webCSEA muscle tissue hit adult fetal.csv"),
  row.names = FALSE
)

primary_tissue_tally = gene_primary_webcsea_hit %>%
  count(direction, primary_tissue, name = "gene_n") %>%
  group_by(direction) %>%
  mutate(
    direction_total_genes_with_hit = sum(gene_n),
    fraction_of_genes_with_hit = gene_n / direction_total_genes_with_hit
  ) %>%
  ungroup() %>%
  arrange(direction, desc(gene_n), primary_tissue)

write.csv(
  primary_tissue_tally,
  file.path(output_dir, "gene primary webCSEA tissue tally.csv"),
  row.names = FALSE
)

primary_tissue_tally_by_stage = gene_primary_webcsea_hit_by_stage %>%
  count(direction, Developmental_stage, primary_tissue, name = "gene_n") %>%
  group_by(direction, Developmental_stage) %>%
  mutate(
    direction_stage_total_genes_with_hit = sum(gene_n),
    fraction_of_direction_stage_genes_with_hit = gene_n / direction_stage_total_genes_with_hit
  ) %>%
  ungroup() %>%
  arrange(direction, Developmental_stage, desc(gene_n), primary_tissue)

primary_cell_type_tally_by_stage = gene_primary_webcsea_hit_by_stage %>%
  mutate(primary_general_cell_type = replace_na(primary_general_cell_type, "Unlabeled")) %>%
  count(direction, Developmental_stage, primary_general_cell_type, name = "gene_n") %>%
  group_by(direction, Developmental_stage) %>%
  mutate(
    direction_stage_total_genes_with_hit = sum(gene_n),
    fraction_of_direction_stage_genes_with_hit = gene_n / direction_stage_total_genes_with_hit
  ) %>%
  ungroup() %>%
  arrange(direction, Developmental_stage, desc(gene_n), primary_general_cell_type)

write.csv(
  primary_tissue_tally_by_stage,
  file.path(output_dir, "gene primary webCSEA tissue tally adult fetal.csv"),
  row.names = FALSE
)

write.csv(
  primary_cell_type_tally_by_stage,
  file.path(output_dir, "gene primary webCSEA general cell type tally adult fetal.csv"),
  row.names = FALSE
)


# ========== 6.0 - Plot primary tissue tallies ==========
primary_tissue_plot_data = primary_tissue_tally %>%
  group_by(direction) %>%
  slice_max(gene_n, n = top_primary_tissues_to_plot, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(direction, gene_n) %>%
  mutate(
    primary_tissue_for_plot = paste(primary_tissue, direction, sep = "___"),
    primary_tissue_for_plot = factor(primary_tissue_for_plot, levels = unique(primary_tissue_for_plot))
  )

primary_tissue_tally_plot = ggplot(
  primary_tissue_plot_data,
  aes(gene_n, primary_tissue_for_plot, fill = direction)
) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.25) +
  geom_text(
    aes(label = gene_n),
    hjust = -0.18,
    size = 3.2,
    color = "black"
  ) +
  facet_wrap(~ direction, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_manual(
    values = c(
      Downregulated = "#B04A5A",
      Upregulated = "#2C7FB8"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Primary webCSEA tissue hit per nominal DEG",
    subtitle = "Each human ortholog DEG is assigned to the strongest matching webCSEA tissue/cell-type result",
    x = "Number of genes",
    y = NULL
  ) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey45"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  here("plots/webCSEA/primary webCSEA tissue hit tally among nominal DEGs.pdf"),
  primary_tissue_tally_plot,
  width = 10,
  height = 7
)

primary_tissue_plot_data_by_stage = primary_tissue_tally_by_stage %>%
  group_by(direction, Developmental_stage) %>%
  slice_max(gene_n, n = top_primary_tissues_to_plot, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Developmental_stage = factor(Developmental_stage, levels = developmental_stages_to_plot),
    direction = factor(direction, levels = c("Upregulated", "Downregulated"))
  ) %>%
  arrange(Developmental_stage, direction, gene_n) %>%
  mutate(
    primary_tissue_for_plot = paste(primary_tissue, Developmental_stage, direction, sep = "___"),
    primary_tissue_for_plot = factor(primary_tissue_for_plot, levels = unique(primary_tissue_for_plot))
  )

primary_tissue_tally_by_stage_plot = ggplot(
  primary_tissue_plot_data_by_stage,
  aes(gene_n, primary_tissue_for_plot, fill = direction)
) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.25) +
  geom_text(
    aes(label = gene_n),
    hjust = -0.18,
    size = 2.8,
    color = "black"
  ) +
  facet_grid(Developmental_stage ~ direction, scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_manual(
    values = c(
      Downregulated = "#B04A5A",
      Upregulated = "#2C7FB8"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "Primary webCSEA organ/tissue hit per nominal DEG",
    subtitle = "Primary hit assigned separately within Adult and Fetal webCSEA results",
    x = "Number of genes",
    y = NULL
  ) +
  theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey45"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  here("plots/webCSEA/primary webCSEA organ tissue hit tally adult fetal nominal DEGs.pdf"),
  primary_tissue_tally_by_stage_plot,
  width = 11,
  height = 10
)

primary_cell_type_plot_data_by_stage = primary_cell_type_tally_by_stage %>%
  group_by(direction, Developmental_stage) %>%
  slice_max(gene_n, n = top_primary_cell_types_to_plot, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Developmental_stage = factor(Developmental_stage, levels = developmental_stages_to_plot),
    direction = factor(direction, levels = c("Upregulated", "Downregulated"))
  ) %>%
  arrange(Developmental_stage, direction, gene_n) %>%
  mutate(
    primary_cell_type_for_plot = paste(
      primary_general_cell_type,
      Developmental_stage,
      direction,
      sep = "___"
    ),
    primary_cell_type_for_plot = factor(
      primary_cell_type_for_plot,
      levels = unique(primary_cell_type_for_plot)
    )
  )

primary_cell_type_tally_by_stage_plot = ggplot(
  primary_cell_type_plot_data_by_stage,
  aes(gene_n, primary_cell_type_for_plot, fill = direction)
) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.25) +
  geom_text(
    aes(label = gene_n),
    hjust = -0.18,
    size = 2.8,
    color = "black"
  ) +
  facet_grid(Developmental_stage ~ direction, scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_manual(
    values = c(
      Downregulated = "#B04A5A",
      Upregulated = "#2C7FB8"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "Primary webCSEA general cell-type hit per nominal DEG",
    subtitle = "Primary hit assigned separately within Adult and Fetal webCSEA results",
    x = "Number of genes",
    y = NULL
  ) +
  theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey45"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  here("plots/webCSEA/primary webCSEA general cell type hit tally adult fetal nominal DEGs.pdf"),
  primary_cell_type_tally_by_stage_plot,
  width = 11,
  height = 10
)


# ========== 7.0 - Console summary ==========
cat("\nTop webCSEA results by direction:\n")
top_webcsea_results %>%
  group_by(direction) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  select(
    direction,
    Tissue_cell_type_name,
    Organ_system,
    General_tissue_name,
    General_cell_type,
    input_list_combined_p,
    input_list_combined_log10p
  ) %>%
  print(n = Inf)

cat("\nTop primary webCSEA tissue assignments by direction:\n")
primary_tissue_tally %>%
  group_by(direction) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  print(n = Inf)

cat("\nTop primary webCSEA organ/tissue assignments by direction and developmental stage:\n")
primary_tissue_tally_by_stage %>%
  group_by(direction, Developmental_stage) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  print(n = Inf)

cat("\nTop primary webCSEA general cell-type assignments by direction and developmental stage:\n")
primary_cell_type_tally_by_stage %>%
  group_by(direction, Developmental_stage) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  print(n = Inf)

cat("\nGenes with significant webCSEA Muscle organ hit by direction:\n")
gene_significant_muscle_webcsea_hits %>%
  count(direction, name = "gene_n") %>%
  print(n = Inf)

cat("\nGenes with significant webCSEA Muscle organ hit by direction and developmental stage:\n")
gene_significant_muscle_webcsea_hits_by_stage %>%
  count(direction, Developmental_stage, name = "gene_n") %>%
  print(n = Inf)
