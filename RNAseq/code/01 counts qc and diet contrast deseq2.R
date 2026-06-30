
library(here)
library(tidyverse)
library(edgeR)
library(cowplot)

here::i_am("code/01 counts qc and diet contrast deseq2.R")

#library(biomaRt)
# First translate ENSMUST ids to gene symbols
# ... biomart is having major issues on this computer so went to marcus sysgen server to make a translation table

# ========== 0.0 - translate ENSMUSG IDs to gene symbols ==========
x = read.delim(here("data/raw/rsem_gene_counts_matrix.tsv"))
# gene_id is ENSMUSG stable ids with .version suffix
length(unique(x$gene_id)) 
length(unique(gsub("\\.\\d+", "", x$gene_id)))
# ... if we remove the suffix the number of unique ENSMUSG stable ids are the same
# Remove suffix and translate to gene symbol

#library(ensembldb)
library(EnsDb.Mmusculus.v79)

x = read.delim(here("data/raw/rsem_gene_counts_matrix.tsv")) %>%
  mutate(
    gene_id = gsub("\\.\\d+", "", gene_id)
  )

ensembl_ids = unique(x$gene_id)

gene_names <- ensembldb::select(
  EnsDb.Mmusculus.v79,
  keys = ensembl_ids,
  keytype = "GENEID",  # Transcript IDs
  columns = c("GENEID", "SYMBOL")
)
names(gene_names) = c("ensembl_id","gene_id")

# note a couple genes have more than one ENSMUSG stable ID
ens = unique(gene_names$ensembl_id)
sym = unique(gene_names$gene_id)
n = gene_names %>%
  group_by(gene_id) %>%
  summarise(n=n()) %>%
  dplyr::filter(n>1)

cnts = x %>%
  pivot_longer(names_to = "sample", cols = !c(gene_id), values_to = "counts") %>%
  dplyr::rename(ensembl_id = "gene_id") %>%
  left_join(., gene_names) %>%
  group_by(sample, gene_id) %>%
  summarise(
    # for the handful of genes with more than one ENSMUSG, sum the counts
    counts = as.integer(sum(counts, na.rm = TRUE))
  ) %>%
  dplyr::filter(gene_id != "") %>%
  pivot_wider(names_from = "sample", values_from = "counts") 
dir.create(here("data/processed/counts"), recursive = TRUE, showWarnings = FALSE)
write.csv(cnts,here("data/processed/counts/raw counts matrix with mouse gene symbols.csv"),row.names = FALSE)





# ========== 1.0 - Filter and PCA ==========
meta = read.csv(here("data/metadata/metadata.csv")) %>%
  mutate(sample = paste0("r",sample)) %>%
  arrange(diet, age.grp) %>%
  column_to_rownames(var = "sample")

x = read.csv(here("data/processed/counts/raw counts matrix with mouse gene symbols.csv"))
cnts = x %>% column_to_rownames(var="gene_id")
cnts <- cnts[, rownames(meta)] 

# counts: genes x samples
dge <- DGEList(cnts, group = meta$diet)
dge$samples <- cbind(dge$samples, meta) # embed all metadata cols into dge obj

# Keep genes with CPM above threshold in enough samples
keep <- filterByExpr(dge, group = meta[["diet"]])   # design-aware, better than manual variance filtering
dge <- dge[keep, , keep.lib.sizes = FALSE]

# see if filtering low variance genes matters
gene_var <- apply(dge$counts, 1, var)
hist(log10(gene_var), breaks = 50)
summary(gene_var)
threshold <- quantile(gene_var, 0.10)  # 10th percentile
dge = dge[gene_var > threshold, , keep.lib.sizes = FALSE]

gene_var <- apply(dge$counts, 1, var)
hist(log10(gene_var), breaks = 50)
summary(gene_var)

dge <- calcNormFactors(dge) 
df = dge$samples


# PCA of unfiltered
mx = x %>%
  column_to_rownames(var="gene_id") %>%
  as.matrix() %>%
  t()
# remove zero variance cols (colsums == 0)
mx = mx[, -which(colSums(mx) == 0)]

meta = meta %>% rownames_to_column(var="sample")

pca = prcomp(mx, scale. = TRUE)
pca_var = round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc_contributions = as.data.frame(pca$x)
pc = pc_contributions %>%
  rownames_to_column(var="sample") %>%
  dplyr::select(sample, PC1, PC2, PC3, PC4, PC5, PC6, PC7,PC8) %>%
  left_join(., meta) %>%
  mutate(
    pca.label = paste0(sample, "\n", age.grp, " ", diet, ", C", cohort),
    pca_set = paste0(
      "Unfiltered counts\n",
      "PC1 = ", pca_var[1], "%, PC2 = ", pca_var[2], "%"
    )
  )
pca_unfiltered = pc

pc %>%
  ggplot(., aes(PC1, PC2, color = diet, shape = age.grp)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    aes(label = pca.label),
    color = "grey20",
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf,
    segment.alpha = 0.4,
    show.legend = FALSE
  ) +
  labs(
    title = "Pre-filtered PCA by diet, age, and cohort",
    color = "Diet",
    shape = "Age"
  ) +
  theme_bw()

pc %>%
  ggplot(.,aes(PC1, PC2, color = diet) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  ggrepel::geom_text_repel(
    aes(label = pca.label),
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf,
    segment.alpha = 0.4,
    show.legend = FALSE
  ) +
  theme_bw()
pc %>%
  ggplot(.,aes(PC1, PC2, color = age.grp) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  ggrepel::geom_text_repel(
    aes(label = pca.label),
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf,
    segment.alpha = 0.4,
    show.legend = FALSE
  ) +
  theme_bw()
pc %>%
  ggplot(.,aes(PC1, PC2, color = factor(cohort)) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  ggrepel::geom_text_repel(
    aes(label = pca.label),
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf,
    segment.alpha = 0.4,
    show.legend = FALSE
  ) +
  theme_bw()
# the variation over PC1/PC2 by diet seems low. The variation seems to be driven by age

# PCA of filtered
mx_filt = dge$counts %>%
  as.matrix() %>%
  t()

pca = prcomp(mx_filt, scale. = TRUE)
pca_var = round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
filtered_pca_var = pca_var
pc_contributions = as.data.frame(pca$x)
pc = pc_contributions %>%
  rownames_to_column(var="sample") %>%
  dplyr::select(sample, PC1, PC2, PC3, PC4, PC5, PC6, PC7,PC8) %>%
  left_join(., meta) %>%
  mutate(
    pca.label = paste0(sample, "\n", age.grp, " ", diet, ", C", cohort),
    pca_set = paste0(
      "Filtered counts\n",
      "PC1 = ", pca_var[1], "%, PC2 = ", pca_var[2], "%"
    )
  )
pca_filtered = pc
pc %>%
  ggplot(.,aes(PC1, PC2, color = diet) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  ggrepel::geom_text_repel(
    aes(label = pca.label),
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf,
    segment.alpha = 0.4,
    show.legend = FALSE
  ) +
  theme_bw()
pc %>%
  ggplot(.,aes(PC1, PC2, color = age.grp) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  ggrepel::geom_text_repel(
    aes(label = pca.label),
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf,
    segment.alpha = 0.4,
    show.legend = FALSE
  ) +
  theme_bw()
pc %>%
  ggplot(.,aes(PC1, PC2, color = factor(cohort)) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  ggrepel::geom_text_repel(
    aes(label = pca.label),
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf,
    segment.alpha = 0.4,
    show.legend = FALSE
  ) +
  theme_bw()
# not substantially different from unfiltered

pca_plot_data = bind_rows(pca_unfiltered, pca_filtered) %>%
  mutate(
    pca_set = factor(pca_set, levels = unique(pca_set)),
    age.grp = factor(age.grp, levels = c("young", "old")),
    diet = factor(diet, levels = c("CD", "ID")),
    cohort = factor(cohort)
  )

rnaseq_pca_plot = ggplot(
  pca_plot_data,
  aes(PC1, PC2, color = age.grp, shape = diet)
) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  stat_ellipse(
    aes(group = interaction(age.grp, diet), linetype = diet),
    type = "norm",
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(size = 3.1, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = paste0(sample, "\nC", cohort)),
    color = "grey20",
    size = 2.9,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    segment.alpha = 0.45,
    show.legend = FALSE
  ) +
  facet_wrap(~ pca_set, scales = "free", nrow = 1) +
  scale_color_manual(values = c("young" = "#2C7FB8", "old" = "#D95F02"), name = "Age") +
  scale_shape_manual(values = c("CD" = 16, "ID" = 17), name = "Diet") +
  scale_linetype_manual(values = c("CD" = "solid", "ID" = "22"), name = "Diet ellipse") +
  labs(
    title = "RNA-seq sample PCA",
    subtitle = "Unfiltered and filtered gene-count matrices",
    x = "PC1",
    y = "PC2"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

rnaseq_pca_plot

dir.create(here("plots/qc"), recursive = TRUE, showWarnings = FALSE)
ggsave(
  here("plots/qc/rnaseq pca unfiltered vs filtered counts.pdf"),
  rnaseq_pca_plot,
  width = 11,
  height = 5.8
)

filtered_group_pca_data = pca_filtered %>%
  mutate(
    analysis_group = case_when(
      age.grp == "young" & diet == "CD" ~ "Young CD",
      age.grp == "old" & diet == "CD" ~ "Old CD",
      age.grp == "old" & diet == "ID" ~ "Old ID",
      TRUE ~ paste(age.grp, diet)
    ),
    analysis_group = factor(analysis_group, levels = c("Young CD", "Old CD", "Old ID")),
    cohort = factor(cohort)
  )

filtered_group_pca_colors = c(
  "Young CD" = "#2C7FB8",
  "Old CD" = "#4D9221",
  "Old ID" = "#D95F02"
)

filtered_group_pca_plot = ggplot(
  filtered_group_pca_data,
  aes(PC1, PC2, color = analysis_group)
) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  stat_ellipse(type = "norm", linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 3.3, alpha = 0.95) +
  ggrepel::geom_text_repel(
    aes(label = paste0(sample, "\nC", cohort)),
    color = "grey20",
    size = 3,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    segment.alpha = 0.45,
    show.legend = FALSE
  ) +
  scale_color_manual(values = filtered_group_pca_colors, name = NULL) +
  labs(
    title = "RNA-seq PCA: filtered counts",
    subtitle = "Major study groups with normal ellipses",
    x = paste0("PC1 (", filtered_pca_var[1], "%)"),
    y = paste0("PC2 (", filtered_pca_var[2], "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

filtered_group_pca_plot

ggsave(
  here("plots/qc/rnaseq pca filtered counts.pdf"),
  filtered_group_pca_plot,
  width = 7.6,
  height = 6.2
)


# ========== 2.0 - SVA for diet contrast ==========
library(edgeR)
library(limma)
library(sva)

# SVA estimates hidden expression-wide structure, such as unrecorded batch or
# sample quality effects. Because there are no young ID samples, estimate SVs
# only within old mice, where the ID vs CD comparison is actually balanced
# across cohorts 1 and 2.
meta_sva <- meta %>%
  dplyr::filter(age.grp == "old") %>%
  mutate(
    cohort = factor(cohort),
    diet = factor(diet, levels = c("CD", "ID"))
  ) %>%
  arrange(cohort, diet) %>%
  column_to_rownames(var = "sample")

dge_sva <- dge[, rownames(meta_sva), keep.lib.sizes = FALSE]
dge_sva <- calcNormFactors(dge_sva)

# Preserve the known biological/model structure while estimating SVs.
# Full model: cohort + diet, where diet is the effect of interest.
# Null model: cohort only, so SVA can estimate hidden variation without
# intentionally removing the diet contrast.
mod_sva  <- model.matrix(~ cohort + diet, data = meta_sva)
mod0_sva <- model.matrix(~ cohort, data = meta_sva)

k_est <- num.sv(dge_sva$counts, mod_sva, method = "leek")

# With 13 old samples, using every estimated SV can overfit the model. Keep this
# conservative for the primary exploratory comparison; increase only after
# inspecting SV correlations and residual degrees of freedom.
max_sv_to_use <- 2
min_residual_df <- 4
k_safe <- max(0, min(
  k_est,
  max_sv_to_use,
  nrow(meta_sva) - qr(mod_sva)$rank - min_residual_df
))
cat("\nSVA estimated ", k_est, " surrogate variables; using ", k_safe, ".\n", sep = "")

if (k_safe > 0) {
  svobj <- svaseq(dat = dge_sva$counts, mod = mod_sva, mod0 = mod0_sva, n.sv = k_safe)
  sv_df <- as.data.frame(svobj$sv)
  colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))
  rownames(sv_df) <- rownames(meta_sva)
  
  # Inspect whether SVs are tracking known covariates. If an SV is strongly
  # correlated with diet, treat the adjusted DE results cautiously.
  pairs(sv_df)
  sv_covariate_cor <- sapply(sv_df, function(z) {
    c(
      cohort = cor(z, as.numeric(meta_sva$cohort)),
      diet = cor(z, as.numeric(meta_sva$diet))
    )
  })
  sv_covariate_cor
  
} else {
  svobj <- NULL
  sv_df <- data.frame(row.names = rownames(meta_sva))
}

meta_with_SV <- cbind(meta_sva, sv_df)

# Build the final DE model dynamically so the code works whether SVA estimates
# zero, one, two, or more surrogate variables.
sv_terms <- colnames(sv_df)
if (length(sv_terms) > 0) {
  sva_formula <- as.formula(paste("~ cohort + diet +", paste(sv_terms, collapse = " + ")))
} else {
  sva_formula <- ~ cohort + diet
}

design_sva <- model.matrix(sva_formula, data = meta_with_SV)

# voom + limma fit for old ID vs old CD, adjusted for cohort and any SVs.
v_sva <- voom(dge_sva, design_sva, plot = FALSE)
fit_sva <- lmFit(v_sva, design_sva)
fit_sva <- eBayes(fit_sva)

dex_sva <- topTable(fit_sva, coef = "dietID", number = Inf) %>%
  rownames_to_column(var = "gene.id")
dex_sva_signif <- dex_sva %>% dplyr::filter(P.Value <= 0.05)

# Compare the diet effect with and without SVA. Both models use the same old
# samples and adjust for known cohort; the only difference is whether SVs are
# included as additional covariates.
design_no_sva <- model.matrix(~ cohort + diet, data = meta_sva)
v_no_sva <- voom(dge_sva, design_no_sva, plot = FALSE)
fit_no_sva <- lmFit(v_no_sva, design_no_sva)
fit_no_sva <- eBayes(fit_no_sva)

dex_no_sva <- topTable(fit_no_sva, coef = "dietID", number = Inf) %>%
  rownames_to_column(var = "gene.id")
dex_no_sva_signif <- dex_no_sva %>% dplyr::filter(P.Value <= 0.05)

sva_diet_comparison <- dex_no_sva %>%
  dplyr::select(
    gene.id,
    logFC_no_sva = logFC,
    P.Value_no_sva = P.Value
  ) %>%
  inner_join(
    dex_sva %>%
      dplyr::select(
        gene.id,
        logFC_sva = logFC,
        P.Value_sva = P.Value
      ),
    by = "gene.id"
  ) %>%
  mutate(
    sig_no_sva = P.Value_no_sva <= 0.05,
    sig_sva = P.Value_sva <= 0.05,
    sig_status = case_when(
      sig_no_sva & sig_sva ~ "Significant in both",
      sig_no_sva & !sig_sva ~ "No SVA only",
      !sig_no_sva & sig_sva ~ "SVA only",
      TRUE ~ "Not significant"
    )
  ) %>%
  dplyr::select(
    gene.id,
    logFC_no_sva,
    logFC_sva,
    P.Value_no_sva,
    P.Value_sva,
    sig_status
  )

sva_diet_summary <- sva_diet_comparison %>%
  count(sig_status)

sva_diet_summary

genes_sva_only <- sva_diet_comparison %>%
  dplyr::filter(sig_status == "SVA only") %>%
  arrange(P.Value_sva, desc(abs(logFC_sva)))

genes_no_sva_only <- sva_diet_comparison %>%
  dplyr::filter(sig_status == "No SVA only") %>%
  arrange(P.Value_no_sva, desc(abs(logFC_no_sva)))

genes_sva_only
genes_no_sva_only

# ----- GO BP ORA on method-specific genes ---
# Split unique genes by direction before ORA. This avoids mixing pathways from
# genes that increase in ID with pathways from genes that decrease in ID.
library(enrichR)
go_bp_db <- "GO_Biological_Process_2025"

genes_sva_only_up <- genes_sva_only %>%
  dplyr::filter(logFC_sva > 0)

genes_sva_only_down <- genes_sva_only %>%
  dplyr::filter(logFC_sva < 0)

genes_no_sva_only_up <- genes_no_sva_only %>%
  dplyr::filter(logFC_no_sva > 0)

genes_no_sva_only_down <- genes_no_sva_only %>%
  dplyr::filter(logFC_no_sva < 0)

genes_sva_only_up_vec <- genes_sva_only_up %>%
  pull(gene.id)

genes_sva_only_down_vec <- genes_sva_only_down %>%
  pull(gene.id)

genes_no_sva_only_up_vec <- genes_no_sva_only_up %>%
  pull(gene.id)

genes_no_sva_only_down_vec <- genes_no_sva_only_down %>%
  pull(gene.id)

# ----- GO BP ORA: SVA-only up genes ---
if (length(genes_sva_only_up_vec) > 0) {
  go_bp_sva_only_up <- enrichr(genes_sva_only_up_vec, databases = go_bp_db)
  go_bp_sva_only_up_tbl <- go_bp_sva_only_up[[go_bp_db]] %>%
    dplyr::filter(P.value <= 0.05)
} else {
  go_bp_sva_only_up <- NULL
  go_bp_sva_only_up_tbl <- tibble()
}

# ----- GO BP ORA: SVA-only down genes ---
if (length(genes_sva_only_down_vec) > 0) {
  go_bp_sva_only_down <- enrichr(genes_sva_only_down_vec, databases = go_bp_db)
  go_bp_sva_only_down_tbl <- go_bp_sva_only_down[[go_bp_db]] %>%
    dplyr::filter(P.value <= 0.05)
} else {
  go_bp_sva_only_down <- NULL
  go_bp_sva_only_down_tbl <- tibble()
}

# ----- GO BP ORA: no-SVA-only up genes ---
if (length(genes_no_sva_only_up_vec) > 0) {
  go_bp_no_sva_only_up <- enrichr(genes_no_sva_only_up_vec, databases = go_bp_db)
  go_bp_no_sva_only_up_tbl <- go_bp_no_sva_only_up[[go_bp_db]] %>%
    dplyr::filter(P.value <= 0.05)
} else {
  go_bp_no_sva_only_up <- NULL
  go_bp_no_sva_only_up_tbl <- tibble()
}

# ----- GO BP ORA: no-SVA-only down genes ---
if (length(genes_no_sva_only_down_vec) > 0) {
  go_bp_no_sva_only_down <- enrichr(genes_no_sva_only_down_vec, databases = go_bp_db)
  go_bp_no_sva_only_down_tbl <- go_bp_no_sva_only_down[[go_bp_db]] %>%
    dplyr::filter(P.value <= 0.05)
} else {
  go_bp_no_sva_only_down <- NULL
  go_bp_no_sva_only_down_tbl <- tibble()
}

go_bp_sva_only_up_tbl
go_bp_sva_only_down_tbl
go_bp_no_sva_only_up_tbl
go_bp_no_sva_only_down_tbl

# ----- Plot GO BP ORA results as a compact dot plot ---
# Dot size is the gene overlap count, and the x-axis/color show enrichment
# strength as -log10(P-value). Facets keep the four method/direction groups
# separate without compressing bars.
go_bp_unique_tbl <- bind_rows(
  go_bp_sva_only_up_tbl %>%
    mutate(comparison = "SVA only up"),
  go_bp_sva_only_down_tbl %>%
    mutate(comparison = "SVA only down"),
  go_bp_no_sva_only_up_tbl %>%
    mutate(comparison = "No SVA only up"),
  go_bp_no_sva_only_down_tbl %>%
    mutate(comparison = "No SVA only down")
)

if (nrow(go_bp_unique_tbl) > 0) {
  go_bp_unique_dotplot_data <- go_bp_unique_tbl %>%
    group_by(comparison) %>%
    slice_min(P.value, n = 8, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      comparison = factor(
        comparison,
        levels = c(
          "No SVA only up",
          "SVA only up",
          "No SVA only down",
          "SVA only down"
        )
      ),
      overlap_count = as.numeric(str_extract(Overlap, "^[0-9]+")),
      neg_log10_p = -log10(P.value),
      term_wrapped = str_wrap(Term, width = 55)
    ) %>%
    arrange(comparison, P.value) %>%
    mutate(
      term_for_plot = paste(term_wrapped, comparison, sep = "___"),
      term_for_plot = factor(term_for_plot, levels = rev(unique(term_for_plot)))
    )
  
  go_bp_unique_plot <- ggplot(
    go_bp_unique_dotplot_data,
    aes(neg_log10_p, term_for_plot)
  ) +
    geom_point(aes(size = overlap_count, color = neg_log10_p), alpha = 0.85) +
    facet_wrap(~comparison, scales = "free_y", ncol = 2) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_color_gradient(low = "#4575B4", high = "#D73027", name = "-log10(P)") +
    labs(
      title = "GO Biological Process ORA for SVA/no-SVA Unique Genes",
      x = "-log10(P-value)",
      y = "",
      size = "Gene count"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )
} else {
  go_bp_unique_plot <- ggdraw() +
    draw_label("No GO BP ORA plots available for unique SVA/no-SVA gene sets", size = 12)
}

go_bp_unique_plot

fc_sva_comparison_plot <- ggplot(sva_diet_comparison, aes(logFC_no_sva, logFC_sva, color = sig_status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.5, size = 0.8) +
  labs(
    title = "Diet Effect With vs Without SVA",
    x = "Log2FC ID/CD without SVA",
    y = "Log2FC ID/CD with SVA",
    color = ""
  ) +
  theme_bw()

pvalue_axis_max <- ceiling(max(
  -log10(sva_diet_comparison$P.Value_no_sva),
  -log10(sva_diet_comparison$P.Value_sva),
  na.rm = TRUE
))

pvalue_sva_comparison_plot <- ggplot(sva_diet_comparison, aes(-log10(P.Value_no_sva), -log10(P.Value_sva), color = sig_status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.5, size = 0.8) +
  coord_equal(xlim = c(0, pvalue_axis_max), ylim = c(0, pvalue_axis_max)) +
  labs(
    title = "Diet P-values With vs Without SVA",
    x = "-log10(P-value) without SVA",
    y = "-log10(P-value) with SVA",
    color = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    aspect.ratio = 1
  )

fc_sva_comparison_plot
pvalue_sva_comparison_plot

# Volcano plots for the same old ID vs CD diet effect.
# Green genes are significant in both approaches. Blue genes are significant
# only in the approach shown in that panel. Everything else is grey.
volcano_no_sva_data <- sva_diet_comparison %>%
  mutate(
    volcano_group = case_when(
      sig_status == "Significant in both" ~ "Both",
      sig_status == "No SVA only" ~ "Contrasting",
      TRUE ~ "Other"
    )
  )

volcano_sva_data <- sva_diet_comparison %>%
  mutate(
    volcano_group = case_when(
      sig_status == "Significant in both" ~ "Both",
      sig_status == "SVA only" ~ "Contrasting",
      TRUE ~ "Other"
    )
  )

volcano_colors <- c(
  "Both" = "#1B9E77",
  "Contrasting" = "#D73027",
  "Other" = "grey80"
)

volcano_no_sva_labels <- volcano_no_sva_data %>%
  dplyr::filter(volcano_group == "Contrasting") %>%
  slice_min(P.Value_no_sva, n = 60, with_ties = FALSE)

volcano_sva_labels <- volcano_sva_data %>%
  dplyr::filter(volcano_group == "Contrasting") %>%
  slice_min(P.Value_sva, n = 60, with_ties = FALSE)

volcano_no_sva <- ggplot(volcano_no_sva_data, aes(logFC_no_sva, -log10(P.Value_no_sva))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(aes(color = volcano_group), alpha = 0.55, size = 0.7) +
  ggrepel::geom_text_repel(
    data = volcano_no_sva_labels,
    aes(label = gene.id),
    size = 2.7,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.35,
    show.legend = FALSE
  ) +
  scale_color_manual(values = volcano_colors, name = "") +
  labs(
    title = "No SVA",
    x = "Log2FC ID/CD",
    y = "-log10(P-value)"
  ) +
  theme_bw() +
  theme(
    legend.position = "top"
  )

volcano_sva <- ggplot(volcano_sva_data, aes(logFC_sva, -log10(P.Value_sva))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(aes(color = volcano_group), alpha = 0.55, size = 0.7) +
  ggrepel::geom_text_repel(
    data = volcano_sva_labels,
    aes(label = gene.id),
    size = 2.7,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.35,
    show.legend = FALSE
  ) +
  scale_color_manual(values = volcano_colors, name = "") +
  labs(
    title = "SVA",
    x = "Log2FC ID/CD",
    y = "-log10(P-value)"
  ) +
  theme_bw() +
  theme(
    legend.position = "top"
  )

volcano_sva_comparison_plot <- plot_grid(volcano_no_sva, volcano_sva, nrow = 1)
volcano_sva_comparison_plot

# ----- Save SVA vs no-SVA summary figure ---
# This figure combines the p-value comparison, method-specific volcano plots,
# and GO BP ORA plots for genes unique to each method/direction.
sva_no_sva_summary_figure <- plot_grid(
  plot_grid(NULL, pvalue_sva_comparison_plot, NULL, nrow = 1, rel_widths = c(0.7, 1, 0.7)),
  volcano_sva_comparison_plot,
  go_bp_unique_plot,
  ncol = 1,
  rel_heights = c(1.15, 1.2, 1.8)
)

sva_no_sva_summary_figure

dir.create(here("plots/sva"), recursive = TRUE, showWarnings = FALSE)
ggsave(
  here("plots/sva/sva vs no sva comparison pvalues volcano ora.pdf"),
  sva_no_sva_summary_figure,
  height = 16,
  width = 12
)

# Carry the selected primary result forward for downstream interactive code.
# Proceed with the cohort-adjusted no-SVA model for old ID vs CD.
dex <- dex_no_sva
dex_signif <- dex_no_sva_signif





# ========== 3.0 - DESeq2 ~cohort + diet (old only) ==========
# --
# Proceeding with filtered dge object; not applying any SVA adjustment
meta = read.csv(here("data/metadata/metadata.csv")) %>%
  mutate(sample = paste0("r",sample)) %>%
  arrange(diet, age.grp) %>%
  column_to_rownames(var = "sample")

filt = dge$counts
length(unique(rownames(filt)))
# genes after edgeR expression and variance filtering

library(DESeq2)

samples = meta %>%
  dplyr::filter(age.grp == "old") %>%
  mutate(
    cohort = factor(cohort),
    diet = factor(diet, levels = c("CD", "ID"))
  ) %>%
  arrange(cohort, diet)

x = filt[, rownames(samples)]
x[is.na(x)] = 0
x = round(as.matrix(x))
storage.mode(x) = "integer"
stopifnot(identical(colnames(x), rownames(samples)))

cat("\nUsing old samples with design formula ~ cohort + diet\n")
dds <- DESeqDataSetFromMatrix(countData = x, colData = samples, design = ~ cohort + diet)

dds <- DESeq(dds)
res <- results(dds, contrast = c("diet", "ID", "CD"))

res2 = data.frame(res) %>%
  arrange(padj, pvalue) %>%
  rownames_to_column(var = "gene.id") %>%
  mutate(
    design = "~cohort + diet",
    comparison = "old ID vs CD"
  )

tmp2 = res2 %>% dplyr::filter(pvalue <= 0.05)
dir.create(here("data/processed/deseq2"), recursive = TRUE, showWarnings = FALSE)
write.csv(res2, here("data/processed/deseq2/deseq2 ID vs CD old only cohort adjusted.csv"), row.names = FALSE)

dex2 = res2 %>%
  dplyr::filter(pvalue <= 0.05) %>%
  mutate(score = -log10(pvalue)*abs(log2FoldChange)) %>%
  arrange(desc(score))

unique(dex2$gene.id)

saveRDS(
  list(
    meta = meta,
    dge_filtered = dge,
    filtered_counts = filt,
    pca_unfiltered = pca_unfiltered,
    pca_filtered = pca_filtered,
    filtered_pca_var = filtered_pca_var,
    pca_plot_data = pca_plot_data
  ),
  here("data/processed/counts/rnaseq filtered counts qc objects.rds")
)

dir.create(here("data/processed/sva"), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    k_est = k_est,
    k_safe = k_safe,
    sv_df = sv_df,
    meta_with_SV = meta_with_SV,
    dex_sva = dex_sva,
    dex_no_sva = dex_no_sva,
    dex_sva_signif = dex_sva_signif,
    dex_no_sva_signif = dex_no_sva_signif,
    sva_diet_comparison = sva_diet_comparison,
    sva_diet_summary = sva_diet_summary
  ),
  here("data/processed/sva/sva old ID vs CD sensitivity objects.rds")
)

saveRDS(
  list(
    dds = dds,
    res = res,
    res2 = res2,
    dex2 = dex2,
    samples = samples,
    old_filtered_counts = x,
    design = "~ cohort + diet",
    comparison = "old ID vs CD"
  ),
  here("data/processed/deseq2/deseq2 ID vs CD old only cohort adjusted objects.rds")
)
