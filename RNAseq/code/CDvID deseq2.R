
library(here)
library(tidyverse)
library(edgeR)
library(cowplot)

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
write.csv(cnts,here("data/processed/counts matrix with mouse gene symbols.csv"),row.names = FALSE)





# ========== 1.0 - Filter and PCA ==========
meta = read.csv(here("data/metadata/metadata.csv")) %>%
  mutate(sample = paste0("r",sample)) %>%
  arrange(diet, age.grp) %>%
  column_to_rownames(var = "sample")

x = read.csv(here("data/processed/counts matrix with mouse gene symbols.csv"))
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
pc_contributions = as.data.frame(pca$x)
pc = pc_contributions %>%
  rownames_to_column(var="sample") %>%
  dplyr::select(sample, PC1, PC2, PC3, PC4, PC5, PC6, PC7,PC8) %>%
  left_join(., meta) %>%
  mutate(
    pca.label = paste0(sample, "\n", age.grp, " ", diet, ", C", cohort)
  )

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
pc_contributions = as.data.frame(pca$x)
pc = pc_contributions %>%
  rownames_to_column(var="sample") %>%
  dplyr::select(sample, PC1, PC2, PC3, PC4, PC5, PC6, PC7,PC8) %>%
  left_join(., meta) %>%
  mutate(
    pca.label = paste0(sample, "\n", age.grp, " ", diet, ", C", cohort)
  )
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


# ========== 2.0 - SVA ==========
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

# ----- GO BP ORA on method-specific genes -----
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

# ----- GO BP ORA: SVA-only up genes -----
if (length(genes_sva_only_up_vec) > 0) {
  go_bp_sva_only_up <- enrichr(genes_sva_only_up_vec, databases = go_bp_db)
  go_bp_sva_only_up_tbl <- go_bp_sva_only_up[[go_bp_db]] %>%
    dplyr::filter(P.value <= 0.05)
} else {
  go_bp_sva_only_up <- NULL
  go_bp_sva_only_up_tbl <- tibble()
}

# ----- GO BP ORA: SVA-only down genes -----
if (length(genes_sva_only_down_vec) > 0) {
  go_bp_sva_only_down <- enrichr(genes_sva_only_down_vec, databases = go_bp_db)
  go_bp_sva_only_down_tbl <- go_bp_sva_only_down[[go_bp_db]] %>%
    dplyr::filter(P.value <= 0.05)
} else {
  go_bp_sva_only_down <- NULL
  go_bp_sva_only_down_tbl <- tibble()
}

# ----- GO BP ORA: no-SVA-only up genes -----
if (length(genes_no_sva_only_up_vec) > 0) {
  go_bp_no_sva_only_up <- enrichr(genes_no_sva_only_up_vec, databases = go_bp_db)
  go_bp_no_sva_only_up_tbl <- go_bp_no_sva_only_up[[go_bp_db]] %>%
    dplyr::filter(P.value <= 0.05)
} else {
  go_bp_no_sva_only_up <- NULL
  go_bp_no_sva_only_up_tbl <- tibble()
}

# ----- GO BP ORA: no-SVA-only down genes -----
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

# ----- Plot GO BP ORA results as a compact dot plot -----
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

# ----- Save SVA vs no-SVA summary figure -----
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

ggsave(
  here("plots/sva vs no sva comparison pvalues volcano ora.pdf"),
  sva_no_sva_summary_figure,
  height = 16,
  width = 12
)

# Carry the selected primary result forward for downstream interactive code.
# Proceed with the cohort-adjusted no-SVA model for old ID vs CD.
dex <- dex_no_sva
dex_signif <- dex_no_sva_signif





# ========== 3.0 - Deseq2 ==========
# --

meta = read.csv(here("data/metadata.csv")) %>%
  mutate(sample = paste0("r",sample)) %>%
  arrange(diet, age.grp) %>%
  column_to_rownames(var = "sample")

filt = dge$counts
length(unique(rownames(filt)))
# 12768 genes after filtering

library(DESeq2)
run_deseq = function(age = c("old","young"), numerator, denominator){
  
  cat(paste0("...\n...\nStarting ", numerator, "/", denominator))
  
  # denominator = "CD"
  # numerator   = "ID"
  # age         = c("old", "young")
  
  samples = meta %>%
    filter(age.grp %in% age) %>%
    arrange(diet, cohort)
  
  x = cnts[, rownames(samples)]
  x[is.na(x)] = 0
  
  # -------------------------------
  # 1. Run DESeq2 with diet as sole design factor
  # -------------------------------
  cat("\nUsing design formula ~diet\n")
  ddsA <- DESeqDataSetFromMatrix(countData = x, colData = samples, design = ~ diet)
  ddsA2 <- DESeq(ddsA)
  resA  <- results(ddsA2, contrast = c("diet", numerator, denominator))
  
  lfcA = data.frame(resA) %>%
    arrange(-log2FoldChange) %>%
    rownames_to_column(var = "gene.id") %>%
    mutate(design = "~diet")
  
  # -------------------------------
  # 2. If multiple age groups exist, add age as a covariate
  # -------------------------------
  if(length(unique(samples$age.grp)) > 1){
    cat("\nUsing design formula ~age.grp + diet\n")
    ddsB  <- DESeqDataSetFromMatrix(countData = x, colData = samples, design = ~ age.grp + diet)
    ddsB2 <- DESeq(ddsB)
    resB  <- results(ddsB2, contrast = c("diet", numerator, denominator))
    
    lfcB = data.frame(resB) %>%
      arrange(-log2FoldChange) %>%
      rownames_to_column(var = "gene.id") %>%
      mutate(design = "~age.grp + diet")
    
    lfc = rbind(lfcA, lfcB)
    
  } else {
    lfc = lfcA
  }
  
  cat(paste0("\nCompleted ", numerator, "/", denominator, " !\n"))
  return(lfc)
}


res1 = run_deseq(age = c("old","young"), "ID", "CD")
tmp1 = res1 %>% filter(padj<=0.05)
#write.csv(res1, here("data/processed/deseq2 ID vs CD old and young design w and wo age grp covariate.csv"))


res2 = run_deseq(age = c("old"), "ID", "CD")
tmp2 = res2 %>% filter(pvalue<=0.05)
# write.csv(res2, here("data/processed/deseq2 ID vs CD old only.csv"))

dex1 = res1 %>%
  filter(pvalue<=0.05, design == "~age.grp + diet") %>%
  mutate(score = -log10(pvalue)*abs(log2FoldChange)) %>%
  arrange(desc(score))

dex2 = res2 %>%
  filter(pvalue<=0.05) %>%
  mutate(score = -log10(pvalue)*abs(log2FoldChange)) %>%
  arrange(desc(score))
  
unique(dex$gene.id)


# ========== preliminary Volcano ==========

# -- Some aging genes I got from GPT
geneset = c(
  # Epigenetic aging / methylation markers
  "Elovl2", "Hsf4", "Prima1", "Aspa", 
  # Sirtuins & chromatin modifiers
  "Sirt1", "Sirt3", "Sirt6", "Sirt7", "Kat7", 
  # DNA repair / genome stability
  "Xrcc6", "Xrcc5", "Trp53", "Lmna", "Tert", "Wrn", 
  # Senescence / cell cycle arrest
  "Cdkn2a", "Cdkn1a", "Rb1", 
  # Nutrient sensing / metabolism
  "Mtor", "Rps6kb1", "Igf1r", "Irs1", "Akt1", "Pdk1", "Foxo1", "Foxo3", 
  # Stress response / proteostasis
  "Hspa1a", "Hspa1b", "Hsp90aa1", "Canx", "Ddost", 
  # Mitochondrial / ROS
  "Prkn", "Pink1", "Sod1", "Sod2", "Gpx1", 
  # Muscle & extracellular matrix / aging tissue function
  "Wnt3a", "Hspg2", "Fgd6", "Apod", "Gprc5b", "Tpp1"
)

res = read.csv(here("data/processed/deseq2 ID vs CD old and young design w and wo age grp covariate.csv")) %>%
  filter(
    design == "~age.grp + diet"
  ) %>%
  mutate(X = NULL)

lfc = data.frame(res) %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(sig.label = ifelse(pvalue <=0.05, gene.id, NA))

vol = ggplot(lfc, aes(log2FoldChange, -log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black", alpha = .5) +
  geom_vline(xintercept = -1, linetype = "dotted", color = "black", alpha = .5) +
  geom_point(size = .2, alpha = .4, color = "grey80") +
  geom_point(size = .2,
             data = subset(lfc, is.na(sig.label) == FALSE),
             color = "black"
  ) +
  labs( x="Log2FC(ID/CD)", "Muscle Gene Expression Inulin Diet vs Control") +
  xlim(-7,7) + ylim(0,6)+
  ggrepel::geom_text_repel( aes( label = sig.label ),
                            vjust = 1.0,
                            box.padding = 0.5,
                            size = 2.0,
                            max.overlaps = 50, alpha = 0.7, segment.alpha = .2 ) +
  theme_bw()





# ========== Pathway Enrichment ==========

library(enrichR)
dbs = as.data.frame(listEnrichrDbs())
mydb = "MSigDB_Hallmark_2020"
mydb = "GO_Cellular_Component_2025" #Collagen-Containing Extracellular Matrix (GO:0062023)


up = subset(lfc, log2FoldChange > 0 & pvalue <= 0.05 )$gene.id
enriched_up = enrichr(up,databases = mydb)
df1 = as.data.frame(enriched)
colnames(df1) = sub(".*?\\.", "", colnames(df1))
df1 = subset(df1, P.value <= 0.05)
enr1= plotEnrich(enriched_up[[1]], showTerms = 20, numChar = 100, y = "Count", orderBy = "P.value", title="Upregulated Pathways in ID")

dn = subset(lfc, log2FoldChange < 0 & pvalue <= 0.05 )$gene.id
enriched_dn = enrichr(dn,databases = mydb)
df1 = as.data.frame(enriched)
colnames(df1) = sub(".*?\\.", "", colnames(df1))
df1 = subset(df1, P.value <= 0.05)
enr2= plotEnrich(enriched_dn[[1]], showTerms = 20, numChar = 100, y = "Count", orderBy = "P.value", title="Downregulated Pathways in ID")


library(cowplot)
enr = plot_grid(enr1, enr2, nrow=1)
plot_grid(vol, enr, ncol=1, rel_heights = c(1.5,1))


# ========== Plot fig for only ID vs CD ==========
# --

res = read.csv(here("data/processed/deseq2 ID vs CD old only.csv")) %>%
  mutate(X = NULL)

lfc = data.frame(res) %>%
  filter(!is.na(log2FoldChange))

# ========== (A) Volcano ====

# categories
# -- A. Macrophage / immune --
macrophage_genes <- c(
  "Il6", "Csf3r", "Ccl8", "Chil3", "Retnlg", "Slfn1",
  "Gbp4", "Gbp8", "Tgtp1", "Tgtp2", "Ifi204", "Psmb8", "Sp110",
  "Csf1r", "Cd68", "Aif1", "Cd33", "Fcgr3", "Mrc1", "Clec7a",
  "Clec4n", "Siglec1", "Slamf8", "Cd209f", "Stab1", "Ackr1",
  "Cx3cr1", "Socs2", "Cish", "Nlrp1a", "Susd1", "Tap2",
  "Gpr141", "H2-T10", "H2-Q5", "H2-Q7", "H2-M3", "Gbp9"
)

# -- B. ECM / fibrosis --
ecm_fibrosis_genes <- c(
  "Col2a1", "Col6a6", "Col14a1", "Frem2", "Loxl3", "F13a1",
  "Adamts10", "Adam33", "Agrn", "Cspg4", "Pdgfrl",
  "Itga8", "Itgb4", "Pkd1", "Pi16", "Tgfbi", "Ltbp4", "Tgfbr2"
)

# -- C. Vascular / endothelial / pericyte --
vascular_genes <- c(
  "Vwf", "Ackr1", "Cspg4", "Notch1", "Acvr2b", "Pdgfrl",
  "Stab1", "Ace"
)

# -- D. Neural / glial / nerve-associated --
neural_genes <- c(
  "Gfap", "Gad1", "Nell2", "Slitrk1", "Kcnf1", "Kcna1",
  "Kcna2", "Atp1a3", "Cngb1", "Pcdh20", "Pcdhb18", "Pcdhb16",
  "Pcdhga10", "Ptprz1", "Gpm6a", "Galr2", "Gpr85",
  "Sim1", "Efnb3", "Gprin2"
)

# -- E. Metabolic / mitochondrial / oxidative --
metabolic_genes <- c(
  "Cox10", "mt-Cytb", "mt-Nd2", "Pecr", "Acox2", "Aldh3a2",
  "Aass", "Gck", "Pemt", "Lcat", "Pon1", "Msra",
  "Slc7a11", "Lsr"
)

# (Optional) secreted / liver-like proteins – could be merged into metabolic if you prefer
secreted_genes <- c(
  "Apoc1", "Hp", "Fetub", "Apoh", "Igfbp2", "Proz"
)

voldata = lfc %>%
  mutate(
    category = case_when(
      gene.id %in% macrophage_genes ~ "Macrophage / immune",
      gene.id %in% ecm_fibrosis_genes ~ "ECM / fibrosis",
      gene.id %in% vascular_genes ~ "Vascular / endothelial",
      gene.id %in% neural_genes ~ "Neural / glial / nerve",
      gene.id %in% c(metabolic_genes, secreted_genes) |
        str_detect(gene.id, "^mt-") ~ "Metabolic / mitochondrial",
      TRUE ~ "Other"
    ),
  )

# Define color palette (high-contrast, colorblind-safe)
cat_colors <- c(
  "Macrophage / immune"       = "#D73027",  # red
  "ECM / fibrosis"            = "#4575B4",  # blue
  "Vascular / endothelial"    = "#1A9850",  # green
  "Neural / glial / nerve"    = "#E66101",  # gold
  "Metabolic / mitochondrial" = "#984EA3",  # purple
  "Other"                     = "grey80"
)

# Ensure factor order
voldata$category <- factor(
  voldata$category,
  levels = names(cat_colors)
)

vol_cat = ggplot(voldata, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", alpha = .5) +
  geom_point(
    data = subset(voldata, category == "Other"),
    color = "grey80", size = 0.3, alpha = 0.3
  ) +
  geom_point(
    data = subset(voldata, category != "Other"),
    aes(color = category),
    size = 1, alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    data = subset(voldata, !is.na(sig.label) & category != "Other"),
    aes(label = sig.label, color = category),
    size = 4,
    box.padding = 0.6,
    max.overlaps = 60,
    segment.alpha = 0.25,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cat_colors, name = "Category") +
  labs(
    x = "Log2FC (ID/CD)",
    y = "-log10(P-value)",
    title = "Muscle DEx Inulin vs Control (Aged Mice)"
  ) +
  xlim(-3.5, 3.5) + ylim(0, 3.5) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  guides(
    col = guide_legend(ncol = 3)
  )



# ========== (B) Bar plots ==========
# --
# Start from your lfc object
df2 <- lfc %>%
  mutate(
    ID = 2^log2FoldChange,  # fold-change vs CD
    CD = 1                  # reference group
  ) %>%
  select(gene.id, ID, CD, pvalue) %>%
  pivot_longer(
    cols      = c("ID", "CD"),
    names_to  = "group",
    values_to = "FC"
  )


plot_gene_category_bar <- function(df2, genes, title = "Gene category", p_thresh = NULL) {
  # Subset to genes of interest
  df_sub <- df2 %>%
    filter(gene.id %in% genes)
  
  # Optional p-value filter for plotting
  if (!is.null(p_thresh)) {
    df_sub <- df_sub %>% filter(pvalue <= p_thresh)
  }
  
  # If nothing left, bail nicely
  if (nrow(df_sub) == 0) {
    warning("No rows to plot after filtering for genes/p_thresh")
    return(NULL)
  }
  
  # Annotation data frame for significance stars
  p_annot <- df_sub %>%
    group_by(gene.id) %>%
    summarise(
      pvalue = min(pvalue, na.rm = TRUE),
      y_pos  = max(FC, na.rm = TRUE) + 0.2
    ) %>%
    mutate(
      sig = case_when(
        pvalue < 0.001 ~ "***",
        pvalue < 0.01  ~ "**",
        pvalue < 0.05  ~ "*",
        TRUE           ~ "ns"
      )
    )
  
  # Order genes nicely on x-axis
  df_sub <- df_sub %>%
    mutate(
      gene.id = factor(gene.id, levels = sort(unique(gene.id)))
    )
  
  # Plot
  gg <- ggplot(df_sub, aes(gene.id, FC, fill = group)) +
    geom_bar(
      stat     = "identity",
      position = position_dodge(width = 0.8),
      width    = 0.7,
      color    = "grey30"
    ) +
    geom_text(
      data        = p_annot,
      aes(x = gene.id, y = y_pos, label = sig),
      inherit.aes = FALSE,
      size        = 4
    ) +
    labs(
      y     = "Fold Change (ID/CD)",
      x     = "",
      title = title
    ) +
    scale_fill_manual(
      values = c("CD" = "#BFBFBF", "ID" = "#A52A2A"),
      name   = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(gg)
}

# -- Collagen genes
collagen_genes <- lfc %>%
  filter(str_detect(gene.id, "^Col")) %>%
  distinct(gene.id) %>%
  pull(gene.id)

gg_collagen <- plot_gene_category_bar(
  df2,
  genes  = collagen_genes,
  title  = "Significant Collagen Genes",
  p_thresh = 0.05         # only plot genes with p ≤ 0.05
)

gg_collagen

# -- Markers from Graber et al 2023
Graber_2023 <- c(
  "Psmb8","Tspo","Irf7","Ms4a6a","Isg15", "Chrna9", "Gadd45a", "Gadd45g",
  "Shroom4", "P2ry1", "Ccl5","Sln","Arhgdig","Cdkn1a","Chrng","Slc23a3","Cxcl10"
)

gg_graber <- plot_gene_category_bar(
  df2,
  genes  = Graber_2023,
  title  = "Graber 2023 Signature",
  p_thresh = NULL         # show all, regardless of p
)

gg_graber

# -- Notable top LFC genes
mygenes <- c("Cyp2c54","Serpine1","Serpina3n","Serpina11","Smim22")

gg_top <- plot_gene_category_bar(
  df2,
  genes  = mygenes,
  title  = "Top LFC Genes",
  p_thresh = 0.05
)

gg_top

# -- Macrophage infiltration
mac_markers <- c("Mrc1","Cd68","S100a9","Il6","Ccl8")

gg_mac <- plot_gene_category_bar(
  df2,
  genes  = mac_markers,
  title  = "Macrophage Marker Genes",
  p_thresh = NULL    # maybe keep all, or set 0.05 if you want only sig.
)

gg_mac


plot_grid(vol_cat, gg_graber, plot_grid(gg_collagen, gg_top, gg_mac, nrow=1),ncol=1, rel_heights = c(1,.6,.6))
ggsave(here("plots/volcano and lfc barplots.pdf"),h=13,w=9)


# ========== (C) Pwy enrichment ==========
# --

library(enrichR)
dbs = as.data.frame(listEnrichrDbs())
mydb = "MSigDB_Hallmark_2020"
#mydb = "GO_Cellular_Component_2025" #Collagen-Containing Extracellular Matrix (GO:0062023)


up = subset(lfc, log2FoldChange > 0 & pvalue <= 0.05 )$gene.id
enriched_up = enrichr(up,databases = mydb)
df1 = as.data.frame(enriched_up)
colnames(df1) = sub(".*?\\.", "", colnames(df1))
df1 = subset(df1, P.value <= 0.05)
enr1= plotEnrich(enriched_up[[1]], showTerms = 20, numChar = 100, y = "Count", orderBy = "P.value", title="Upregulated Pathways in ID")

dn = subset(lfc, log2FoldChange < 0 & pvalue <= 0.05 )$gene.id
enriched_dn = enrichr(dn,databases = mydb)
df1 = as.data.frame(enriched_dn)
colnames(df1) = sub(".*?\\.", "", colnames(df1))
df1 = subset(df1, P.value <= 0.05)
enr2= plotEnrich(enriched_dn[[1]], showTerms = 20, numChar = 100, y = "Count", orderBy = "P.value", title="Downregulated Pathways in ID")


library(cowplot)
plot_grid(enr1, enr2, nrow=1)
