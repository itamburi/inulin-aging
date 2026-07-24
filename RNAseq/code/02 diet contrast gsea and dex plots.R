library(here)
library(tidyverse)
library(cowplot)

here::i_am("code/02 diet contrast gsea and dex plots.R")

deseq2_objects = readRDS(
  here("data/processed/deseq2/deseq2 ID vs CD old only cohort adjusted objects.rds")
)

dds = deseq2_objects$dds
res = deseq2_objects$res
res2 = deseq2_objects$res2
dex2 = deseq2_objects$dex2
samples = deseq2_objects$samples

dir.create(here("data/processed/pathway_enrichment/go_gsea"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("data/processed/pathway_enrichment/msigdb_hallmark"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("plots/diet DEx/gsea"), recursive = TRUE, showWarnings = FALSE)

# ========== 1.0 - Volcano ~cohort + diet (old only) ====

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

lfc = res2 %>%
  dplyr::filter(!is.na(log2FoldChange), !is.na(pvalue)) %>%
  mutate(
    sig.label = ifelse(gene.id %in% dex2$gene.id, gene.id, NA_character_),
    neg_log10_pvalue = -log10(pmax(pvalue, .Machine$double.xmin))
  )

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

# Muted, print-friendly palette for manuscript-style figures.
cat_colors <- c(
  "Macrophage / immune"       = "#B2182B",
  "ECM / fibrosis"            = "#2166AC",
  "Vascular / endothelial"    = "#1B7837",
  "Neural / glial / nerve"    = "#D95F02",
  "Metabolic / mitochondrial" = "#762A83",
  "Other"                     = "#D9D9D9"
)

# Ensure factor order
voldata$category <- factor(
  voldata$category,
  levels = names(cat_colors)
)

volcano_labels = voldata %>%
  dplyr::filter(!is.na(sig.label), category != "Other") %>%
  arrange(pvalue, desc(abs(log2FoldChange))) %>%
  slice_head(n = 20)

x_axis_limit = max(2, quantile(abs(voldata$log2FoldChange), 0.995, na.rm = TRUE))
x_axis_limit = min(3, ceiling(x_axis_limit * 2) / 2)
y_axis_limit = max(2, quantile(voldata$neg_log10_pvalue, 0.995, na.rm = TRUE))
y_axis_limit = min(4, ceiling(y_axis_limit * 2) / 2)

vol_cat = ggplot(voldata, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "22",
    color = "grey55",
    linewidth = 0.35
  ) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.25) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "22",
    color = "grey55",
    linewidth = 0.35
  ) +
  geom_point(
    data = subset(voldata, category == "Other"),
    color = "#D0D0D0",
    size = 1.1,
    alpha = 0.55,
    stroke = 0
  ) +
  geom_point(
    data = subset(voldata, category != "Other"),
    aes(color = category),
    size = 2.2,
    alpha = 0.95,
    stroke = 0
  ) +
  ggrepel::geom_text_repel(
    data = volcano_labels,
    aes(label = sig.label),
    color = "grey15",
    size = 3.6,
    box.padding = 0.45,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.color = "grey35",
    segment.size = 0.45,
    segment.alpha = 0.95,
    segment.linetype = "solid",
    force = 1.5,
    force_pull = 1,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cat_colors, name = NULL, drop = FALSE) +
  labs(
    x = expression(log[2]~fold~change~"(ID/CD)"),
    y = expression(-log[10]~italic(P)),
    title = "Old muscle: inulin diet vs control",
    subtitle = "DESeq2 model: ~ cohort + diet"
  ) +
  coord_cartesian(xlim = c(-x_axis_limit, x_axis_limit), ylim = c(0, y_axis_limit)) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    axis.ticks.length = unit(2, "mm"),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13),
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 12, color = "grey30", hjust = 0),
    plot.margin = margin(6, 8, 6, 6),
    legend.position = "top",
    legend.justification = "left",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.key.width = unit(5, "mm"),
    legend.key.height = unit(4, "mm")
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(size = 3, alpha = 1)
    )
  )

vol_cat



# ========== 2.0 - GO BP GSEA and ORA ==========
# Uses the old ID vs CD DESeq2 result from the ~ cohort + diet model above.

library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(cowplot)

run_gobp_enrichment = function(gsea_rank_by) {
  gsea_rank_column <- case_when(
    gsea_rank_by == "stat" ~ "stat",
    gsea_rank_by == "logFC" ~ "log2FoldChange",
    TRUE ~ NA_character_
  )
  stopifnot(!is.na(gsea_rank_column))
  
  gobp_de = res2 %>%
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(pvalue),
      !is.na(.data[[gsea_rank_column]])
    )
  
  gene_symbol_to_entrez = AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(gobp_de$gene.id),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  ) %>%
    as_tibble() %>%
    dplyr::filter(!is.na(ENTREZID))
  
  gobp_ranked_genes = gobp_de %>%
    mutate(rank_value = .data[[gsea_rank_column]]) %>%
    inner_join(gene_symbol_to_entrez, by = c("gene.id" = "SYMBOL")) %>%
    arrange(desc(abs(rank_value))) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  gobp_gene_ranks = gobp_ranked_genes$rank_value
  names(gobp_gene_ranks) = gobp_ranked_genes$ENTREZID
  gobp_gene_ranks = sort(gobp_gene_ranks, decreasing = TRUE)
  
  gsea_gobp = gseGO(
    geneList = gobp_gene_ranks,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.10,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE
  )
  
  gsea_gobp_tbl = as.data.frame(gsea_gobp) %>%
    as_tibble() %>%
    arrange(p.adjust)
  gsea_gobp_table_file = here(
    "data/processed/pathway_enrichment/go_gsea",
    paste0("go bp gsea old ID vs CD ", gsea_rank_by, " ranking.csv")
  )
  write.csv(gsea_gobp_tbl, gsea_gobp_table_file, row.names = FALSE)
  
  if (nrow(gsea_gobp_tbl) > 0) {
    gsea_gobp_plot_data = gsea_gobp_tbl %>%
      mutate(
        direction = ifelse(NES > 0, "Activated in ID", "Suppressed in ID"),
        direction = factor(direction, levels = c("Activated in ID", "Suppressed in ID")),
        pathway = stringr::str_wrap(Description, width = 60)
      ) %>%
      group_by(direction) %>%
      slice_min(p.adjust, n = 20, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(direction, NES) %>%
      mutate(
        pathway_for_plot = paste(pathway, direction, sep = "___"),
        pathway_for_plot = factor(pathway_for_plot, levels = unique(pathway_for_plot)),
        neg_log10_padj = -log10(pmax(p.adjust, .Machine$double.xmin))
      )
    
    gsea_nes_limit = max(abs(gsea_gobp_plot_data$NES), na.rm = TRUE)
    gsea_nes_limit = ceiling(gsea_nes_limit * 2) / 2
    
    gsea_gobp_plot = ggplot(gsea_gobp_plot_data, aes(NES, pathway_for_plot)) +
      geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
      geom_point(aes(size = setSize, color = neg_log10_padj), alpha = 0.9) +
      facet_grid(
        direction ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
      scale_color_viridis_c(name = expression(-log[10]~adjusted~italic(P))) +
      coord_cartesian(xlim = c(-gsea_nes_limit, gsea_nes_limit)) +
      labs(
        title = paste0("GO BP GSEA ranked by ", gsea_rank_by),
        x = "Normalized enrichment score",
        y = NULL,
        size = "Gene set size"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8, lineheight = 0.9),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 28),
        legend.position = "bottom",
        legend.box = "horizontal"
      ) +
      guides(
        color = guide_colorbar(
          title.position = "top",
          barwidth = grid::unit(35, "mm"),
          barheight = grid::unit(3, "mm")
        ),
        size = guide_legend(title.position = "top")
      )
  } else {
    gsea_gobp_plot_data = tibble()
    gsea_gobp_plot = ggdraw() +
      draw_label("No GO BP GSEA terms passed the selected cutoff", size = 12)
  }
  
  # ORA uses the significant old ID vs CD genes from dex2, split by direction.
  gobp_universe = unique(gobp_ranked_genes$ENTREZID)
  
  ora_up_genes = dex2 %>%
    dplyr::filter(log2FoldChange > 0) %>%
    inner_join(gene_symbol_to_entrez, by = c("gene.id" = "SYMBOL")) %>%
    pull(ENTREZID) %>%
    unique()
  
  ora_down_genes = dex2 %>%
    dplyr::filter(log2FoldChange < 0) %>%
    inner_join(gene_symbol_to_entrez, by = c("gene.id" = "SYMBOL")) %>%
    pull(ENTREZID) %>%
    unique()
  
  ora_up_gobp = enrichGO(
    gene = ora_up_genes,
    universe = gobp_universe,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  ora_down_gobp = enrichGO(
    gene = ora_down_genes,
    universe = gobp_universe,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  ora_up_gobp_tbl = as.data.frame(ora_up_gobp) %>% as_tibble()
  ora_down_gobp_tbl = as.data.frame(ora_down_gobp) %>% as_tibble()
  
  if (nrow(ora_up_gobp_tbl) > 0) {
    ora_up_gobp_plot = dotplot(ora_up_gobp, showCategory = 10) +
      labs(title = "GO BP ORA: upregulated in ID", x = "Gene ratio", y = NULL) +
      theme_bw()
  } else {
    ora_up_gobp_plot = ggdraw() +
      draw_label("No GO BP ORA terms for upregulated genes", size = 12)
  }
  
  if (nrow(ora_down_gobp_tbl) > 0) {
    ora_down_gobp_plot = dotplot(ora_down_gobp, showCategory = 10) +
      labs(title = "GO BP ORA: downregulated in ID", x = "Gene ratio", y = NULL) +
      theme_bw()
  } else {
    ora_down_gobp_plot = ggdraw() +
      draw_label("No GO BP ORA terms for downregulated genes", size = 12)
  }
  
  ora_gobp_plot = plot_grid(ora_up_gobp_plot, ora_down_gobp_plot, nrow = 1)
  gobp_enrichment_plot = plot_grid(
    gsea_gobp_plot,
    ora_gobp_plot,
    ncol = 1,
    rel_heights = c(1.35, 1)
  )
  
  output_file = here(
    "plots/diet DEx/gsea",
    paste0("go bp gsea ora old ID vs CD ", gsea_rank_by, " ranking.pdf")
  )
  ggsave(output_file, gobp_enrichment_plot, width = 13, height = 12)
  
  list(
    gsea = gsea_gobp,
    gsea_tbl = gsea_gobp_tbl,
    gsea_plot_data = gsea_gobp_plot_data,
    gsea_plot = gsea_gobp_plot,
    ora_up = ora_up_gobp,
    ora_up_tbl = ora_up_gobp_tbl,
    ora_down = ora_down_gobp,
    ora_down_tbl = ora_down_gobp_tbl,
    ora_plot = ora_gobp_plot,
    combined_plot = gobp_enrichment_plot,
    gsea_table_file = gsea_gobp_table_file,
    output_file = output_file
  )
}

gsea_rank_modes = c("stat", "logFC")
gobp_enrichment_results = setNames(
  lapply(gsea_rank_modes, run_gobp_enrichment),
  gsea_rank_modes
)

gobp_stat = gobp_enrichment_results$stat
gobp_logFC = gobp_enrichment_results$logFC

# Keep stat-ranked results in the original object names for interactive follow-up.
gsea_gobp = gobp_stat$gsea
gsea_gobp_tbl = gobp_stat$gsea_tbl
gsea_gobp_plot = gobp_stat$gsea_plot

ora_up_gobp = gobp_stat$ora_up
ora_up_gobp_tbl = gobp_stat$ora_up_tbl
ora_down_gobp = gobp_stat$ora_down
ora_down_gobp_tbl = gobp_stat$ora_down_tbl
ora_gobp_plot = gobp_stat$ora_plot

gobp_enrichment_plot = gobp_stat$combined_plot
gobp_enrichment_plot


gsea_gobp_plot
ggsave(here("plots/diet DEx/gsea/go bp gsea.pdf"),h=9,w=7)




# ========== 2.1 - GO CC GSEA ==========
# Uses the same old ID vs CD DESeq2 result and ranking options as the GO BP GSEA.

run_gocc_gsea = function(gsea_rank_by) {
  gsea_rank_column <- case_when(
    gsea_rank_by == "stat" ~ "stat",
    gsea_rank_by == "logFC" ~ "log2FoldChange",
    TRUE ~ NA_character_
  )
  stopifnot(!is.na(gsea_rank_column))
  
  gocc_de = res2 %>%
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(pvalue),
      !is.na(.data[[gsea_rank_column]])
    )
  
  gene_symbol_to_entrez = AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(gocc_de$gene.id),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  ) %>%
    as_tibble() %>%
    dplyr::filter(!is.na(ENTREZID))
  
  gocc_ranked_genes = gocc_de %>%
    mutate(rank_value = .data[[gsea_rank_column]]) %>%
    inner_join(gene_symbol_to_entrez, by = c("gene.id" = "SYMBOL")) %>%
    arrange(desc(abs(rank_value))) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  gocc_gene_ranks = gocc_ranked_genes$rank_value
  names(gocc_gene_ranks) = gocc_ranked_genes$ENTREZID
  gocc_gene_ranks = sort(gocc_gene_ranks, decreasing = TRUE)
  
  gsea_gocc = gseGO(
    geneList = gocc_gene_ranks,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.10,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE
  )
  
  gsea_gocc_tbl = as.data.frame(gsea_gocc) %>%
    as_tibble() %>%
    arrange(p.adjust)
  gsea_gocc_table_file = here(
    "data/processed/pathway_enrichment/go_gsea",
    paste0("go cc gsea old ID vs CD ", gsea_rank_by, " ranking.csv")
  )
  write.csv(gsea_gocc_tbl, gsea_gocc_table_file, row.names = FALSE)
  
  if (nrow(gsea_gocc_tbl) > 0) {
    gsea_gocc_plot_data = gsea_gocc_tbl %>%
      mutate(
        direction = ifelse(NES > 0, "Activated in ID", "Suppressed in ID"),
        direction = factor(direction, levels = c("Activated in ID", "Suppressed in ID")),
        pathway = stringr::str_wrap(Description, width = 60)
      ) %>%
      group_by(direction) %>%
      slice_min(p.adjust, n = 20, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(direction, NES) %>%
      mutate(
        pathway_for_plot = paste(pathway, direction, sep = "___"),
        pathway_for_plot = factor(pathway_for_plot, levels = unique(pathway_for_plot)),
        neg_log10_padj = -log10(pmax(p.adjust, .Machine$double.xmin))
      )
    
    gocc_nes_limit = max(abs(gsea_gocc_plot_data$NES), na.rm = TRUE)
    gocc_nes_limit = ceiling(gocc_nes_limit * 2) / 2
    
    gsea_gocc_plot = ggplot(gsea_gocc_plot_data, aes(NES, pathway_for_plot)) +
      geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
      geom_point(aes(size = setSize, color = neg_log10_padj), alpha = 0.9) +
      facet_grid(
        direction ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
      scale_color_viridis_c(name = expression(-log[10]~adjusted~italic(P))) +
      coord_cartesian(xlim = c(-gocc_nes_limit, gocc_nes_limit)) +
      labs(
        title = paste0("GO CC GSEA ranked by ", gsea_rank_by),
        x = "Normalized enrichment score",
        y = NULL,
        size = "Gene set size"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8, lineheight = 0.9),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 28),
        legend.position = "bottom",
        legend.box = "horizontal"
      ) +
      guides(
        color = guide_colorbar(
          title.position = "top",
          barwidth = grid::unit(35, "mm"),
          barheight = grid::unit(3, "mm")
        ),
        size = guide_legend(title.position = "top")
      )
  } else {
    gsea_gocc_plot_data = tibble()
    gsea_gocc_plot = ggdraw() +
      draw_label("No GO CC GSEA terms passed the selected cutoff", size = 12)
  }
  
  output_file = here(
    "plots/diet DEx/gsea",
    paste0("go cc gsea old ID vs CD ", gsea_rank_by, " ranking.pdf")
  )
  ggsave(output_file, gsea_gocc_plot, width = 10, height = 8)
  
  list(
    gsea = gsea_gocc,
    gsea_tbl = gsea_gocc_tbl,
    gsea_plot_data = gsea_gocc_plot_data,
    gsea_plot = gsea_gocc_plot,
    gsea_table_file = gsea_gocc_table_file,
    output_file = output_file
  )
}

gocc_enrichment_results = setNames(
  lapply(gsea_rank_modes, run_gocc_gsea),
  gsea_rank_modes
)

gocc_stat = gocc_enrichment_results$stat
gocc_logFC = gocc_enrichment_results$logFC
gsea_gocc = gocc_stat$gsea
gsea_gocc_tbl = gocc_stat$gsea_tbl
gsea_gocc_plot = gocc_stat$gsea_plot
gsea_gocc_plot
ggsave(here("plots/diet DEx/gsea/go cc gsea.pdf"),h=9,w=7)

# ========== 2.2 - MSigDB Hallmark GSEA ==========
# Uses the old ID vs CD DESeq2 result ranked by the Wald statistic.

hallmark_term2gene = msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "MH"
) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

run_hallmark_gsea = function(gsea_rank_by) {
  gsea_rank_column <- case_when(
    gsea_rank_by == "stat" ~ "stat",
    gsea_rank_by == "logFC" ~ "log2FoldChange",
    TRUE ~ NA_character_
  )
  stopifnot(!is.na(gsea_rank_column))
  
  hallmark_ranked_genes = res2 %>%
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(pvalue),
      !is.na(.data[[gsea_rank_column]])
    ) %>%
    mutate(rank_value = .data[[gsea_rank_column]]) %>%
    semi_join(hallmark_term2gene, by = c("gene.id" = "gene_symbol")) %>%
    arrange(desc(abs(rank_value))) %>%
    distinct(gene.id, .keep_all = TRUE)
  
  hallmark_gene_ranks = hallmark_ranked_genes$rank_value
  names(hallmark_gene_ranks) = hallmark_ranked_genes$gene.id
  hallmark_gene_ranks = sort(hallmark_gene_ranks, decreasing = TRUE)
  
  hallmark_gsea = GSEA(
    geneList = hallmark_gene_ranks,
    TERM2GENE = hallmark_term2gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE
  )
  
  hallmark_gsea_tbl = as.data.frame(hallmark_gsea) %>%
    as_tibble() %>%
    arrange(p.adjust)
  hallmark_gsea_table_file = here(
    "data/processed/pathway_enrichment/msigdb_hallmark",
    paste0("msigdb hallmark gsea old ID vs CD ", gsea_rank_by, " ranking.csv")
  )
  write.csv(hallmark_gsea_tbl, hallmark_gsea_table_file, row.names = FALSE)
  
  if (nrow(hallmark_gsea_tbl) > 0) {
    hallmark_gsea_plot_data = hallmark_gsea_tbl %>%
      mutate(
        direction = ifelse(NES > 0, "Activated in ID", "Suppressed in ID"),
        direction = factor(direction, levels = c("Activated in ID", "Suppressed in ID")),
        pathway = Description %>%
          str_remove("^HALLMARK_") %>%
          str_replace_all("_", " ") %>%
          str_to_sentence() %>%
          stringr::str_wrap(width = 60)
      ) %>%
      group_by(direction) %>%
      slice_min(p.adjust, n = 20, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(direction, NES) %>%
      mutate(
        pathway_for_plot = paste(pathway, direction, sep = "___"),
        pathway_for_plot = factor(pathway_for_plot, levels = unique(pathway_for_plot)),
        neg_log10_padj = -log10(pmax(p.adjust, .Machine$double.xmin))
      )
    
    hallmark_nes_limit = max(abs(hallmark_gsea_plot_data$NES), na.rm = TRUE)
    hallmark_nes_limit = ceiling(hallmark_nes_limit * 2) / 2
    
    hallmark_gsea_plot = ggplot(hallmark_gsea_plot_data, aes(NES, pathway_for_plot)) +
      geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
      geom_point(aes(size = setSize, color = neg_log10_padj), alpha = 0.9) +
      facet_grid(
        direction ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
      scale_color_viridis_c(name = expression(-log[10]~adjusted~italic(P))) +
      coord_cartesian(xlim = c(-hallmark_nes_limit, hallmark_nes_limit)) +
      labs(
        title = paste0("MSigDB Hallmark GSEA ranked by ", gsea_rank_by),
        x = "Normalized enrichment score",
        y = NULL,
        size = "Gene set size"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 9, lineheight = 0.9),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 28),
        legend.position = "bottom",
        legend.box = "horizontal"
      ) +
      guides(
        color = guide_colorbar(
          title.position = "top",
          barwidth = grid::unit(35, "mm"),
          barheight = grid::unit(3, "mm")
        ),
        size = guide_legend(title.position = "top")
      )
  } else {
    hallmark_gsea_plot_data = tibble()
    hallmark_gsea_plot = ggdraw() +
      draw_label("No Hallmark GSEA terms available", size = 12)
  }
  
  output_file = here(
    "plots/diet DEx/gsea/",
    paste0("msigdb gsea.pdf")
  )
  ggsave(output_file, hallmark_gsea_plot, width = 10, height = 8)
  
  list(
    gsea = hallmark_gsea,
    gsea_tbl = hallmark_gsea_tbl,
    gsea_plot_data = hallmark_gsea_plot_data,
    gsea_plot = hallmark_gsea_plot,
    gsea_table_file = hallmark_gsea_table_file,
    output_file = output_file
  )
}

hallmark_rank_mode = "stat"

hallmark_enrichment_results = setNames(
  lapply(hallmark_rank_mode, run_hallmark_gsea),
  hallmark_rank_mode
)

hallmark_stat = hallmark_enrichment_results$stat
hallmark_gsea = hallmark_stat$gsea
hallmark_gsea_tbl = hallmark_stat$gsea_tbl
hallmark_gsea_plot = hallmark_stat$gsea_plot
hallmark_gsea_plot


# ========== 3.0 - Summary figure: volcano plus stat-ranked GSEA (not final) ==========
# Volcano highlights significant genes in the top activated/suppressed GO BP
# GSEA terms from the stat-ranked analysis.

top_gobp_terms_per_direction = 4
top_gobp_pathways = gsea_gobp_tbl %>%
  dplyr::filter(!is.na(NES)) %>%
  mutate(direction = ifelse(NES > 0, "Activated", "Suppressed")) %>%
  group_by(direction) %>%
  slice_min(p.adjust, n = top_gobp_terms_per_direction, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(factor(direction, levels = c("Activated", "Suppressed")), p.adjust) %>%
  mutate(
    pathway_rank = row_number(),
    pathway_label = paste0(
      direction,
      ": ",
      stringr::str_wrap(Description, width = 34)
    )
  )

# Use full GO term membership for volcano coloring, not only GSEA leading-edge genes.
top_gobp_symbols = AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(top_gobp_pathways$ID),
  keytype = "GOALL",
  columns = c("GOALL", "SYMBOL", "ONTOLOGYALL")
) %>%
  as_tibble() %>%
  dplyr::filter(
    GOALL %in% top_gobp_pathways$ID,
    ONTOLOGYALL == "BP",
    !is.na(SYMBOL)
  ) %>%
  distinct(GOALL, SYMBOL)

top_gobp_gene_categories = top_gobp_symbols %>%
  inner_join(
    top_gobp_pathways %>% dplyr::select(GOALL = ID, pathway_rank, pathway_label),
    by = "GOALL"
  ) %>%
  arrange(pathway_rank) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(gene.id = SYMBOL, pathway_label)

top_up_gene_labels = dex2 %>%
  dplyr::filter(log2FoldChange > 0) %>%
  arrange(pvalue, desc(log2FoldChange)) %>%
  slice_head(n = 15) %>%
  pull(gene.id)

top_down_gene_labels = dex2 %>%
  dplyr::filter(log2FoldChange < 0) %>%
  arrange(pvalue, log2FoldChange) %>%
  slice_head(n = 15) %>%
  pull(gene.id)

summary_volcano_data = lfc %>%
  mutate(is_significant = gene.id %in% dex2$gene.id) %>%
  left_join(top_gobp_gene_categories, by = "gene.id") %>%
  mutate(
    pathway_label = ifelse(!is.na(pathway_label), pathway_label, "Other"),
    pathway_label = factor(pathway_label, levels = c(unique(top_gobp_pathways$pathway_label), "Other")),
    gene_label = ifelse(is_significant & pathway_label != "Other", gene.id, NA_character_)
  )

summary_volcano_labels = summary_volcano_data %>%
  dplyr::filter(!is.na(gene_label), is_significant, pathway_label != "Other")

summary_pathway_palette = c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#7F3C8D", "#11A579"
)

summary_volcano_colors = c(
  setNames(
    summary_pathway_palette[seq_len(nrow(top_gobp_pathways))],
    top_gobp_pathways$pathway_label
  ),
  "Other" = "#D0D0D0"
)

summary_x_limit = max(2, quantile(abs(summary_volcano_data$log2FoldChange), 0.995, na.rm = TRUE))
summary_x_limit = min(3, ceiling(summary_x_limit * 2) / 2)
summary_y_limit = max(2, quantile(summary_volcano_data$neg_log10_pvalue, 0.995, na.rm = TRUE))
summary_y_limit = min(4, ceiling(summary_y_limit * 2) / 2)

summary_volcano = ggplot(summary_volcano_data, aes(log2FoldChange, neg_log10_pvalue)) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", color = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = c(-1, 1), linetype = "22", color = "grey55", linewidth = 0.35) +
  geom_point(
    data = summary_volcano_data %>% dplyr::filter(pathway_label == "Other"),
    color = "#D0D0D0",
    size = 1.1,
    alpha = 0.45,
    stroke = 0
  ) +
  geom_point(
    data = summary_volcano_data %>% dplyr::filter(pathway_label != "Other"),
    aes(color = pathway_label),
    size = 2.2,
    alpha = 0.95,
    stroke = 0
  ) +
  ggrepel::geom_text_repel(
    data = summary_volcano_labels,
    aes(label = gene_label),
    color = "grey15",
    size = 3.5,
    box.padding = 0.45,
    point.padding = 0.35,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.color = "grey70",
    segment.size = 0.3,
    segment.alpha = 0.7,
    show.legend = FALSE
  ) +
  scale_color_manual(values = summary_volcano_colors, name = "Top GO BP pathways", drop = FALSE) +
  coord_cartesian(xlim = c(-summary_x_limit, summary_x_limit), ylim = c(0, summary_y_limit)) +
  labs(
    title = "Old muscle: inulin diet vs control",
    subtitle = "Highlighted points are genes in top stat-ranked GO BP pathways",
    x = expression(log[2]~fold~change~"(ID/CD)"),
    y = expression(-log[10]~italic(P))
  ) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.55),
    axis.text = element_text(color = "black", size = 42),
    axis.title = element_text(color = "black", size = 52),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.position = "top",
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)))

volcano_gsea_summary_figure = plot_grid(
  summary_volcano,
  hallmark_gsea_plot,
  gsea_gobp_plot,
  gsea_gocc_plot,
  ncol = 1,
  rel_heights = c(1.1, 1, 1.25, 1.25)
)

# ========== 3.1 - Summary figure with MLL interesting DEG labels (not final) ==========
# Same top GO BP pathway coloring and stat-ranked GSEA panels as above, but
# label genes with any entry in the MLL_interesting DEGs flag column.

mll_flagged_gene_labels = readr::read_csv(
  here("data/processed/overlaps/MLL_interesting DEGs.csv"),
  show_col_types = FALSE
) %>%
  mutate(flag = stringr::str_trim(as.character(flag))) %>%
  dplyr::filter(!is.na(flag), flag != "") %>%
  distinct(gene.id, .keep_all = TRUE) %>%
  dplyr::select(gene.id, flag)

top_gobp_terms_per_direction = 4
top_gobp_pathways = gsea_gobp_tbl %>%
  dplyr::filter(!is.na(NES)) %>%
  mutate(direction = ifelse(NES > 0, "Activated", "Suppressed")) %>%
  group_by(direction) %>%
  slice_min(p.adjust, n = top_gobp_terms_per_direction, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(factor(direction, levels = c("Activated", "Suppressed")), p.adjust) %>%
  mutate(
    pathway_rank = row_number(),
    pathway_label = paste0(
      direction,
      ": ",
      stringr::str_wrap(Description, width = 34)
    )
  )

# Use full GO term membership for volcano coloring, not only GSEA leading-edge genes.
top_gobp_symbols = AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(top_gobp_pathways$ID),
  keytype = "GOALL",
  columns = c("GOALL", "SYMBOL", "ONTOLOGYALL")
) %>%
  as_tibble() %>%
  dplyr::filter(
    GOALL %in% top_gobp_pathways$ID,
    ONTOLOGYALL == "BP",
    !is.na(SYMBOL)
  ) %>%
  distinct(GOALL, SYMBOL)

top_gobp_gene_categories = top_gobp_symbols %>%
  inner_join(
    top_gobp_pathways %>% dplyr::select(GOALL = ID, pathway_rank, pathway_label),
    by = "GOALL"
  ) %>%
  arrange(pathway_rank) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(gene.id = SYMBOL, pathway_label)

mll_summary_volcano_data = lfc %>%
  mutate(is_significant = gene.id %in% dex2$gene.id) %>%
  left_join(top_gobp_gene_categories, by = "gene.id") %>%
  left_join(mll_flagged_gene_labels, by = "gene.id") %>%
  mutate(
    pathway_label = ifelse(!is.na(pathway_label), pathway_label, "Other"),
    pathway_label = factor(pathway_label, levels = c(unique(top_gobp_pathways$pathway_label), "Other")),
    gene_label = ifelse(is_significant & pathway_label != "Other", gene.id, NA_character_)
  )

mll_summary_volcano_labels = mll_summary_volcano_data %>%
  dplyr::filter(!is.na(gene_label), is_significant, pathway_label != "Other", is.na(flag))

mll_flagged_volcano_labels = mll_summary_volcano_data %>%
  dplyr::filter(!is.na(flag))

mll_summary_pathway_palette = c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#7F3C8D", "#11A579"
)

mll_summary_volcano_colors = c(
  setNames(
    mll_summary_pathway_palette[seq_len(nrow(top_gobp_pathways))],
    top_gobp_pathways$pathway_label
  ),
  "Other" = "#D0D0D0"
)

mll_summary_x_limit = max(2, quantile(abs(mll_summary_volcano_data$log2FoldChange), 0.995, na.rm = TRUE))
mll_summary_x_limit = min(3, ceiling(mll_summary_x_limit * 2) / 2)
mll_summary_y_limit = max(2, quantile(mll_summary_volcano_data$neg_log10_pvalue, 0.995, na.rm = TRUE))
mll_summary_y_limit = min(4, ceiling(mll_summary_y_limit * 2) / 2)

mll_summary_volcano = ggplot(mll_summary_volcano_data, aes(log2FoldChange, neg_log10_pvalue)) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", color = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = c(-1, 1), linetype = "22", color = "grey55", linewidth = 0.35) +
  geom_point(
    data = mll_summary_volcano_data %>% dplyr::filter(pathway_label == "Other"),
    color = "#D0D0D0",
    size = 1.1,
    alpha = 0.45,
    stroke = 0
  ) +
  geom_point(
    data = mll_summary_volcano_data %>% dplyr::filter(pathway_label != "Other"),
    aes(color = pathway_label),
    size = 2.2,
    alpha = 0.95,
    stroke = 0
  ) +
  ggrepel::geom_text_repel(
    data = mll_summary_volcano_labels,
    aes(label = gene_label),
    color = "grey15",
    size = 3.5,
    box.padding = 0.45,
    point.padding = 0.35,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.color = "grey70",
    segment.size = 0.3,
    segment.alpha = 0.7,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = mll_flagged_volcano_labels,
    aes(label = gene.id),
    color = "#B2182B",
    fontface = "bold",
    size = 3.8,
    box.padding = 0.5,
    point.padding = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.color = "#B2182B",
    segment.size = 0.35,
    segment.alpha = 0.75,
    show.legend = FALSE
  ) +
  scale_color_manual(values = mll_summary_volcano_colors, name = "Top GO BP pathways", drop = FALSE) +
  coord_cartesian(xlim = c(-mll_summary_x_limit, mll_summary_x_limit), ylim = c(0, mll_summary_y_limit)) +
  labs(
    title = "Old muscle: inulin diet vs control",
    subtitle = "Bold red labels are genes flagged in MLL_interesting DEGs",
    x = expression(log[2]~fold~change~"(ID/CD)"),
    y = expression(-log[10]~italic(P))
  ) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.position = "top",
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)))

mll_volcano_gsea_summary_figure = plot_grid(
  mll_summary_volcano,
  hallmark_gsea_plot,
  gsea_gobp_plot,
  gsea_gocc_plot,
  ncol = 1,
  rel_heights = c(1.1, 1, 1.25, 1.25)
)

# ========== 4.0 - MLL Figure ==========
# --

mll_superpathway_label_genes = c(
  # ECM remodeling
  "Ddah1", "Vwf",
  # Cell-matrix adhesion
  "Itga2", "Syne2",
  # Cytoskeletal dynamics
  "Kcnc1", "Phactr1",
  # Immune signaling
  "Ccl8", "Spns2",
  # Mitochondrial respiration
  "mt-Rnr2", "Efhd1",
  # Mitochondrial translation
  "Fastkd3",
  # Proteostasis
  "Hspa5", "Socs2",
  # ER-Golgi trafficking
  "Copz1", "Rab3a"
)

mll_selected_label_genes = c(
  "Glcci1", "Clec11a", "Stab1", "Meg3", "Lep", "Hid1",
  "Socs2", "Gpx8", "Cish", "Folh1", "Ccr5", "Serpina3n",
  "Serpina1c", "Ttr", "Ahsg", "Wfdc", "Npr3", "Galr2", "Acsm3", "Chil3"
)

# -- review GOBP pathways for representative/non-redundant terms
gsea_gopb_sig = gsea_gobp_tbl %>%
  dplyr::filter(!is.na(NES), qvalue <= 0.05) %>%
  mutate(direction = ifelse(NES > 0, "Activated", "Suppressed")) %>%
  group_by(direction) %>%
  ungroup() %>%
  arrange(factor(direction, levels = c("Activated", "Suppressed")), p.adjust) %>%
  mutate(
    pathway_rank = row_number(),
    pathway_label = paste0(
      direction,
      ": ",
      stringr::str_wrap(Description, width = 34)
    )
  )

my_gobp_pwys_up = c("GO:0030198","GO:0043062","GO:0050900","GO:0060326","GO:0019882",
                 "GO:0001819","GO:0051251","GO:0043410","GO:0061572","GO:0001525")

my_gobp_pwys_down = c("GO:0006119","GO:0042775","GO:0032543","GO:0007029","GO:0000423",
                      "GO:0031929","GO:0008203","GO:0043161","GO:0048193","GO:0006413")

gsea_gobp_mypathways = gsea_gopb_sig %>% dplyr::filter(ID %in% c(my_gobp_pwys_down, my_gobp_pwys_up))


# -- review GOCC pathways for representative/non-redundant terms
gsea_gocc_sig = gsea_gocc_tbl %>%
  dplyr::filter(!is.na(NES), qvalue <= 0.05) %>%
  mutate(direction = ifelse(NES > 0, "Activated", "Suppressed")) %>%
  group_by(direction) %>%
  ungroup() %>%
  arrange(factor(direction, levels = c("Activated", "Suppressed")), p.adjust) %>%
  mutate(
    pathway_rank = row_number(),
    pathway_label = paste0(
      direction,
      ": ",
      stringr::str_wrap(Description, width = 34)
    )
  )

my_gocc_pwys_up = c("GO:0031012","GO:0005581","GO:0005604","GO:0005925","GO:0008305",
                    "GO:0015629","GO:0031252","GO:0042611","GO:0043235","GO:0045335")

my_gocc_pwys_down = c("GO:0098798","GO:0005743","GO:0005741","GO:0005761","GO:0098803",
                      "GO:0030964","GO:0000502","GO:0005776","GO:0030133","GO:0030134")

gsea_gocc_mypathways = gsea_gocc_sig %>% dplyr::filter(ID %in% c(my_gocc_pwys_down, my_gocc_pwys_up))


# -- TAKE GENES FROM ABOVE TERMS AND MERGE INTO SUPERCATEGORIES --
# Activated
# 
# ECM remodeling
# Cell–matrix adhesion
# Cytoskeletal dynamics
# Immune signaling
# 
# Suppressed
# 
# Mitochondrial respiration
# Mitochondrial translation
# Proteostasis
# ER–Golgi trafficking

mll_superpathway_ids = list(
  "ECM remodeling" = c("GO:0030198", "GO:0043062", "GO:0001525",
                       "GO:0031012", "GO:0005581", "GO:0005604"),
  "Cell-matrix adhesion" = c("GO:0005925", "GO:0008305"),
  "Cytoskeletal dynamics" = c("GO:0061572", "GO:0015629", "GO:0031252"),
  "Immune signaling" = c("GO:0050900", "GO:0060326", "GO:0019882",
                         "GO:0001819", "GO:0051251", "GO:0043410",
                         "GO:0042611", "GO:0043235"),
  "Mitochondrial respiration" = c("GO:0006119", "GO:0042775",
                                  "GO:0098798", "GO:0005743", "GO:0005741",
                                  "GO:0098803", "GO:0030964"),
  "Mitochondrial translation" = c("GO:0032543", "GO:0005761"),
  "Proteostasis" = c("GO:0043161", "GO:0006413", "GO:0000502"),
  "ER-Golgi trafficking" = c("GO:0007029", "GO:0048193",
                             "GO:0030133", "GO:0030134")
)

mll_superpathway_terms = enframe(mll_superpathway_ids, name = "superpathway", value = "ID") %>%
  unnest(ID)

mll_nes_term_group_ids = c(
  mll_superpathway_ids,
  list(
    "Autophagy / TOR signaling" = c("GO:0000423", "GO:0031929", "GO:0005776"),
    "Cholesterol metabolism" = "GO:0008203",
    "Phagocytic activity" = "GO:0045335"
  )
)

mll_nes_term_group_levels = names(mll_nes_term_group_ids)

mll_nes_term_groups = enframe(mll_nes_term_group_ids, name = "superpathway", value = "ID") %>%
  unnest(ID)

mll_all_selected_go_terms = bind_rows(
  gsea_gobp_mypathways %>% mutate(ontology = "GO BP"),
  gsea_gocc_mypathways %>% mutate(ontology = "GO CC")
)

mll_selected_go_terms_missing_superpathway = mll_all_selected_go_terms %>%
  anti_join(mll_nes_term_groups, by = "ID") %>%
  dplyr::select(ontology, direction, ID, Description, NES, p.adjust, qvalue)

mll_selected_go_terms = mll_all_selected_go_terms %>%
  inner_join(mll_nes_term_groups, by = "ID") %>%
  arrange(factor(direction, levels = c("Activated", "Suppressed")), superpathway, ontology, p.adjust)

mll_selected_go_terms
mll_selected_go_terms_missing_superpathway

# -- Build refined MLL figure with superpathway volcano coloring
mll_superpathway_levels = names(mll_superpathway_ids)

mll_visual_category_map = tribble(
  ~superpathway, ~visual_category,
  "ECM remodeling", "ECM / adhesion / cytoskeleton",
  "Cell-matrix adhesion", "ECM / adhesion / cytoskeleton",
  "Cytoskeletal dynamics", "ECM / adhesion / cytoskeleton",
  "Immune signaling", "immune signaling",
  "Mitochondrial respiration", "mitochondrial / translation",
  "Mitochondrial translation", "mitochondrial / translation",
  "Proteostasis", "proteostasis / trafficking",
  "ER-Golgi trafficking", "proteostasis / trafficking"
)

mll_visual_category_levels = c(
  "ECM / adhesion / cytoskeleton",
  "immune signaling",
  "mitochondrial / translation",
  "proteostasis / trafficking"
)

mll_visual_category_colors = c(
  "ECM / adhesion / cytoskeleton" = "#0072B2",
  "immune signaling" = "#CC79A7",
  "mitochondrial / translation" = "#E69F00",
  "proteostasis / trafficking" = "#7F3C8D",
  "Other" = "#D0D0D0"
)

# Use the GSEA core enrichment genes for volcano coloring. These are the genes
# that contribute to the selected significant NES results, not all GO members.
mll_refined_core_gene_terms = mll_selected_go_terms %>%
  semi_join(mll_superpathway_terms, by = "ID") %>%
  dplyr::select(
    ID,
    ontology,
    direction,
    superpathway,
    term_padj = p.adjust,
    core_enrichment
  ) %>%
  dplyr::filter(!is.na(core_enrichment), core_enrichment != "") %>%
  separate_rows(core_enrichment, sep = "/") %>%
  transmute(
    ENTREZID = as.character(core_enrichment),
    ID,
    ontology,
    direction,
    superpathway,
    term_padj
  )

mll_refined_core_gene_symbols = AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(mll_refined_core_gene_terms$ENTREZID),
  keytype = "ENTREZID",
  columns = c("ENTREZID", "SYMBOL")
) %>%
  as_tibble() %>%
  dplyr::filter(!is.na(SYMBOL))

mll_refined_go_symbols = mll_refined_core_gene_terms %>%
  inner_join(mll_refined_core_gene_symbols, by = "ENTREZID") %>%
  mutate(superpathway = factor(superpathway, levels = mll_superpathway_levels)) %>%
  arrange(factor(direction, levels = c("Activated", "Suppressed")), superpathway, term_padj) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(gene.id = SYMBOL, superpathway)

mll_refined_go_symbol_counts = mll_refined_go_symbols %>%
  dplyr::count(superpathway, name = "core_gene_count")

mll_refined_go_symbol_counts

mll_refined_volcano_data = lfc %>%
  mutate(is_significant = gene.id %in% dex2$gene.id) %>%
  left_join(mll_refined_go_symbols, by = "gene.id") %>%
  left_join(mll_visual_category_map, by = "superpathway") %>%
  mutate(
    superpathway = as.character(superpathway),
    superpathway = ifelse(!is.na(superpathway), superpathway, "Other"),
    superpathway = factor(superpathway, levels = c(mll_superpathway_levels, "Other")),
    visual_category = ifelse(!is.na(visual_category), visual_category, "Other"),
    visual_category = factor(visual_category, levels = c(mll_visual_category_levels, "Other")),
    label_group = case_when(
      gene.id %in% mll_selected_label_genes ~ "MLL selected",
      gene.id %in% mll_superpathway_label_genes ~ "Superpathway representative",
      TRUE ~ NA_character_
    ),
    gene_label = ifelse(!is.na(label_group), gene.id, NA_character_)
  )

mll_refined_volcano_labels = mll_refined_volcano_data %>%
  dplyr::filter(!is.na(label_group))

mll_refined_x_limit = max(
  1.5,
  max(abs(mll_refined_volcano_labels$log2FoldChange), na.rm = TRUE) + 0.35
)
mll_refined_x_limit = ceiling(mll_refined_x_limit * 2) / 2
mll_refined_y_limit = max(
  2,
  quantile(mll_refined_volcano_data$neg_log10_pvalue, 0.995, na.rm = TRUE),
  max(mll_refined_volcano_labels$neg_log10_pvalue, na.rm = TRUE)
)
mll_refined_y_limit = ceiling((mll_refined_y_limit + 0.35) * 2) / 2

mll_refined_volcano = ggplot(mll_refined_volcano_data, aes(log2FoldChange, neg_log10_pvalue)) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", color = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = c(-1, 1), linetype = "22", color = "grey55", linewidth = 0.35) +
  geom_point(
    data = mll_refined_volcano_data %>% dplyr::filter(superpathway == "Other"),
    color = "#D0D0D0",
    size = 1.1,
    alpha = 0.3,
    stroke = 0
  ) +
  geom_point(
    data = mll_refined_volcano_data %>% dplyr::filter(
      label_group == "MLL selected",
      superpathway == "Other"
    ),
    color = "black",
    size = 2.2,
    alpha = 0.85,
    stroke = 0,
    show.legend = FALSE
  ) +
  geom_point(
    data = mll_refined_volcano_data %>% dplyr::filter(superpathway != "Other"),
    aes(color = visual_category),
    size = 2.2,
    alpha = 0.7,
    stroke = 0
  ) +
  ggrepel::geom_text_repel(
    data = mll_refined_volcano_labels,
    aes(label = gene_label),
    color = "black",
    fontface = "bold",
    size = 4.8,
    box.padding = 0.45,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.color = "grey35",
    segment.size = 0.45,
    segment.alpha = 0.95,
    segment.linetype = "solid",
    force = 1.5,
    force_pull = 1,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = mll_visual_category_colors,
    breaks = mll_visual_category_levels,
    name = "GSEA core enrichment groups:",
    drop = FALSE
  ) +
  coord_cartesian(xlim = c(-mll_refined_x_limit, mll_refined_x_limit), ylim = c(0, mll_refined_y_limit)) +
  labs(
    x = expression(log[2]~fold~change~"(Inulin Diet/Control Diet)"),
    y = expression(-log[10]~italic(P))
  ) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 20),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = grid::unit(6, "mm")
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4.2, alpha = 0.85)))

ggsave(
  here("plots/diet DEx/volcano with refined GO superpathways MLL labels.pdf"),
  mll_refined_volcano,
  width = 13,
  height = 7
)

mll_refined_go_label_width = 76

mll_refined_go_label_overrides = c(
  "proteasome-mediated ubiquitin-dependent protein catabolic process" =
    "Proteasome/ubiquitin protein catabolism",
  "mitochondrial ATP synthesis coupled electron transport" =
    "Mitochon. ATP synthesis/electron transport"
)

mll_refined_gsea_terms_for_plot = c(
  # -- ECM remodeling
  "GO:0030198", 
  "GO:0043062", 
  "GO:0001525", 
  "GO:0031012", 
  "GO:0005581", 
  #"GO:0005604",
  # -- Cell-matrix adhesion
  #"GO:0005925",
  #"GO:0008305",
  # -- Cytoskeletal dynamics
  "GO:0061572", 
  "GO:0015629", 
  #"GO:0031252",
  # -- Immune signaling
  #"GO:0050900",
  #"GO:0060326",
  #"GO:0019882",
  "GO:0001819", 
  #"GO:0051251",
  #"GO:0043410",
  #"GO:0042611",
  #"GO:0043235",
  # Mitochondrial respiration
  "GO:0006119", 
  "GO:0042775", 
  "GO:0098798", 
  "GO:0005743", 
  "GO:0005741", 
  "GO:0098803", 
  "GO:0030964", 
  # Mitochondrial translation
  "GO:0032543", 
  "GO:0005761", 
  # Proteostasis
  #"GO:0043161",
  #"GO:0006413",
  #"GO:0000502",
  # ER-Golgi trafficking
  "GO:0007029", 
  #"GO:0048193",
  #"GO:0030133",
  #"GO:0030134",
  # Additional NES panel terms
  #"GO:0000423",
  "GO:0031929", 
  #"GO:0005776",
  "GO:0008203" 
  #"GO:0045335"
)

mll_refined_gsea_selected_terms = mll_selected_go_terms %>%
  dplyr::filter(ID %in% mll_refined_gsea_terms_for_plot)

mll_refined_gsea_selected_terms %>%
  dplyr::select(superpathway, ontology, direction, ID, Description, NES, p.adjust, qvalue)

mll_refined_nes_axis_limit = max(abs(mll_refined_gsea_selected_terms$NES), na.rm = TRUE)
mll_refined_nes_axis_limit = ceiling(mll_refined_nes_axis_limit * 2) / 2

mll_refined_nes_color_limits = range(
  -log10(pmax(mll_refined_gsea_selected_terms$p.adjust, .Machine$double.xmin)),
  na.rm = TRUE
)

mll_refined_nes_size_limits = range(mll_refined_gsea_selected_terms$setSize, na.rm = TRUE)

mll_refined_gsea_facet_levels = c(
  "Activated GOBP",
  "Activated GOCC",
  "Suppressed GOBP",
  "Suppressed GOCC"
)

mll_refined_gsea_facet_labels = c(
  "Activated GOBP" = "Activated\nGOBP",
  "Activated GOCC" = "Activated\nGOCC",
  "Suppressed GOBP" = "Suppressed\nGOBP",
  "Suppressed GOCC" = "Suppressed\nGOCC"
)

mll_refined_gsea_plot_data = mll_refined_gsea_selected_terms %>%
  mutate(
    direction = factor(direction, levels = c("Activated", "Suppressed")),
    ontology_label = recode(ontology, "GO BP" = "GOBP", "GO CC" = "GOCC"),
    facet_label = factor(
      paste(direction, ontology_label),
      levels = mll_refined_gsea_facet_levels
    ),
    superpathway = factor(superpathway, levels = mll_nes_term_group_levels),
    display_description = recode(
      Description,
      !!!mll_refined_go_label_overrides,
      .default = paste0(
        stringr::str_to_upper(stringr::str_sub(Description, 1, 1)),
        stringr::str_sub(Description, 2)
      )
    ),
    pathway = stringr::str_wrap(
      paste0(
        display_description,
        " (",
        ID,
        ")"
      ),
      width = mll_refined_go_label_width
    ),
    neg_log10_padj = -log10(pmax(p.adjust, .Machine$double.xmin))
  ) %>%
  arrange(facet_label, dplyr::desc(NES)) %>%
  mutate(
    pathway_for_plot = paste(pathway, facet_label, sep = "___"),
    pathway_for_plot = factor(pathway_for_plot, levels = rev(unique(pathway_for_plot)))
  )

mll_refined_gsea_plot = ggplot(mll_refined_gsea_plot_data, aes(NES, pathway_for_plot)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
  geom_point(aes(size = setSize, color = neg_log10_padj), alpha = 0.9) +
  facet_grid(
    facet_label ~ .,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(facet_label = mll_refined_gsea_facet_labels)
  ) +
  scale_y_discrete(
    labels = function(x) sub("___.*$", "", x),
    expand = expansion(add = c(0.55, 0.55))
  ) +
  scale_color_viridis_c(
    name = expression(-log[10]~adjusted~italic(P)),
    limits = mll_refined_nes_color_limits
  ) +
  scale_size_continuous(
    name = "Gene set size",
    limits = mll_refined_nes_size_limits,
    range = c(2.4, 6.2)
  ) +
  coord_cartesian(xlim = c(-mll_refined_nes_axis_limit, mll_refined_nes_axis_limit)) +
  labs(
    x = "Normalized enrichment score",
    y = NULL
  ) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13, color = "black", margin = margin(3, 0, 3, 0)),
    axis.text.y = element_text(size = 11.5, lineheight = 0.9, color = "black"),
    axis.text.x = element_text(size = 11.5, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(4, 10, 6, 10),
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      barwidth = grid::unit(35, "mm"),
      barheight = grid::unit(3, "mm"),
      order = 1
    ),
    size = guide_legend(title.position = "top", order = 2)
  )

mll_refined_volcano_gsea_summary_figure = plot_grid(
  mll_refined_volcano,
  mll_refined_gsea_plot,
  ncol = 1,
  rel_heights = c(1, 2)
)

mll_refined_volcano_gsea_summary_figure

ggsave(
  here("plots/diet DEx/volcano and gsea with refined GO superpathways MLL labels.pdf"),
  mll_refined_volcano_gsea_summary_figure,
  width = 13,
  height = 15
)

mll_refined_gsea_plot
ggsave(
  here("plots/diet DEx/gsea/go bp and cc refined pathways.pdf"),
  mll_refined_gsea_plot,
  width = 12,
  height = 8
)


# ========== z-old - Bar plots ==========
# --
# Start from your lfc object
df2 <- lfc %>%
  mutate(
    ID = 2^log2FoldChange,  # fold-change vs CD
    CD = 1                  # reference group
  ) %>%
  dplyr::select(gene.id, ID, CD, pvalue) %>%
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


volcano_lfc_barplot = plot_grid(
  vol_cat,
  gg_graber,
  plot_grid(gg_collagen, gg_top, gg_mac, nrow = 1),
  ncol = 1,
  rel_heights = c(1, .6, .6)
)
volcano_lfc_barplot
