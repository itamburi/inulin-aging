
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

x = read.delim(here("rsem_gene_counts_matrix.tsv")) %>%
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
meta = read.csv(here("data/metadata.csv")) %>%
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
  left_join(., meta)

pc %>%
  ggplot(.,aes(PC1, PC2, color = diet) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw()
pc %>%
  ggplot(.,aes(PC1, PC2, color = age.grp) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw()
pc %>%
  ggplot(.,aes(PC1, PC2, color = factor(cohort)) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
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
  left_join(., meta %>% rownames_to_column(var="sample"))
pc %>%
  ggplot(.,aes(PC1, PC2, color = diet) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw()
pc %>%
  ggplot(.,aes(PC1, PC2, color = age.grp) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw()
pc %>%
  ggplot(.,aes(PC1, PC2, color = factor(cohort)) ) +
  geom_point() +
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw()
# not substantially different from unfiltered


# ========== 2.0 - SVA ==========
library(edgeR)
library(limma)
library(sva)

# Model matrix
# Full model: diet (main effect) + age group + cohort
mod  <- model.matrix(~ age.grp + diet, data = meta)
# Null model: exclude diet
mod0 <- model.matrix(~ age.grp, data = meta)

k_est <- num.sv(dge$counts, mod, method = "leek")
k_safe <- min(k_est, nrow(meta) - qr(mod)$rank - 1)
svobj <- svaseq(dat = dge$counts, mod = mod, mod0 = mod0, n.sv = k_safe)
dim(svobj$sv)


# As a data frame
sv_df <- as.data.frame(svobj$sv)
colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))

# Add to your meta / DGEList
dge$samples <- cbind(dge$samples, sv_df)
meta_with_SV <- cbind(meta, sv_df)

# Pairwise scatter of first few SVs
pairs(svobj$sv[, 1:min(5, ncol(svobj$sv))])

# Or check correlation with known covariates
cor(svobj$sv[,1], as.numeric(factor(meta$cohort)))

design <- model.matrix(~ age.grp + diet + SV1 + SV2, data = meta_with_SV[, c("age.grp", "diet", "SV1","SV2")])

# voom transform and fit
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extract DE for the group coefficient
dex = topTable(fit, coef = "dietID", number = Inf)
dex_signif = dex %>% filter(P.Value <= 0.05)


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

