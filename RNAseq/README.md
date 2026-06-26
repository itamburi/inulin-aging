# RNAseq Analysis

This directory contains an interactive R workflow for RNA-seq analysis in the inulin aging study. The current scripts focus on mouse gene-count processing, DESeq2 differential expression, pathway enrichment, and overlap figures.

Open `RNAseq.Rproj` before running scripts so `here::here()` resolves paths from this directory.

## Directory Layout

- `code/`
  Main analysis scripts, ordered by their numeric prefixes.
- `code/alignment_scripts/`
  Cluster submission scripts used for upstream alignment and count matrix generation.
- `data/raw/`
  Raw RSEM count matrices. Treat these as source data.
- `data/metadata/`
  Sample metadata used by the R analysis scripts.
- `data/extraction/`
  RNA extraction and Nanodrop records.
- `data/processed/counts/`
  Derived count matrices.
- `data/processed/deseq2/`
  DESeq2 result tables.
- `data/processed/pathway_enrichment/`
  GO GSEA, MSigDB Hallmark GSEA, ORA-related outputs, and curated pathway selections.
- `data/processed/overlaps/`
  DEG overlap summaries, gene lists, and manually flagged DEG tables.
- `data/processed/youth_shift/`
  DESeq2 contrast-concordance scores for genes that shift old ID toward or away from young CD relative to old CD.
- `data/processed/plsda/`
  PLS-DA sample scores and VIP-ranked gene tables.
- `data/processed/archive/`
  Older processed outputs kept for reference.
- `plots/`
  Generated analysis figures, organized by analysis type.
- `plots/qc/`
  Sample-level quality-control figures, including PCA plots.

## Current Script Order

1. `code/01_prepare_counts_and_diet_deseq2.R`
   Converts Ensembl gene IDs to mouse gene symbols, writes the processed count matrix, explores PCA/SVA, runs old ID vs old CD DESeq2 analyses, and generates pathway and volcano outputs.

2. `code/02_age_diet_gsea_overlap.R`
   Compares young CD vs old CD and old ID vs old CD, then generates GO BP/CC GSEA overlap outputs.

3. `code/03_deg_overlap_and_youth_shift.R`
   Builds DEG overlap summaries, Venn-style figures, and youth-like old ID shift concordance outputs from saved DESeq2 result tables.

4. `code/04_plsda_major_groups.R`
   Runs PLS-DA across young CD, old CD, and old ID samples using TMM-normalized log2 CPM values, then outputs ordination and top VIP gene heatmap figures.

## Main Inputs

- `data/raw/rsem_gene_counts_matrix.tsv`
- `data/raw/rsem_isoform_counts_matrix.tsv`
- `data/metadata/metadata.csv`
- `ensmust transcript ids and gene names.csv`

## Main Outputs

- `data/processed/counts/counts matrix with mouse gene symbols.csv`
- `data/processed/deseq2/deseq2 ID vs CD old only cohort adjusted.csv`
- `data/processed/deseq2/deseq2 young CD vs old CD.csv`
- `data/processed/deseq2/deseq2 old ID vs old CD for age overlap.csv`
- `data/processed/pathway_enrichment/go_gsea/go * gsea *.csv`
- `data/processed/pathway_enrichment/msigdb_hallmark/msigdb hallmark gsea *.csv`
- `data/processed/overlaps/deseq2 DEG venn overlap *.csv`
- `data/processed/youth_shift/deseq2 youth-like old ID shift *.csv`
- `data/processed/plsda/plsda major groups *.csv`
- `plots/qc/rnaseq pca unfiltered filtered counts.pdf`
- `plots/qc/rnaseq pca filtered counts major groups.pdf`
- `plots/plsda/plsda major groups *.pdf`
- `plots/youth_shift/deseq2 youth-like old ID shift *.pdf`
- `plots/deseq2/*.pdf`
- `plots/pathway_enrichment/**/*.pdf`
- `plots/overlaps/*.pdf`

## Workflow Notes

- Keep scripts linear and easy to step through interactively in RStudio.
- Keep `here::i_am("code/...")` near the top of analysis scripts.
- Do not modify files in `data/raw/`.
- Write derived tables to `data/processed/`.
- Write figures to `plots/`.
- Prefer clear section headers and explicit intermediate objects over dense helper abstractions.

## R Package Dependencies

The current scripts use packages from CRAN and Bioconductor, including:

- `tidyverse`
- `here`
- `edgeR`
- `limma`
- `sva`
- `DESeq2`
- `clusterProfiler`
- `org.Mm.eg.db`
- `EnsDb.Mmusculus.v79`
- `msigdbr`
- `enrichR`
- `cowplot`
- `ggforce`
