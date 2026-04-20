# Inulin Aging Analysis

This directory contains exploratory analysis for a mouse study examining how lifelong inulin feeding relates to tissue metabolism across aging. The project includes both metabolomics and RNA-seq. We will focus on **metabolomics only** at first and **ignore `RNAseq/` for now**.

## Experimental Background

Inulin is a pre-biotic fiber that is becoming appreciated for benefits toward gut health, microbiome, and overall health. The positive effects of long-term inulin consumption have yet to be fully characterized. The broad biological goal in this project is to understand how long-term inulin exposure shapes host metabolism over the lifespan. In this analysis, the immediate emphasis is tissue metabolomics, with data present for **liver** and **muscle**. A practical working research question is:

**Which metabolites differ across aging in inulin-fed mice relative to control diet, and how are those patterns reflected across tissues?**

## Groups Being Analyzed

The main study groups are:

- **Age group:** `young` and `old`
- **Diet:** `CD` and `ID` (control or inulin diet)
- **Study cohorts:** `1`, `2`, `4`, and `5`
- **Tissue:** `musle` or `liver`

The group combinations by cohort are:

- cohort `1`: `old`, with both `CD` and `ID`
- cohort `2`: `old`, with both `CD` and `ID`
- cohort `4`: `young`, `CD`
- cohort `5`: `young`, `CD`


## Directory Orientation

The top level currently has two parallel projects:

- `Metabolomics `
- `RNAseq`

We will stay in `Metabolomics `.

Within `Metabolomics `:

- `Metabolomics/Metabolomics.Rproj`
  The project file. Open this first so `here()` resolves paths from the metabolomics project root.

- `Metabolomics/code`
  Data cleanup and analysis scripts live here.

- `Metabolomics/code/01 Process tissue metabolomics.R`
  Main processing script. It reshapes Compound Discoverer exports, filters annotations, resolves repeated/isomeric features, assigns cohorts, and writes processed metabolite tables.

- `Metabolomics/code/z-source.R`
  Shared helper functions used by the main script, including functions for resolving repeated/isomeric annotations.

- `Metabolomics/data/raw tissue`
  Raw Compound Discoverer Excel exports. Current files are liver and muscle runs in positive and negative mode.

- `Metabolomics/data/processed`
  Intended home for clean data tables and processed outputs.

- `Metabolomics/plots`
  Plot output area.

## Metabolomics Workflow, At A High Level

The current tissue metabolomics script `01 Process tissue metabolomics.R` does the following cleanup steps:

1. Read raw Compound Discoverer exports.
2. Convert/reformat tables into a long format with one ion-count observation per row.
3. Filter annotations.
4. Resolve repeated annotations and isomers within ion mode and cohort.
5. Compare annotation agreement across cohorts.
6. Consolidate positive and negative ion mode calls.
7. Internal standard and median normalization.

We will make subsequent scripts for statistical analysis and plotting

## Syntax And Housekeeping Conventions

- Open the metabolomics R project and use `here::here()` for paths.
- Keep `here::i_am("code/...")` near the top of scripts so the project root is explicit.
- Do **not** use absolute file paths tied to one computer.
- Do **not** modify raw files in `data/raw tissue/`.
- Write derived tables, cleaned data, and intermediate outputs to `data/processed/`.
- Put new analysis scripts in `Metabolomics /code/` and give them informative, ordered names.
- Keep reusable helper functions in `z-source.R` or in another clearly named helper file.


