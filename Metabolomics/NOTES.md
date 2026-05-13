# Notes

## Setup Coding Tools

Follow these steps to open the project and start working with the analysis scripts.

1. In Finder, navigate through Google Drive to the `Inulin Aging/` directory.
2. Open `Metabolomics.Rproj` in the `Metabolomics/` folder to launch RStudio. This sets the project working directory correctly, so relative file paths in the scripts resolve from the project root.
3. In RStudio, open the analysis script you want to work on from the `code/` folder.
4. In Finder, right-click the `Metabolomics` directory and select **Open New Terminal at Folder** to launch a terminal in this project directory.
5. In that terminal window, type `codex` to open the Codex command-line interface.

## Housekeeping

- `AGENTS.md` is the project instructions file for Codex. It explains the workflow conventions for this repository, including how the R scripts should be organized and what coding style to follow.
- `README.md` is the top-level project overview file. If one is added, it should give a short summary of the project, explain the directory structure, and list the main steps for running the analysis.
- In the `code/` folder, keep analysis work organized into clear sections with headings and comments. Each major step should be separated so someone can run and inspect the script one stage at a time.

## Next Step

Once the setup is ready, a good first task is to ask Codex to make a volcano plot. Be specific about what you want the plot to show, how it should be built, and what should be labeled.

Example prompt:

> In the next section, make a volcano plot using `ggplot`. Make the x axis the log2 fold change and the y axis the -log10 p value. Add p value cutoffs, draw threshold lines, and label the points that pass a chosen significance threshold.

If needed, add more detail to the prompt, such as which comparison to use, which groups to compare, whether to color significant points differently, and how many points should be labeled.

## Upcoming

- Make a PCA for the untargeted metabolomics data, separated by tissue.
- Correlate RNA-seq genes with metabolites.
- Correlate metabolites with cytokines.
- Look for other cross-modal correlations that may be biologically interesting.
