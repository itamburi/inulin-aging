# AGENTS.md

## Project Guidance

This repository is primarily an interactive R analysis workflow for metabolomics data.

## Coding Style

- Prefer linear, top-to-bottom R scripts that are easy to run interactively.
- Minimize helper functions unless they clearly improve readability without hiding important analysis logic.
- Keep code human-readable and easy to follow for someone stepping through the script manually.
- Prefer `dplyr` and broader `tidyverse` conventions for data manipulation.
- Favor explicit transformation steps over dense abstraction.
- Use clear section headers and comments to separate major analysis stages.

## Script Expectations

- Assume scripts may be run interactively in RStudio rather than only as fully automated pipelines.
- Keep intermediate objects and processing steps visible when that improves traceability.
- Avoid unnecessary indirection, metaprogramming, or framework-style structure.
- When adding code, optimize for clarity and inspection over reuse for its own sake.

## File Organization

- Main analysis scripts live in `code/`.
- Preserve the existing script-oriented workflow unless there is a strong reason to change it.
