# Differences in pathogen resistance between diploid and polyploid plants: a systematic review and meta-analysis
---

This repository contains the data and code used in our meta-analysis. All analyses are contained in Hagen_Mason_2023_code.R, which must be run using R. Running the code requires two inputs included in this repository: the data file containing all information about included effect sizes, HagenandMason2023_data.xlsx, and the phylogeny file, family_level_angiosperm_tree.txt. The Excel file possesses two sheets: Sheet 1, which contains the data necessary for running analyses, and Sheet 2 ("Column Explanations"), which contains information about each column in Sheet 1. Finally, the repository includes a bibliography of every paper included in our meta-analysis, included_papers_biblio.docx.

## Description of the data and file structure
All effect sizes contained in the Excel file were derived from published papers. Information about those publications can be found in included_papers_biblio.docx, and information about the procedures used in the R code file can be found in the methods of our meta-analysis.

##Sharing/Access Information
Data were derived from the following sources:
	*Phylogeny file: Qian and Zhang 2014: "Using an updated time-calibrated family-level phylogeny of seed plants to test for non-random patterns of life forms across the phylogeny"
	*Effect sizes: Papers listed in included_papers_biblio.docx

## Code/Software
Prior to running any analyses, users must write in the directory locations of the necessary input files (the effect size data and the phylogeny) at Lines 18 and 19 in Hagen_Mason_2023_code.R. In our original meta-analysis, we used R versions 4.0.3 and 4.2.1 to run analyses, as well as the following package versions: metafor version 3.8-1, metagear version 0.7, dplyr version 1.0.10, ape version 5.6-2, ggplot2 versions 3.3.6 and 3.4.0, ggpubr versions 0.4.0 and 0.5.0, emmeans version 1.8.1-1, xlsx version 0.6.5, orchaRd (downloaded from daniel1noble/orchaRd), and rlang version 1.0.6.
