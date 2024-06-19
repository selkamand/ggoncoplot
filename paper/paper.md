---
title: 'ggoncoplot: an R package for visualising somatic mutation data from cancer patient cohorts'
tags:
  - R
  - cancer
  - genomics
  - visualisation
  - oncoplot
authors:
  - name: Sam El-Kamand
    orcid: 0000-0003-2270-8088
    affiliation: 1 
  - name: Julian M.W. Quinn
    orcid: 0000-0001-9674-9646
    affiliation: "1, 2"
  - name: Mark J. Cowley
    affiliation: "1, 3"
    orcid: 0000-0002-9519-5714
    corresponding: true 
affiliations:
 - name: Children’s Cancer Institute, Australia
   index: 1
 - name: Skeletal Research Program, Garvan Institute of Medical Research, Australia
   index: 2
 - name: School of Clinical Medicine, UNSW Medicine & Health, Australia
   index: 3
date: 17 June 2024
bibliography: paper.bib
---

# Summary

The ggoncoplot R package generates interactive oncoplots (also called oncoprints) that visualize mutational patterns across patient cancer cohorts (\autoref{fig:oncoplot}). These plots reveal patterns of mutation co-occurrence in a cohort, with marginal plots that indicate correlations between gene mutations and specific tumour characteristics. These tumour characteristics can include user-supplied annotations such as patient clinical data, cancer subtypes and histological features.

ggoncoplot offers several features to enhance utility and user experience. These include automatic colour palette selection for common mutation impact dictionaries, customizable tooltips, and automatic rendering of clinical annotations as barplots or heatmaps for quantitative or qualitative data, respectively. ggoncoplot supports visualisation of mutation-level data in tidy, tabular formats, making it easy to run on existing large somatic mutation datasets stored as MAF files or in relational databases. The ggoncoplot package is available from github at <https://github.com/selkamand/ggoncoplot>.

![ggoncoplot output visualising mutational trends in the TCGA breast carcinoma cohort. Individual patient samples are plotted on the x-axis, ordered by ggoncoplot. The plot indicates (y-axis, sorted by genes mutation frequency) that PIK3CA is the most recurrently mutated gene, followed by TP53. Marginal plots indicate the total number of mutations per sample (top), and the number of samples showing mutations in each gene, coloured by mutation type (right). A range of clinical features, including progesterone and estrogen receptor status are shown on the marginal plot at the bottom. \label{fig:oncoplot}](oncoplot.pdf)

# Statement of Need

Oncoplots effectively visualize cohort-level mutation but are challenging to generate with the major R plotting systems (base, lattice, or ggplot2) due to their algorithmic and graphical complexity. Simplifying the process would make oncoplots more accessible to researchers. Packages like ComplexHeatmap [@Gu:2022], maftools [@Mayakonda:2018], and genVisR [@Skidmore:2016] make static oncoplots easier to create, but there remains a significant and unmet need to easily create oncoplots with:

-   **Interactive features**: Customizable tooltips, cross-selection of samples across different plots, and auto-copying of sample identifiers on click.

-   **Support for tidy datasets**: Compatibility with tidy, tabular mutation-level formats (MAF files or relational databases), typical of cancer cohort datasets.

-   **Auto colouring**: Automatic selection of color palettes for datasets with consequence annotations aligned with standard variant effect dictionaries (PAVE, SO, or MAF).

-   **Versatility**: The ability to visualize entities beyond gene mutations, including noncoding features (e.g., enhancers) and non-genomic entities (e.g., microbial presence in microbiome datasets).

We developed ggoncoplot as the first R package that addresses all these challenges simultaneously (\autoref{fig:comparison}). Examples of all key features are available in the [ggoncoplot manual](https://selkamand.github.io/ggoncoplot/articles/manual.html).


![Comparison of R packages for creating oncoplots. ^1^Requires the shiny and interactiveComplexHeatmap packages. ^2^Requires the user to first summarise mutations at the gene level as a sample by gene character matrix with mutations separated by semicolons (wide format). ^3^For MAF inputs the most severe consequence is chosen, however for non-MAF datasets users must manually define the mutation impact hierarchy. ^4^Non-unique mutation types are treated as one observation, however if different mutation types affect one gene, the indiviual mutations can be plotted with different shapes/configurations in a user-configured manner. \label{fig:comparison}](ggoncoplot_comparision.pdf)


# Acknowledgements

We thank the developers of the packages integral to ggoncoplot, especially David Gohel for ggiraph [@gohel:2024], which enables its interactivity, and Thomas Lin Pedersen for patchwork [@pedersen:2024] and ggplot2 maintenance. We also acknowledge Hadley Wickham and all contributors to ggplot2 [@wickham:2016]. 
Additionally, we thank Dr. Marion Mateos for her insightful feedback during the development of ggoncoplot.

# References
