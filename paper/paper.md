---
title: 'ggoncoplot: an R package for visualisation of somatic mutation data from cancer patient cohorts  '
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

ggoncoplot offers several features to enhance utility and user experience. These include automatic colour palette selection for common mutation impact dictionaries, customizable tooltips, and automatic rendering of clinical annotations as barplots or heatmaps for quantitative or qualitative data, respectively. ggoncoplot supports visualisation of mutation-level data in tidy, tabular formats, making it easy to run on existing large somatic mutation datasets stored as MAF files or in relational databases. The ggoncoplot package is available from github at https://github.com/selkamand/ggoncoplot. 

![Default ggoncoplot output visualising mutational trends in the TCGA glioblastoma multiforme cohort. Individual patient samples are plotted on the x-axis, ordered by ggoncoplot.The plot indicates (y-axis, sorted by genes mutation frequency) that PTEN is the most recurrently mutated gene in the cohort, followed by TP53. Marginal plots indicate the total number of mutations per sample (top), and the number of samples showing mutations in each gene, coloured by mutation type (right). A range of clinical features, including histological type and gender, are shown on the marginal plot at the bottom. This default settings output can be altered as required. \label{fig:oncoplot}](oncoplot.png)

# Statement of Need

Oncoplots effectively visualize cohort-level mutation but are challenging to generate with the major R plotting systems (base, lattice, or ggplot2) due to their algorithmic and graphical complexity. Simplifying the process would make oncoplots more accessible to researchers. Packages like ComplexHeatmap [@Gu:2022], maftools [@Mayakonda:2018], and genVisR [@Skidmore:2016] make static oncoplots easier to create, but there remains a significant and unmet need to easily create oncoplots with:

### Interactive Features
- Customizable tooltips
- Cross-selection of samples across different plots
- Auto-copying of sample identifiers on click

### Support for Tidy Datasets
- Compatibility with tidy, tabular mutation-level formats (MAF files or relational databases), typical of cancer cohort datasets

### Auto Colouring
- Automatic selection of color palettes for datasets with consequence annotations aligned with standard variant effect dictionaries (PAVE, SO, or MAF)

### Versatility
- The ability to visualize entities beyond gene mutations, including noncoding features (e.g., enhancers) and non-genomic entities (e.g., microbial presence in microbiome datasets)

We developed ggoncoplot as the first R package that addresses all these challenges simultaneously (Table 1). Examples of all key features are available in the [ggoncplot manual](https://selkamand.github.io/ggoncoplot/articles/manual.html).

# Table 1. Comparison of R packages for creating oncoplots

| Property                                                                                        | complexheatmap                              | maftools           | GenVisR                                                     | ggoncoplot             |
| ----------------------------------------------------------------------------------------------- | ------------------------------------------- | ------------------ | ----------------------------------------------------------- | ---------------------- |
| Sample sorting algorithm                                                                        | memo sort                                   | hierarchical sort  | hierarchical sort                                           | hierarchical sort      |
| Underlying plotting system                                                                      | BaseR                                       | BaseR              | ggplot2                                                     | ggplot2                |
| Automatic rendering of clinical annotations as bar or tile plots based on datatype              | No                                          | No                 | No                                                          | Yes                    |
| Works on tabular, tidy, long-form input data as would be stored in large databases              | No                                          | Yes                | Yes                                                         | Yes                    |
| Interactive                                                                                     | Yes<sup>1</sup>                             | No                 | No                                                          | Yes                    |
| Customisable tooltips                                                                           | No                                          | No                 | No                                                          | Yes                    |
| Allows any mutation dictionary to be used?                                                      | Yes                                         | No                 | Yes                                                         | Yes                    |
| Automatic colour palette selection when mutation impact dictionary conforms to known ontologies | No                                          | Yes (MAF only)     | No                                                          | Yes (MAF, SO, or PAVE) |
| Approach for resolving genes with multiple mutations                                            | Different Visualisation on Plot<sup>3</sup> | Flags as Multi-Hit | Picks more severe consequence or leaves to user<sup>3</sup> | Flags as Multi-Hit     |
| Supports a mutation level dataset as input                                                      | No<sup>4</sup>                              | Yes                | Yes                                                         | Yes                    |
| Native Support for Faceting by Pathway                                                          | No                                          | Yes                | No                                                          | Yes                    |
| Supports marginal plots describing TMB, gene mutation recurrence, and clinical annotations      | Yes                                         | Yes                | Yes                                                         | Yes                    |

<sup>1</sup> Can be made interactive and displayed in a shiny app using the interactiveComplexHeatmap package.  
<sup>2</sup> Mutations must already be summarised at the gene level. Expects a sample X genes character matrix with different mutations separated by semicolons.  
<sup>3</sup> If a MAF is supplied it will choose the most severe consequence. For non-MAF dataset user can choose to define the mutation impact hierarchy.  
<sup>4</sup> If multiple mutations are of different types, it can be rendered in different ways on plot (user-controlled) - if identical, non-unique mutation types are treated as one observation.


# Acknowledgements

We extend our gratitude to the developers of the packages integral to ggoncoplot. We owe special thanks to David Gohel for his work on ggiraph [@gohel:2024], which enables the interactivity of ggoncoplot, and to Thomas Lin Pedersen for his contributions to patchwork [@pedersen:2024] and the maintenance of ggplot2. We also acknowledge Hadley Wickham and all contributors to the ggplot2 package [@wickam:2016]. Additionally, we thank Dr. Marion Mateos for her insightful feedback during the early stages of ggoncoplot development. 

# References