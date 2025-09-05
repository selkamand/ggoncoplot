---
title: 'ggoncoplot: an R package for interactive visualisation of somatic mutation data from cancer patient cohorts'
tags:
- R
- cancer
- genomics
- visualisation
- oncoplot
date: "17 June 2024"
output: word_document
authors:
- name: "Sam El-Kamand"
  orcid: "0000-0003-2270-8088"
  affiliation: 1
- name: Julian M.W. Quinn
  orcid: "0000-0001-9674-9646"
  affiliation: 1
- name: Mark J. Cowley
  affiliation: 1, 2
  orcid: "0000-0002-9519-5714"
  corresponding: true
bibliography: paper.bib
affiliations:
- name: Children’s Cancer Institute, Australia
  index: 1
- name: School of Clinical Medicine, UNSW Medicine & Health, Australia
  index: 2
---

# Summary

The ggoncoplot R package generates interactive oncoplots to visualize mutational patterns across patient cancer cohorts (\autoref{fig:oncoplot}). Oncoplots, also called oncoprints, reveal patterns of gene co-mutation and include marginal plots that indicate co-occurrence of gene mutations with tumour and clinical features. It is useful to relate gene mutation patterns seen in an oncoplot to patterns in other plot types, including gene expression t-SNE plots or methylation UMAPs. The simplest and most intuitive approach to examining such relations is to link plots dynamically such that samples selected in an oncoplot can be highlighted in other plots, and vice versa. There are, however, no existing oncoplot-generating R packages that support dynamic data linkage between different plots. To address this gap and enable rapid exploration of a variety of data types we constructed the ggoncoplot package for the production of oncoplots that are easily integrated with custom visualisations and that support synchronised data-selections across plots (\autoref{fig:multimodal_selection}). The intended audience for ggoncoplot includes genomics researchers, clinical data analysts, and computational biologists studying the mutation patterns in cancer cohort datasets. ggoncoplot is available on GitHub at <https://github.com/selkamand/ggoncoplot>. 

![ggoncoplot output visualising mutational trends in the TCGA breast carcinoma cohort. Individual patient samples are plotted on the x-axis, hierarchically sorted so that samples with the most frequent gene mutations appear on the leftmost side. The plot indicates that PIK3CA is the most frequently mutated gene, followed by TP53. Marginal plots indicate the total number of mutations per sample (top), and the number of samples showing mutations in each gene, coloured by mutation type (right). A range of clinical features, including progesterone and estrogen receptor status are shown on the marginal plot at the bottom. A detailed description of the ggoncoplot sorting algorithm is available [here](https://selkamand.github.io/ggoncoplot/articles/sorting_algorithm.html). \label{fig:oncoplot}](oncoplot.pdf)


![Example of the ggoncoplot shown in Figure 1, where the oncoplot has been dynamically cross-linked to a gene expression t-SNE plot (top left) and a methylation UMAP (top right). Here, the lasso tool was used to select a cluster of gene expression data points (i.e., individual samples) in the t-SNE plot. Selected samples were automatically highlighted on the UMAP and oncoplot. This reveals that samples which cluster on the left of the t-SNE plot also cluster in the oncoplot, chiefly containing mutations in TP53 and wild type PIK3CA. The plots of progesterone, estrogen, HER2 status and triple negative classification show that the samples selected in the t-SNE are enriched for triple negative breast cancers.  In contrast to the oncoplot, the methylation UMAP shows no strong clustering, consistent with knowledge of methylation patterns in triple negative breast cancer. \label{fig:multimodal_selection}](multimodal_selection_with_lasso.png)



# Statement of Need

Oncoplots are highly effective for visualising mutation data in cancer cohorts but are challenging to generate with the major R plotting systems (base, lattice, or ggplot2) due to their algorithmic and graphical complexity.  Simplifying the process of generating oncoplots would make them more accessible to researchers. Existing packages including ComplexHeatmap [@Gu:2022], maftools [@Mayakonda:2018], and genVisR [@Skidmore:2016] all make static oncoplots easier to create, but there is still a significant unmet need for a user-friendly method of creating oncoplots with the following features: 

-	**Interactive plots**: Customizable tooltips, cross-selection of samples across different plots, and auto-copying of sample identifiers on click. This enables exploration of multiomic datasets as shown in \autoref{fig:multimodal_selection}.

-	**Support for tidy datasets**: Compatibility with tidy, tabular mutation-level formats that cancer cohort datasets are typically stored in. This greatly improves the range of datasets that can be quickly and easily visualised in an oncoplot since genomic data in Mutation Annotation Format (MAF) files and relational databases usually follow this structure.

-	**Auto-colouring**: Automatic selection of accessible colour palettes for datasets where the consequence annotations are aligned with standard variant effect dictionaries including Prediction and Annotation of Variant Effects (PAVE), Sequence Ontology (SO) and MAF Variant Classifications.

- **Versatility**: The ability to visualize entities other than gene mutations, such as noncoding features (e.g., promoter or enhancer mutations) and non-genomic entities (e.g., microbial presence in microbiome datasets). 


We developed ggoncoplot as the first R package to address all these challenges together (\autoref{fig:comparison}). Examples of all key features are available in the [ggoncoplot manual](https://selkamand.github.io/ggoncoplot/articles/manual.html).


![Comparison of R packages for creating oncoplots. ^1^Requires the shiny and interactiveComplexHeatmap packages. ^2^Exclusively colours tiles based on mutation impact which must be described using valid MAF variant classification terms. ^3^Requires the user to first summarise mutations at the gene level and format as a sample by gene matrix with mutations separated by semicolons (wide format). ^4^For MAF inputs the most severe consequence is chosen, however for non-MAF datasets users must manually define the mutation impact hierarchy. ^5^Non-unique mutation types are treated as one observation, however if different mutation types affect one gene, the individual mutations can be plotted with different shapes or sizes in a user-configured manner. \label{fig:comparison}](ggoncoplot_comparison.pdf)

# Acknowledgements

The results shown here are in whole or part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga. Methylation, expression, and somatic mutation datasets were obtained from the Xena TCGA Pan-Cancer Atlas Hub [@goldman:2020; @ellrott:2018]. 

We thank the developers of the packages integral to ggoncoplot, especially David Gohel for ggiraph [@gohel:2024], which enables its interactivity, and Thomas Lin Pedersen for patchwork [@pedersen:2024] and ggplot2 maintenance. We also acknowledge Hadley Wickham and all contributors to ggplot2 [@wickham:2016]. 
Additionally, we thank Dr. Marion Mateos for her insightful feedback during the development of ggoncoplot.

# References
