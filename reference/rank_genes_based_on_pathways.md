# Calculate Pathway-informed Genes Rankings

Which genes should appear at the top of the oncoplot? This function
takes pathway and gene ranks and returns a list of genes sorted first by
pathway then by gene rank. Gene & pathway rankings can be calculated
upstream. By default will use their order in gene_pathway_map.

## Usage

``` r
rank_genes_based_on_pathways(
  gene_pathway_map,
  generanks = unique(as.character(gene_pathway_map[[1]])),
  pathwayranks = unique(as.character(gene_pathway_map[[2]]))
)
```

## Arguments

- gene_pathway_map:

  dataframe where column 1 = gene names and column 2 = pathway names

- generanks:

  gene names in the order they should be ranked, where earlier in vector
  = further up in oncoplot. (character)

- pathwayranks:

  pathway names in the order they should be ranked, where earlier in
  vector = further up in oncoplot (character)

## Value

gene names, sorted based on order they should appear in oncoplot (first
= top). Only returns genes present in generanks (character)
