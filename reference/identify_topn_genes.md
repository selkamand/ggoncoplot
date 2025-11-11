# Identify top genes from a mutation df

Identify top genes from a mutation df

## Usage

``` r
identify_topn_genes(
  data,
  col_samples,
  col_genes,
  topn,
  genes_to_ignore = NULL,
  return_extra_genes_if_tied = FALSE,
  verbose = TRUE
)
```

## Arguments

- data:

  data for oncoplot. A data.frame with 1 row per mutation in your
  cohort. Must contain columns describing gene_symbols and
  sample_identifiers (data.frame)

- col_samples:

  name of **data** column containing sample identifiers (string)

- col_genes:

  name of **data** column containing gene names/symbols (string)

- topn:

  how many of the top genes to visualize. Ignored if `genes_to_include`
  is supplied (number, default 10)

- genes_to_ignore:

  names of the genes that should be ignored (character, optional)

- return_extra_genes_if_tied:

  instead of strictly returning `topn` genes, in the case of ties (where
  multiple genes are mutated in the exact same number of samples,
  complicating selection of top n genes), return all tied genes
  (potentially more than topn). If FALSE, will return strictly `topn`
  genes, breaking ties based on order of appearance in dataset (flag,
  default FALSE)

- verbose:

  verbose mode (flag, default TRUE)

## Value

vector of topn genes. Their order will be their rank (most mutated =
first) (character)
