# Generate score based on genes

Score used to sort samples based on which genes are mutated. Make sure
to run on one sample at once (use grouping)

## Usage

``` r
score_based_on_gene_rank(
  mutated_genes,
  genes_informing_score,
  gene_rank,
  debug_mode = FALSE
)
```

## Arguments

- mutated_genes:

  vector of genes that are mutated for a single sample (character)

- genes_informing_score:

  which genes determine the sort order? (character)

- gene_rank:

  what is the order of importance of genes used to determine sort order.
  Higher number = higher in sort order (character)

- debug_mode:

  debug mode (flag)

## Value

a score (higher = should be higher in the sorting order) (number)

## Examples

``` r
if (FALSE) { # \dontrun{
# First set of genes has a high rank since both BRCA2 and EGFR are mutated
score_based_on_gene_rank(c("TERT", "EGFR", "PTEN", "BRCA2"), c("EGFR", "BRCA2"), gene_rank = 1:2)

# If EGFR is mutated without BRCA2, we get a lower score
score_based_on_gene_rank(c("TERT", "EGFR", "PTEN", "IDH1"), c("EGFR", "BRCA2"), gene_rank = 1:2)

# If BRCA2 is mutated without EGFR,
# we get a score lower than BRCA2+EGFR but higher than EGFR alone due to higher gene_rank of BRCA2
score_based_on_gene_rank(c("TERT", "IDH1", "PTEN", "BRCA2"), c("EGFR", "BRCA2"), gene_rank = 1:2)
} # }
```
