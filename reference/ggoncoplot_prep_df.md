# Prep data for oncoplot

Prep data for oncoplot

## Usage

``` r
ggoncoplot_prep_df(
  data,
  col_genes,
  col_samples,
  genes_for_oncoplot,
  col_mutation_type = NULL,
  col_tooltip = col_samples,
  pathway = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  data for oncoplot. A data.frame with 1 row per mutation in your
  cohort. Must contain columns describing gene_symbols and
  sample_identifiers (data.frame)

- col_genes:

  name of **data** column containing gene names/symbols (string)

- col_samples:

  name of **data** column containing sample identifiers (string)

- genes_for_oncoplot:

  a list of genes to include in the oncoplot (character).

- col_mutation_type:

  name of **data** column describing mutation types (string)

- col_tooltip:

  name of **data** column containing whatever information you want to
  display in (string)

- pathway:

  a two column dataframe describing pathway. The column containing gene
  names should have the same name as **col_gene** (data.frame, optional)

- verbose:

  verbose mode (flag, default TRUE)

## Value

dataframe with the following columns: 'Gene', 'Sample', 'MutationType',
'Tooltip'. Sample is a factor with levels sorted in appropriate order
for oncoplot vis. Genes represents either topn genes or specific genes
set by `genes_to_include`

## Examples

``` r
#' # ===== GBM =====
gbm_csv <- system.file(
  package = "ggoncoplot",
  "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
)

gbm_df <- read.csv(file = gbm_csv, header = TRUE)

# Get genes in appropriate order for oncoplot
genes_for_oncoplot <- ggoncoplot:::get_genes_for_oncoplot(
  data = gbm_df,
  col_samples = "Tumor_Sample_Barcode",
  col_genes = "Hugo_Symbol",
  topn = 20,
  verbose = FALSE
)

# Create dataframe basis of oncoplot (1 row per sample-gene combo)
ggoncoplot:::ggoncoplot_prep_df(
  gbm_df,
  col_genes = "Hugo_Symbol",
  col_samples = "Tumor_Sample_Barcode",
  col_mutation_type = "Variant_Classification",
  genes_for_oncoplot = genes_for_oncoplot
)
#> # A tibble: 1,057 × 5
#>    Sample          Gene   MutationType      MutationCount Tooltip        
#>    <fct>           <fct>  <chr>                     <int> <chr>          
#>  1 TCGA-19-5956-01 PTEN   Multi_Hit                     3 TCGA-19-5956-01
#>  2 TCGA-19-5956-01 TP53   Multi_Hit                     2 TCGA-19-5956-01
#>  3 TCGA-19-5956-01 TTN    Multi_Hit                    28 TCGA-19-5956-01
#>  4 TCGA-19-5956-01 EGFR   Missense_Mutation             1 TCGA-19-5956-01
#>  5 TCGA-19-5956-01 MUC16  Multi_Hit                    15 TCGA-19-5956-01
#>  6 TCGA-19-5956-01 FLG    Multi_Hit                     8 TCGA-19-5956-01
#>  7 TCGA-19-5956-01 NF1    Multi_Hit                     3 TCGA-19-5956-01
#>  8 TCGA-19-5956-01 RYR2   Missense_Mutation             1 TCGA-19-5956-01
#>  9 TCGA-19-5956-01 ATRX   Multi_Hit                     5 TCGA-19-5956-01
#> 10 TCGA-19-5956-01 PIK3R1 Multi_Hit                     2 TCGA-19-5956-01
#> # ℹ 1,047 more rows
```
