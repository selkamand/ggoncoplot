
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggoncoplot <a href="https://selkamand.github.io/ggoncoplot/"><img src="man/figures/logo.png" align="right" height="104" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/ggoncoplot)](https://CRAN.R-project.org/package=ggoncoplot)
[![Codecov test
coverage](https://codecov.io/gh/selkamand/ggoncoplot/branch/master/graph/badge.svg)](https://app.codecov.io/gh/selkamand/ggoncoplot?branch=master)
[![R-CMD-check](https://github.com/selkamand/ggoncoplot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/selkamand/ggoncoplot/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**ggoncoplot** creates interactive oncoplots from mutation level
datasets

## Installation

You can install the development version of ggoncoplot like so:

``` r
remotes::install_github('selkamand/ggoncoplot')
```

## Usage

For complete usage, see
[manual](https://selkamand.github.io/ggoncoplot/articles/quick_start.html)

### Input

The input for ggoncoplot is a data.frame with 1 row per mutation in
cohort and columns describing the following:

- Gene Symbol

- Sample Identifier

- (optional) mutation type

- (optional) tooltip (character string: what we show on mouse hover over
  a particular mutation)

These columns can be in any order, and named anything. You define the
mapping of your input dataset columns to the required features in the
call to **ggoncoplot**

### Basic Example

``` r
library(ggoncoplot)

# TCGA GBM dataset from TCGAmuations package
gbm_csv <- system.file(package='ggoncoplot', "testdata/GBM_tcgamutations_mc3_maf.csv.gz")
gbm_df <- read.csv(file = gbm_csv, header=TRUE)

gbm_df |> 
  ggoncoplot(
    col_genes = 'Hugo_Symbol', 
    col_samples = 'Tumor_Sample_Barcode', 
    col_mutation_type = 'Variant_Classification', 
    topn = 10, 
    interactive = FALSE
  )
#> 
#> ── Identify Class ──────────────────────────────────────────────────────────────
#> ℹ Found 7 unique mutation types in input set
#> ℹ 0/7 mutation types were valid SO terms
#> ℹ 7/7 mutation types were valid MAF terms
#> ✔ Mutation Types are described using valid MAF terms ... using MAF palete
```

<img src="man/figures/README-example-1.png" width="100%" />

### Add marginal plots

``` r

gbm_df |> 
  ggoncoplot(
    col_genes = 'Hugo_Symbol', 
    col_samples = 'Tumor_Sample_Barcode', 
    col_mutation_type = 'Variant_Classification', 
    topn = 10, 
    draw_gene_barplot = TRUE, 
    draw_tmb_barplot = TRUE,
    interactive = FALSE
  )
#> 
#> ── Identify Class ──────────────────────────────────────────────────────────────
#> ℹ Found 7 unique mutation types in input set
#> ℹ 0/7 mutation types were valid SO terms
#> ℹ 7/7 mutation types were valid MAF terms
#> ✔ Mutation Types are described using valid MAF terms ... using MAF palete
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />
