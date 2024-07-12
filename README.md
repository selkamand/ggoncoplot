
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggoncoplot <a href="https://selkamand.github.io/ggoncoplot/"><img src="man/figures/logo.png" align="right" height="104"/></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/ggoncoplot)](https://CRAN.R-project.org/package=ggoncoplot)
[![Codecov test
coverage](https://codecov.io/gh/selkamand/ggoncoplot/branch/master/graph/badge.svg)](https://app.codecov.io/gh/selkamand/ggoncoplot?branch=master)
[![R-CMD-check](https://github.com/selkamand/ggoncoplot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/selkamand/ggoncoplot/actions/workflows/R-CMD-check.yaml)
![GitHub Issues or Pull
Requests](https://img.shields.io/github/issues-closed/selkamand/ggoncoplot)
[![](https://img.shields.io/github/languages/code-size/selkamand/ggoncoplot.svg)](https://github.com/selkamand/ggoncoplot)
![GitHub last
commit](https://img.shields.io/github/last-commit/selkamand/ggoncoplot)

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
[manual](https://selkamand.github.io/ggoncoplot/articles/manual.html)

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
#> ── Identify Class ──
#> 
#> ℹ Found 7 unique mutation types in input set
#> ℹ 0/7 mutation types were valid PAVE terms
#> ℹ 0/7 mutation types were valid SO terms
#> ℹ 7/7 mutation types were valid MAF terms
#> ✔ Mutation Types are described using valid MAF terms ... using MAF palete
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

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
#> ── Identify Class ──
#> 
#> ℹ Found 7 unique mutation types in input set
#> ℹ 0/7 mutation types were valid PAVE terms
#> ℹ 0/7 mutation types were valid SO terms
#> ℹ 7/7 mutation types were valid MAF terms
#> ✔ Mutation Types are described using valid MAF terms ... using MAF palete
#> ! TMB plot: Ignoring `col_mutation_type` since `log10_transform = TRUE`.
#> This is because you cannot accurately plot stacked bars on a logarithmic scale
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

### Add clinical metadata

``` r
gbm_clinical_csv <- system.file(package = "ggoncoplot", "testdata/GBM_tcgamutations_mc3_clinical.csv")
gbm_clinical_df <- read.csv(file = gbm_clinical_csv, header = TRUE)

gbm_df |> 
  ggoncoplot(
   col_genes = "Hugo_Symbol",
   col_samples = "Tumor_Sample_Barcode",
   col_mutation_type = "Variant_Classification",
   metadata = gbm_clinical_df,
   cols_to_plot_metadata = c('gender', 'histological_type', 'prior_glioma', 'tumor_tissue_site'),
   draw_tmb_barplot = TRUE, 
   draw_gene_barplot = TRUE, 
   show_all_samples = TRUE,
   interactive = FALSE
  )
#> ℹ 2 samples with metadata have no mutations. Fitering these out
#> ℹ To keep these samples, set `metadata_require_mutations = FALSE`. To view them in the oncoplot ensure you additionally set `show_all_samples = TRUE`
#> → TCGA-06-0165-01
#> → TCGA-06-0167-01
#> 
#> ── Identify Class ──
#> 
#> ℹ Found 7 unique mutation types in input set
#> ℹ 0/7 mutation types were valid PAVE terms
#> ℹ 0/7 mutation types were valid SO terms
#> ℹ 7/7 mutation types were valid MAF terms
#> ✔ Mutation Types are described using valid MAF terms ... using MAF palete
#> ! TMB plot: Ignoring `col_mutation_type` since `log10_transform = TRUE`.
#> This is because you cannot accurately plot stacked bars on a logarithmic scale
#> 
#> ── Plotting Sample Metadata ────────────────────────────────────────────────────
#> ! Categorical columns must have <= 6 unique values to be visualised. Columns with too many unique values:  (20),  (388), and  (388)
#> 
#> ── Sorting
#> ℹ Sorting X axis by: Order of appearance
#> 
#> ── Generating Plot
#> ℹ Found 4 plottable columns in data
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

## Acknowledgements

We acknowledge the developers and contributors whose packages and
efforts were integral to the development of ggoncoplot:

- **David Gohel** for the `ggiraph` package, which enables the
  interactivity of ggoncoplot.
- **Thomas Lin Pedersen** for his contributions to the `patchwork`
  package and the maintenance of `ggplot2`.
- **Hadley Wickham** and all contributors to the `ggplot2` package,
  which provides a robust foundation for data visualization in R.

Additionally, we thank **Dr. Marion Mateos** for her insightful feedback
during the early stages of ggoncoplot development.

## Community Contributions

### Contribute to the Software

There are many ways to contribute to ggoncoplot.

1.  Request features you would like to by [creating new issues on
    github](https://github.com/selkamand/ggoncoplot/issues)
2.  [Make your visualisation packages
    ggoncoplot-compatible](#make-your-visualisation-packages-ggoncoplot-compatible)
3.  [Directly contribute to the ggoncoplot
    codebase](#directly-contribute-to-the-ggoncoplot-codebase)

#### Make your visualisation packages ggoncoplot-compatible

If your package produces ggplots that you would like to interactively
link with ggoncoplot, consider converting your geoms to their ggiraph
interactive equivalents and adding a data_id based on a sample
identifier column in the dataset. That way end-users can create a
data-linked oncoplot composed with your packages plots using patchwork
([example](https://selkamand.github.io/ggoncoplot/articles/manual.html#interaction-with-other-packages))

#### Directly contribute to the ggoncoplot codebase

We welcome contributions from the community to enhance and expand the
functionality of `ggoncoplot`. Whether you want to fix a bug, add new
features, improve documentation, or optimize performance, your efforts
are highly valued. To get started:

1.  **Fork the Repository**: Click on the ‘Fork’ button at the top right
    of this page to create a copy of the repository in your GitHub
    account.

2.  **Clone the Repository**: Use `git clone` to clone your forked
    repository to your local machine.

    ``` bash
    git clone https://github.com/selkamand/ggoncoplot.git
    ```

3.  **Create a Branch**

    ``` bash
    git checkout -b feature-name
    ```

4.  **Make Changes**: Implement your changes in the new branch

5.  **Commit and Push**: Commit your changes and push the branch to your
    forked repository.

6.  **Create a Pull Request**: Go to the original repository and open a
    pull request from your branch. Please provide a clear description of
    your changes and any relevant issues or discussions.

### Report Issues or Problems with the Software

If you encounter any issues, bugs, or have suggestions for improvements,
please report them using the [GitHub Issues
Tab](https://github.com/selkamand/ggoncoplot/issues/).

### Seek Support

For any questions or support regarding the use of ggoncoplot you can:

- **Check the Documentation**: Comprehensive documentation is available
  [here](https://selkamand.github.io/ggoncoplot/index.html).

- **Create a** [new
  issue](https://github.com/selkamand/ggoncoplot/issues/new) with your
  query.

- **Browse Existing Issues**: Check the
  [Issues](https://github.com/selkamand/ggoncoplot/issues) page to see
  if your query has been addressed.

- **Contact Us**: If you need direct assistance, please [contact the
  maintainers directly](mailto:selkamand@ccia.org.au?subject=ggoncoplot)
