# Simulated Cancer Dataset Metadata

A sample‐level metadata table for the `oncosim` simulated cancer
dataset. Contains assorted numeric, categorical, clinical, and logical
features for each sample.

## Usage

``` r
oncosim_metadata
```

## Format

### `oncosim_metadata`

A data frame with 11 rows and 6 columns:

- Samples:

  Unique sample identifiers

- numeric_feature:

  Numeric variable including zeros, positive and negative values, `NA`,
  and `Inf`/`-Inf`

- categorical_feature4levels:

  Categorical variable with four levels (`"cat"`, `"dog"`, `"magpie"`,
  `"giraffe"`), may contain empty strings or `NA`

- clinical_feature2levels:

  Clinical categorical variable indicating biological sex with two
  levels (`"male"`, `"female"`), may contain `NA`

- logical_feature:

  Logical variable with `TRUE`, `FALSE`, or `NA`

- numeric_that_could_be_logical:

  Integer variable coded as `0` or `1` (and `NA`) that could be
  interpreted as logical

## Source

not applicable, simulated
