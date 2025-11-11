# Prepare dataset for plotting

Take a dataframe containing a column describing sample IDs
(`col_sample`) Filter on `col_sample` %in% samples_to_show. Add any
missing samples_to_show not present DF as levels of `col_sample`. This
way, when plotting we can use scale_x_discrete(drop=FALSE) to display
all the samples we care about

## Usage

``` r
unify_samples(data, col_samples, samples_to_show)
```

## Arguments

- data:

  dataframe with a column describing sample IDs (data.frame)

- col_samples:

  name of column in `data` containing sample IDs (character)

- samples_to_show:

  the samples we want to show in plots. These samples should be the only
  ones represented in data.frame content, and any missing ones will be
  added as factor levels (character)

## Value

data.frame
