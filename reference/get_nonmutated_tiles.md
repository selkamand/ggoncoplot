# Get data.frame o

Takes same data input as ggoncoplot and returns a dataframe with
'Sample' and 'Gene' columns ONLY for sample-gene pairs that are
unmutated. This lets us colour render them separately (as grey)

## Usage

``` r
get_nonmutated_tiles(data)
```

## Arguments

- data:

  transformed data from
  [`ggoncoplot_prep_df()`](https://selkamand.github.io/ggoncoplot/reference/ggoncoplot_prep_df.md)
  (data.frame)

## Value

a dataframe with 'Sample' and 'Gene' columns ONLY for sample-gene pairs
that are unmutated. This lets us colour render them separately (as grey)
(data.frame)
