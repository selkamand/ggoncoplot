# Gene barplot

Gene barplot

## Usage

``` r
ggoncoplot_gene_barplot(
  data,
  fontsize_count = 14,
  palette = NULL,
  colour_mutation_type_unspecified = "grey10",
  show_axis,
  total_samples,
  show_genebar_labels = TRUE,
  genebar_label_nudge = 2,
  genebar_label_padding = 0.2,
  only_pad_if_labels_shown = TRUE,
  digits_to_round_to = 0,
  genebar_scale_n_breaks = 3,
  genebar_scale_breaks = ggplot2::waiver()
)
```

## Arguments

- data:

  data frame output by ggoncoplot_prep_df

- fontsize_count:

  fontsize of gene mutation count x axis (number)

- palette:

  a named vector mapping all possible mutation types (vector names) to
  colors (vector values, optional)

- colour_mutation_type_unspecified:

  colour of mutations in oncoplot and margin plots if
  `col_mutation_type` is not supplied (string)

- show_axis:

  show axis text/ticks/line (flag)

- total_samples:

  Strategy for calculating the total number of samples. This value is
  used to compute the proportion of mutation recurrence displayed in the
  tooltip when hovering over the gene barplot, or as a text annotation
  when `ggoncoplot_options(show_genebar_labels = TRUE)` is set to TRUE.
  Possible values:

  - **any_mutations**: All the samples that are in `data` (the mutation
    dataset), irrespective of whether they are on the oncoplot or not.

  - **oncoplot**: Only the samples that are present on the oncoplot.

  - **all**: All the samples in either `data` or `metadata`.

- show_genebar_labels:

  should gene barplot be labelled with % of samples the gene is mutated
  in (flag)

- genebar_label_nudge:

  how much padding to add between the gene barplot and bar annotations
  (number)

- genebar_label_padding:

  how much padding to add to the x axis of the gene barplot (number)

- only_pad_if_labels_shown:

  should expansion to x axis be applied if bar labels aren't shown?

- digits_to_round_to:

  how many digits to round recurrence proportions to

- genebar_scale_n_breaks:

  an integer guiding the number of breaks The algorithm may choose a
  slightly different number to ensure nice break labels. Will only have
  an effect if `genebar_scale_breaks = ggplot2::waiver()`. Use `NULL` to
  use the default

- genebar_scale_breaks:

  fine-grained control over the x axis breaks on the gene barplot. One
  of:

  - `NULL` for no minor breaks

  - `waiver()` for the default breaks (none for discrete, one minor
    break between each major break for continuous)

  - A numeric vector of positions

  - A function that given the limits returns a vector of minor breaks.
    When the function has two arguments, it will be given the limits and
    major break positions.

## Value

ggplot showing gene mutation counts
