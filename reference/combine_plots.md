# Combine margin plots with main plot

Combine margin plots with main plot

## Usage

``` r
combine_plots(
  gg_main,
  gg_tmb = NULL,
  gg_gene = NULL,
  gg_metadata = NULL,
  gg_tmb_height,
  gg_gene_width,
  gg_metadata_height,
  metadata_position,
  buffer_metadata,
  buffer_tmb
)
```

## Arguments

- gg_main:

  main oncoplot tileplot (ggplot)

- gg_tmb:

  barplot describing total mutations. Set to NULL to not draw barplot
  (ggplot)

- gg_gene:

  barplot describing number of mutated samples per gene. Set to NULL to
  not draw barplot (ggplot)

- gg_metadata:

  tile plot describing sample-level metadata

- gg_tmb_height:

  percentage of plot height taken up by TMB plot (should be between
  5-95) (number)

- gg_gene_width:

  percentage of plot width taken up by genebar plot (should be between
  5-95) (number)

- gg_metadata_height:

  percentage of plot height taken up by metadata plot (should be between
  5-95) (number)

- metadata_position:

  should metadata plot be on the 'top' or the 'bottom' of the oncoplot?

- buffer_metadata, buffer_tmb:

  amount of space to add between the main oncoplot and tmb/metadata
  marginal plots (number)

## Value

patchwork object (or ggplot obj if both `gg_tmb` and `gg_gene` are NULL)
