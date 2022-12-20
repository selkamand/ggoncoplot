test_that("ggoncoplot_plot_gene_barplot works", {
  gbm_csv <- system.file(
    package = "ggoncoplot",
    "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
  )

  gbm_df <- read.csv(file = gbm_csv, header = TRUE)


  # Throws no error
  expect_error(
    gg <- ggoncoplot(
    .data = gbm_df,
    col_genes = 'Hugo_Symbol',
    col_samples = 'Tumor_Sample_Barcode',
    col_mutation_type = 'Variant_Classification',
    interactive = TRUE,
    draw_gene_barplot = TRUE,
    verbose = FALSE
  ),
  NA
  )

  # Create oncoplot with gene_barplot
  ggiraph_gene_barplot <- ggoncoplot(
    .data = gbm_df,
    col_genes = 'Hugo_Symbol',
    col_samples = 'Tumor_Sample_Barcode',
    col_mutation_type = 'Variant_Classification',
    interactive = TRUE,
    draw_gene_barplot = TRUE,
    verbose = FALSE
  )


  # Test ggiraph object looks good (same as previous snapshot)
  # Otherwise will throw error
  expect_snapshot(ggiraph_gene_barplot)

})
