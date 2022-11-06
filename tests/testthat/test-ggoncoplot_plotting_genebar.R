test_that("ggoncoplot_plot_gene_barplot works", {
  gbm_csv <- system.file(
    package = "ggoncoplot",
    "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
  )

  gbm_df <- read.csv(file = gbm_csv, header = TRUE)

  genes_for_oncoplot_top5 <- get_genes_for_oncoplot(
    .data = gbm_df,
    col_genes = "Hugo_Symbol",
    col_samples ="Tumor_Sample_Barcode", topn = 5,
    verbose = FALSE
  )

  df_top_genes <- ggoncoplot_prep_df(
    gbm_df,
    col_genes = "Hugo_Symbol",
    col_samples ="Tumor_Sample_Barcode",
    col_mutation_type = "Variant_Classification",
    genes_for_oncoplot = genes_for_oncoplot_top5
  )

  ggoncoplot_plot_gene_barplot(df_top_genes)
})
