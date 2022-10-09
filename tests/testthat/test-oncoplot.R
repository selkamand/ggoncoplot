test_that("ggoncoplot runs without error", {
  gbm_csv <- system.file(
    package = "ggoncoplot",
    "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
  )

  gbm_df <- read.csv(file = gbm_csv, header = TRUE)

  # Runs without error (static)
  expect_error(
    ggoncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      col_mutation_type = "Variant_Classification",
      interactive = FALSE
    ),
    regexp = NA
  )

  # Runs without error (interactive)
  expect_error(
    ggoncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      col_mutation_type = "Variant_Classification",
      interactive = TRUE
    ),
    regexp = NA
  )

  # Produces interactive ggiraph obj
  expect_s3_class(
    ggoncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      col_mutation_type = "Variant_Classification",
      interactive = TRUE
    ),
    class = "girafe"
  )


  # Produces interactive ggiraph obj
  expect_s3_class(
    ggoncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      col_mutation_type = "Variant_Classification",
      interactive = TRUE
    ),
    class = "girafe"
  )


  # Produces static ggplot
  expect_s3_class(
    ggoncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      col_mutation_type = "Variant_Classification",
      interactive = FALSE
    ),
    class = "ggplot"
  )

})

test_that('ggoncoplot_prep_df works', {

  gbm_csv <- system.file(
    package = "ggoncoplot",
    "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
  )

  gbm_df <- read.csv(file = gbm_csv, header = TRUE)

  # Runs without error
  expect_error(
    ggoncoplot_prep_df(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples ="Tumor_Sample_Barcode",
      col_mutation_type = "Variant_Classification"
    ),
    NA
  )

  # Create dataframe to use as input for downstream tests
  prepped_df <- ggoncoplot_prep_df(
    gbm_df,
    col_genes = "Hugo_Symbol",
    col_samples ="Tumor_Sample_Barcode",
    col_mutation_type = "Variant_Classification",
    topn = 10
  )

  # same dataframe without col_mutation_type specified
  prepped_df_no_mutation_type <- ggoncoplot_prep_df(
    gbm_df,
    col_genes = "Hugo_Symbol",
    col_samples ="Tumor_Sample_Barcode",
    topn = 10
  )


  # Check result is dataframe
  expect_s3_class(
    prepped_df, 'data.frame'
  )
  expect_s3_class(
    prepped_df_no_mutation_type, 'data.frame'
  )

  # Check dataframe has required names
  expect_named(prepped_df, expected = c('Sample', 'Gene', 'MutationType', 'Tooltip'), ignore.order = TRUE)
  expect_named(prepped_df_no_mutation_type, expected = c('Sample', 'Gene', 'MutationType', 'Tooltip'), ignore.order = TRUE)

  # Check that there is one row per sample-gene, never two
  # if one sample has multiple mutations in a gene we'd expect them to be collapse into a single row with MutationType == 'multiple')
  expect_true(nrow(dplyr::distinct(prepped_df, Sample, Gene)) == nrow(prepped_df))
  expect_true(nrow(dplyr::distinct(prepped_df_no_mutation_type, Sample, Gene)) == nrow(prepped_df_no_mutation_type))


  # Check returned dataframe is not grouped
  expect_false(dplyr::is.grouped_df(prepped_df))
  expect_false(dplyr::is.grouped_df(prepped_df_no_mutation_type))


  # Diltering for specific genes properly filters for these genes
  # First define 11 genes to filter on
  genes_to_include <- c(
    "ATRX", "PIK3R1", "SPTA1",
    "RYR2", "NF1", "FLG",
    "MUC16", "EGFR", "TTN",
    "TP53", "PTEN"
  )

  # Expect no error when requesting genes
  expect_error(
    ggoncoplot_prep_df(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples ="Tumor_Sample_Barcode",
      col_mutation_type = "Variant_Classification",
      genes_to_include = genes_to_include
    ),
    NA
  )

  # non mutation type version
  expect_error(
    ggoncoplot_prep_df(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples ="Tumor_Sample_Barcode",
      #col_mutation_type = "Variant_Classification",
      genes_to_include = genes_to_include
    ),
    NA
  )

  # Create df for future tests
  prepped_df_custom_genes <- ggoncoplot_prep_df(
    gbm_df,
    col_genes = "Hugo_Symbol",
    col_samples ="Tumor_Sample_Barcode",
    col_mutation_type = "Variant_Classification",
    genes_to_include = genes_to_include
  )

  # non mutation type version
  prepped_df_custom_genes_no_mutation_type <- ggoncoplot_prep_df(
    gbm_df,
    col_genes = "Hugo_Symbol",
    col_samples ="Tumor_Sample_Barcode",
    #col_mutation_type = "Variant_Classification",
    genes_to_include = genes_to_include
  )



  # Expect included genes to be those specifically requested
  expect_identical(
    sort(as.character(unique(prepped_df_custom_genes[['Gene']]))),
    sort(genes_to_include)
  )

  expect_identical(
    sort(as.character(unique(prepped_df_custom_genes_no_mutation_type[['Gene']]))),
    sort(genes_to_include)
  )

  # Specifying Custom tooltip leads to no errors
  expect_error(
  gbm_df |>
    mutate(tooltip = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) |>
    ggoncoplot(
      col_genes = 'Hugo_Symbol',
      col_samples = 'Tumor_Sample_Barcode',
      col_mutation_type = 'Variant_Classification',
      col_tooltip = 'tooltip' # We'll specify a custom tooltip based on our new 'tooltip' column
    ),
  NA
  )
})
