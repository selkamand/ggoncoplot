test_that("identify_topn_genes works", {
  gbm_csv <- system.file(
    package = "ggoncoplot",
    "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
  )

  gbm_df <- read.csv(file = gbm_csv, header = TRUE)

  # Throws no errors
  expect_error(
    identify_topn_genes(
      .data = gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      topn = 10,
      return_extra_genes_if_tied = FALSE
      ),
    regexp = NA
  )

  # Throws error on no topn specified error
  expect_error(
    identify_topn_genes(
      .data = gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      return_extra_genes_if_tied = FALSE
    ),
    regexp = 'argument "topn" is missing, with no default'
  )

  # Produces expected result
  expect_equal(
    identify_topn_genes(
      .data = gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      topn = 10,
      return_extra_genes_if_tied = FALSE
    ),
    c("PTEN", "TP53", "TTN", "EGFR", "MUC16", "FLG", "NF1", "RYR2",
      "ATRX", "PIK3R1")
  )

  # Produces expected result when return_extra_genes_if_tied = TRUE
  expect_equal(
    identify_topn_genes(
      .data = gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      topn = 10,
      return_extra_genes_if_tied = TRUE
    ) |> suppressMessages(),
    c("PTEN", "TP53", "TTN", "EGFR", "MUC16", "FLG", "NF1", "RYR2",
      "ATRX", "PIK3R1", "SPTA1")
  )

  # Produces expected warning when return_extra_genes_if_tied = TRUE
  expect_message(
    identify_topn_genes(
      .data = gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      topn = 10,
      return_extra_genes_if_tied = TRUE
    ),
    "due to ties, the top 11 genes were returned"
  )

  # Warning message not produced if verbose = FALSE  return_extra_genes_if_tied = TRUE
  expect_message(
    identify_topn_genes(
      .data = gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      topn = 10,
      return_extra_genes_if_tied = TRUE,
      verbose = FALSE
    ),
    NA
  )

 # Ignore genes works
  expect_equal(
    identify_topn_genes(
      .data = gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples = "Tumor_Sample_Barcode",
      topn = 10,
      return_extra_genes_if_tied = FALSE,
      genes_to_ignore = c("TTN", "NF1")
    ),
    c("PTEN", "TP53", "EGFR", "MUC16", "FLG", "RYR2",
      "ATRX", "PIK3R1", "SPTA1", "PIK3CA")
  )

  # If you only have 20 genes in input data but request top 100 be displayed, should return all 20 (and throw no error)
  expect_length(
    gbm_df |>
      dplyr::distinct(Hugo_Symbol, .keep_all = TRUE) |>
      head(n=20) |>
      identify_topn_genes(
        col_genes = "Hugo_Symbol",
        col_samples = "Tumor_Sample_Barcode",
        topn = 100,
        return_extra_genes_if_tied = FALSE,
        verbose = FALSE
      ),
    20
    )

  # Proce
  expect_message(
    gbm_df |>
      dplyr::distinct(Hugo_Symbol, .keep_all = TRUE) |>
      head(n=20) |>
      identify_topn_genes(
        col_genes = "Hugo_Symbol",
        col_samples = "Tumor_Sample_Barcode",
        topn = 100,
        return_extra_genes_if_tied = FALSE,
        verbose = TRUE
      ),
    "not enough genes are in data"
  )


})

# This function wraps the above, offering additional functionality to handle user-specified genelists
test_that("get_genes_for_oncoplot works", {
  gbm_csv <- system.file(
    package = "ggoncoplot",
    "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
  )

  gbm_df <- read.csv(file = gbm_csv, header = TRUE)


  ## Topn --------------------------------------------------------------------
  # Runs without error
  expect_error(
    get_genes_for_oncoplot(
      .data = gbm_df,
      col_samples = "Tumor_Sample_Barcode",
      col_genes = "Hugo_Symbol",
      topn = 10
    ),
    NA
  )

  # Ensure verbose working
  expect_message(
      get_genes_for_oncoplot(
        .data = gbm_df,
        col_samples = "Tumor_Sample_Barcode",
        col_genes = "Hugo_Symbol",
        topn = 10,
        return_extra_genes_if_tied = TRUE,
        verbose = TRUE
      ),
    "due to ties, the top 11 genes were returned"
  )

  # Ignoring genes works
  expect_equal(
    get_genes_for_oncoplot(
      .data = gbm_df,
      col_samples = "Tumor_Sample_Barcode",
      col_genes = "Hugo_Symbol",
      topn = 10,
      genes_to_ignore = c("TTN", "NF1"),
      verbose = TRUE
    ),
    c("PTEN", "TP53", "EGFR", "MUC16", "FLG", "RYR2",
      "ATRX", "PIK3R1", "SPTA1", "PIK3CA")
  )



  ## Custom Geneset ----------------------------------------------------------
  # Filtering for specific genes properly filters for these genes
  # First define 11 genes to filter on
  genes_to_include <- c(
    "ATRX", "PIK3R1", "SPTA1",
    "RYR2", "NF1", "FLG",
    "MUC16", "EGFR", "TTN",
    "TP53", "PTEN"
  )


  # # Expect no error when requesting genes
  expect_error(
    get_genes_for_oncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples ="Tumor_Sample_Barcode",
      genes_to_include = genes_to_include
    ),
    NA
  )




  # Expected values returned when specifying genes
  expect_equal(
    get_genes_for_oncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples ="Tumor_Sample_Barcode",
      genes_to_include = genes_to_include
    ),
    genes_to_include
  )

  # Works when excluding genes + specifying geneset to include
  genes_to_exclude = c("ATRX", "PIK3R1", "RandomGene")
  expect_equal(
    get_genes_for_oncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples ="Tumor_Sample_Barcode",
      genes_to_ignore = genes_to_exclude,
      genes_to_include = genes_to_include
    ),
    genes_to_include[!genes_to_include %in% genes_to_exclude]
  )

 # Return zero length char vector if genes to ignore covers all genes to include
  expect_equal(
    get_genes_for_oncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples ="Tumor_Sample_Barcode",
      genes_to_ignore = genes_to_include,
      genes_to_include = genes_to_include
    ),
    character(0)
  )

 # Throws error if no genes in genes_to_include are in dataset
  expect_error(
    get_genes_for_oncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples ="Tumor_Sample_Barcode",
      genes_to_include = c("asda"),
      verbose = FALSE
    ),
    "Couldn't find any of the genes you supplied in your dataset"
  )

  # Warn user if only some of the genes in genes_to_include are in the dataset, and just return those
  expect_warning(
    get_genes_for_oncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples ="Tumor_Sample_Barcode",
      genes_to_include = c("asda", "PTEN"),
      verbose = TRUE
    )
  )

  expect_equal(
    get_genes_for_oncoplot(
      gbm_df,
      col_genes = "Hugo_Symbol",
      col_samples ="Tumor_Sample_Barcode",
      genes_to_include = c("asda", "PTEN"),
      verbose = FALSE
    ),
    "PTEN"
  )
})

