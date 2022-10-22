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

  # Verbose test?
})
