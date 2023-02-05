# test_that('ggoncoplot_prep_df works', {
#
#
#   # Load data ---------------------------------------------------------------
#   gbm_csv <- system.file(
#     package = "ggoncoplot",
#     "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
#   )
#
#   gbm_df <- read.csv(file = gbm_csv, header = TRUE)
#
#   # Create some genelists for later
#   genes_for_oncoplot_top5 <- get_genes_for_oncoplot(
#     .data = gbm_df,
#     col_genes = "Hugo_Symbol",
#     col_samples ="Tumor_Sample_Barcode", topn = 5,
#     verbose = FALSE
#   )
#
#   genes_for_oncoplot_top10 <- get_genes_for_oncoplot(
#     .data = gbm_df,
#     col_genes = "Hugo_Symbol",
#     col_samples ="Tumor_Sample_Barcode", topn = 10,
#     verbose = FALSE
#   )
#
#   genes_for_oncoplot_top10_with_ties <- get_genes_for_oncoplot(
#     .data = gbm_df,
#     col_genes = "Hugo_Symbol",
#     col_samples ="Tumor_Sample_Barcode", topn = 10,
#     return_extra_genes_if_tied = TRUE,
#     verbose = FALSE
#   )
#
#   genes_for_oncoplot_top10_ignore_ttn_nf1 <- get_genes_for_oncoplot(
#     .data = gbm_df,
#     col_genes = "Hugo_Symbol",
#     col_samples ="Tumor_Sample_Barcode", topn = 10,
#     genes_to_ignore = c("TTN", "NF1"),
#     verbose = FALSE
#   )
#
#   genes_to_include <- c(
#     "ATRX", "PIK3R1", "SPTA1",
#     "RYR2", "NF1", "FLG",
#     "MUC16", "EGFR", "TTN",
#     "TP53", "PTEN"
#   )
#
#   genes_for_oncoplot_custom_genelist <- get_genes_for_oncoplot(
#     .data = gbm_df,
#     col_genes = "Hugo_Symbol",
#     col_samples ="Tumor_Sample_Barcode", topn = 10,
#     genes_to_include = genes_to_include,
#     verbose = FALSE
#   )
#
# # Actual Testing ----------------------------------------------------------
#   # Runs without error
#   expect_error(
#     ggoncoplot_prep_df(
#       gbm_df,
#       col_genes = "Hugo_Symbol",
#       col_samples ="Tumor_Sample_Barcode",
#       col_mutation_type = "Variant_Classification",
#       genes_for_oncoplot = genes_for_oncoplot_top5
#     ),
#     NA
#   )
#
#   # Create dataframe to use as input for downstream tests
#   prepped_df <- ggoncoplot_prep_df(
#     gbm_df,
#     col_genes = "Hugo_Symbol",
#     col_samples ="Tumor_Sample_Barcode",
#     col_mutation_type = "Variant_Classification",
#     genes_for_oncoplot = genes_for_oncoplot_top5
#   )
#
#
#
#   # same dataframe without col_mutation_type specified
#   prepped_df_no_mutation_type <- ggoncoplot_prep_df(
#     gbm_df,
#     col_genes = "Hugo_Symbol",
#     col_samples ="Tumor_Sample_Barcode",
#     genes_for_oncoplot = genes_for_oncoplot_top5
#   )
#
#
#
#   # Testing Default  ----------------------------------------------------------
#
#   ## Outputs Correct ----------------------------------------------------------
#   # Check result is dataframe
#   expect_s3_class(
#     prepped_df, 'data.frame'
#   )
#   expect_s3_class(
#     prepped_df_no_mutation_type, 'data.frame'
#   )
#
#   # Check dataframe has required names
#   expect_named(prepped_df, expected = c('Sample', 'Gene', 'MutationType', 'MutationCount', 'Tooltip'), ignore.order = TRUE)
#   expect_named(prepped_df_no_mutation_type, expected = c('Sample', 'Gene', 'MutationType', 'MutationCount', 'Tooltip'), ignore.order = TRUE)
#
#   # Check that there is one row per sample-gene, never two
#   # if one sample has multiple mutations in a gene we'd expect them to be collapse into a single row with MutationType == 'multiple')
#   expect_true(nrow(dplyr::distinct(prepped_df, Sample, Gene)) == nrow(prepped_df))
#   expect_true(nrow(dplyr::distinct(prepped_df_no_mutation_type, Sample, Gene)) == nrow(prepped_df_no_mutation_type))
#
#
#   # Check returned dataframe is not grouped
#   expect_false(dplyr::is.grouped_df(prepped_df))
#   expect_false(dplyr::is.grouped_df(prepped_df_no_mutation_type))
#
#   # Sample should be a factor
#   expect_s3_class(prepped_df[['Sample']], "factor")
#   expect_s3_class(prepped_df_no_mutation_type[['Sample']], "factor")
#
#   # Expect sample levels be ordered appropriately for oncoplot
#   expect_snapshot(levels(prepped_df[['Sample']]))
#   expect_snapshot(levels(prepped_df_no_mutation_type[['Sample']]))
#
#
#   # Mutation count as expected
#   expect_snapshot(prepped_df[["MutationCount"]])
#   expect_snapshot(prepped_df_no_mutation_type[["MutationCount"]])
#
#
#   # Expect full output of ggoncoplot_prep_df to be as expected
#   # We do the individual tests first since output is likely to be more informative but this
#   # full test gives us absolute confidence we haven't broken anything
#   expect_snapshot(prepped_df)
#   expect_snapshot(prepped_df_no_mutation_type)
#
#   # Specifying Custom tooltip leads to no errors
#   expect_error(
#     gbm_df |>
#       dplyr::mutate(tooltip = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) |>
#       ggoncoplot_prep_df(
#         col_genes = 'Hugo_Symbol',
#         col_samples = 'Tumor_Sample_Barcode',
#         col_mutation_type = 'Variant_Classification',
#         col_tooltip = 'tooltip', # We'll specify a custom tooltip based on our new 'tooltip' column
#         genes_for_oncoplot = genes_for_oncoplot_top10
#       ),
#     NA
#   )
#
#
#   # Testing Gene Selection Works  ----------------------------------------------
#
#   ## Top N ------------------------------------------------------------------
#
#   #Produces expected result when getting top10 genes return_extra_genes_if_tied = FALSE
#   expect_equal(
#     gbm_df |>
#       dplyr::mutate(tooltip = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) |>
#       ggoncoplot_prep_df(
#         col_genes = 'Hugo_Symbol',
#         col_samples = 'Tumor_Sample_Barcode',
#         col_mutation_type = 'Variant_Classification',
#         genes_for_oncoplot = genes_for_oncoplot_top10,
#         verbose = FALSE
#       ) |>
#       dplyr::pull(Gene) |>
#       levels(),
#     c("PTEN", "TP53", "TTN", "EGFR", "MUC16", "FLG", "NF1", "RYR2",
#       "ATRX", "PIK3R1")
#   )
#
#
#   # Works even if return_extra_genes_if_tied = TRUE
#   expect_error(
#     gbm_df |>
#       dplyr::mutate(tooltip = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) |>
#       ggoncoplot_prep_df(
#         col_genes = 'Hugo_Symbol',
#         col_samples = 'Tumor_Sample_Barcode',
#         col_mutation_type = 'Variant_Classification',
#         genes_for_oncoplot = genes_for_oncoplot_top10_with_ties,
#         verbose = FALSE
#       ),
#     NA
#   )
#
#
#   # Returns expected genes when `return_extra_genes_if_tied == TRUE`
#   expect_equal(
#     gbm_df |>
#       dplyr::mutate(tooltip = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) |>
#       ggoncoplot_prep_df(
#         col_genes = 'Hugo_Symbol',
#         col_samples = 'Tumor_Sample_Barcode',
#         col_mutation_type = 'Variant_Classification',
#         genes_for_oncoplot = genes_for_oncoplot_top10_with_ties,
#         verbose = FALSE
#       ) |>
#       dplyr::pull(Gene) |>
#       levels(),
#     c("PTEN", "TP53", "TTN", "EGFR", "MUC16", "FLG", "NF1", "RYR2",
#       "ATRX", "PIK3R1", "SPTA1")
#   )
#
#   # Ensure verbose working
#   # expect_message(
#   #   gbm_df |>
#   #     dplyr::mutate(tooltip = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) |>
#   #     ggoncoplot_prep_df(
#   #       col_genes = 'Hugo_Symbol',
#   #       col_samples = 'Tumor_Sample_Barcode',
#   #       col_mutation_type = 'Variant_Classification',
#   #       genes_for_oncoplot = genes_for_oncoplot_top10_with_ties,
#   #       verbose = TRUE
#   #     ),
#   #   "due to ties, the top 11 genes were returned"
#   # )
#
#   # Ignoring genes works
#   expect_equal(
#     gbm_df |>
#       dplyr::mutate(tooltip = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) |>
#       ggoncoplot_prep_df(
#         col_genes = 'Hugo_Symbol',
#         col_samples = 'Tumor_Sample_Barcode',
#         col_mutation_type = 'Variant_Classification',
#         genes_for_oncoplot = genes_for_oncoplot_top10_ignore_ttn_nf1,
#         verbose = FALSE
#       ) |>
#       dplyr::pull(Gene) |>
#       levels(),
#     c("PTEN", "TP53", "EGFR", "MUC16", "FLG", "RYR2",
#       "ATRX", "PIK3R1", "SPTA1", "PIK3CA")
#   )
#
#
# ## Custom Geneset ----------------------------------------------------------
#   # Diltering for specific genes properly filters for these genes
#   # First define 11 genes to filter on
#   genes_to_include <- c(
#     "ATRX", "PIK3R1", "SPTA1",
#     "RYR2", "NF1", "FLG",
#     "MUC16", "EGFR", "TTN",
#     "TP53", "PTEN"
#   )
#   #
#   # Expect no error when requesting genes
#   expect_error(
#     ggoncoplot_prep_df(
#       gbm_df,
#       col_genes = "Hugo_Symbol",
#       col_samples ="Tumor_Sample_Barcode",
#       col_mutation_type = "Variant_Classification",
#       genes_for_oncoplot = genes_for_oncoplot_custom_genelist
#     ),
#     NA
#   )
#   #
#   # # non mutation type version
#   expect_error(
#     ggoncoplot_prep_df(
#       gbm_df,
#       col_genes = "Hugo_Symbol",
#       col_samples ="Tumor_Sample_Barcode",
#       #col_mutation_type = "Variant_Classification",
#       genes_for_oncoplot = genes_for_oncoplot_custom_genelist
#     ),
#     NA
#   )
#
#
#   # # Create df for future tests
#   prepped_df_custom_genes <- ggoncoplot_prep_df(
#     gbm_df,
#     col_genes = "Hugo_Symbol",
#     col_samples ="Tumor_Sample_Barcode",
#     col_mutation_type = "Variant_Classification",
#     genes_for_oncoplot = genes_for_oncoplot_custom_genelist
#   )
#
#
#   # non mutation type version
#   prepped_df_custom_genes_no_mutation_type <- ggoncoplot_prep_df(
#     gbm_df,
#     col_genes = "Hugo_Symbol",
#     col_samples ="Tumor_Sample_Barcode",
#     #col_mutation_type = "Variant_Classification",
#     genes_for_oncoplot = genes_for_oncoplot_custom_genelist
#   )
#
#
#   # Expect included genes to be level'd appropriately
#   expect_equal(
#     levels(prepped_df_custom_genes[["Gene"]]),
#     genes_to_include
#   )
#
#
# })
