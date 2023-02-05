# test_that("ggoncoplot runs without error", {
#   gbm_csv <- system.file(
#     package = "ggoncoplot",
#     "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
#   )
#
#   gbm_df <- read.csv(file = gbm_csv, header = TRUE)
#
#   # Runs without error (static)
#   expect_error(
#     gg <- ggoncoplot(
#       gbm_df,
#       col_genes = "Hugo_Symbol",
#       col_samples = "Tumor_Sample_Barcode",
#       col_mutation_type = "Variant_Classification",
#       interactive = FALSE
#     ) |> suppressMessages(),
#     regexp = NA
#   )
#
#   # Runs without error (interactive)
#   expect_error(
#     gg <- ggoncoplot(
#       gbm_df,
#       col_genes = "Hugo_Symbol",
#       col_samples = "Tumor_Sample_Barcode",
#       col_mutation_type = "Variant_Classification",
#       interactive = TRUE
#     ) |> suppressMessages(),
#     regexp = NA
#   )
#
#   # Produces interactive ggiraph obj
#   expect_s3_class(
#     gg <- ggoncoplot(
#       gbm_df,
#       col_genes = "Hugo_Symbol",
#       col_samples = "Tumor_Sample_Barcode",
#       col_mutation_type = "Variant_Classification",
#       interactive = TRUE
#     ) |> suppressMessages(),
#     class = "girafe"
#   )
#
#
#   # Produces interactive ggiraph obj
#   expect_s3_class(
#     gg <- ggoncoplot(
#       gbm_df,
#       col_genes = "Hugo_Symbol",
#       col_samples = "Tumor_Sample_Barcode",
#       col_mutation_type = "Variant_Classification",
#       interactive = TRUE
#     ) |> suppressMessages(),
#     class = "girafe"
#   )
#
#
#   # Produces static ggplot
#   expect_s3_class(
#     gg <- ggoncoplot(
#       gbm_df,
#       col_genes = "Hugo_Symbol",
#       col_samples = "Tumor_Sample_Barcode",
#       col_mutation_type = "Variant_Classification",
#       interactive = FALSE
#     ) |> suppressMessages(),
#     class = "ggplot"
#   )
#
#   # Runs without error when specifying custom tooltip
#   expect_error(
#     gg <- gbm_df |>
#       dplyr::mutate(tooltip = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) |>
#       ggoncoplot(
#         col_genes = 'Hugo_Symbol',
#         col_samples = 'Tumor_Sample_Barcode',
#         col_mutation_type = 'Variant_Classification',
#         col_tooltip = 'tooltip', draw_gene_barplot = FALSE # We'll specify a custom tooltip based on our new 'tooltip' column
#       ) |> suppressMessages(),
#     NA
#   )
#
# })
#
# # Create function for testing ggplot output
# get_ggplot_axis_text <- function(gg, axis = c('x', 'y')){
#   axis <- rlang::arg_match(axis)
#   ggplot2::ggplot_build(gg)$layout$panel_params[[1]][[axis]]$get_labels()
# }
#
#
# test_that("ggoncoplot axis text are appropriate", {
#   gbm_csv <- system.file(
#     package = "ggoncoplot",
#     "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
#   )
#
#   gbm_df <- read.csv(file = gbm_csv, header = TRUE)
#
#   ggtest <- ggoncoplot(
#     gbm_df,
#     col_genes = "Hugo_Symbol",
#     col_samples = "Tumor_Sample_Barcode",
#     col_mutation_type = "Variant_Classification",
#     topn = 10,
#     interactive = FALSE,
#     draw_gene_barplot = FALSE
#   ) |> suppressMessages()
#
#   # Test y axis hasn't changed
#   expect_snapshot(get_ggplot_axis_text(ggtest, 'y'))
#
#   # Test x axis hasn't changed
#   expect_snapshot(get_ggplot_axis_text(ggtest, 'x'))
# })
#
# test_that("ggoncoplot metadata works", {
#
#   # Paths
#   gbm_csv <- system.file(
#     package = "ggoncoplot",
#     "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
#   )
#
#   gbm_clinical_csv <- system.file(
#     package = "ggoncoplot",
#     "testdata/GBM_tcgamutations_mc3_clinical.csv"
#   )
#
#   gbm_clinical_duplicated_csv <- system.file(
#     package = "ggoncoplot",
#     "testdata/GBM_tcgamutations_mc3_clinical.duplicates.csv"
#   )
#
#   # Read Data
#   df_gbm <- read.csv(file = gbm_csv, header = TRUE)
#   df_gbm_clinical <- read.csv(file = gbm_clinical_csv, header = TRUE)
#   df_gbm_clinical_duplicates <- read.csv(file = gbm_clinical_duplicated_csv, header = TRUE)
#
#
#   # Expect no error
#   expect_error(
#     gg <- ggoncoplot(
#       df_gbm,
#       col_samples = "Tumor_Sample_Barcode",
#       col_genes = "Hugo_Symbol",
#       col_mutation_type = "Variant_Classification",
#       metadata = df_gbm_clinical,
#       cols_to_plot_metadata = c('gender'),
#       verbose = FALSE
#     ),
#     regexp = NA
#   )
#
#   fig <- ggoncoplot(
#     df_gbm,
#     col_samples = "Tumor_Sample_Barcode",
#     col_genes = "Hugo_Symbol",
#     col_mutation_type = "Variant_Classification",
#     metadata = df_gbm_clinical,
#     cols_to_plot_metadata = c('gender'),
#     verbose = FALSE
#   )
#
#   # Expect general doppelganger
#   vdiffr::expect_doppelganger(title = "GBM clinical metadata", fig = fig)
#
#   # Expect error due to duplicates
#   expect_snapshot(
#     gg <- ggoncoplot(
#       df_gbm,
#       col_samples = "Tumor_Sample_Barcode",
#       col_genes = "Hugo_Symbol",
#       col_mutation_type = "Variant_Classification",
#       metadata = df_gbm_clinical_duplicates,
#       cols_to_plot_metadata = c('gender')
#     ),
#     error = TRUE
#   )
#
# })
#
