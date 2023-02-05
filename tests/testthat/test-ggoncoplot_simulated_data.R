# Illegal versions (missing data / for clinical data: duplicated samples)
df_mutations_invalid <- read.csv(system.file(package = "ggoncoplot", "testdata/simulated_mutations.tsv"), sep = "\t", header = TRUE, check.names = FALSE)
df_clinical_invalid <- read.csv(system.file(package = "ggoncoplot", "testdata/simulated_mutations.clinical.tsv"), sep = "\t", header = TRUE, check.names = FALSE)

# Ensure no empty Sample / Gene Entries (filter them out)
df_mutations_valid_sample_genes <- subset(
  x = df_mutations_invalid,
  subset = !is.na(Samples) &
    Samples != "" &
    !is.na(Genes) &
    Genes != ""
)

# Ensure no empty VariantType entries  (mutate into Missense_Mutation)
df_mutations <- transform(
  `_data` = df_mutations_valid_sample_genes,
  VariantType = ifelse(
    is.na(VariantType) | nchar(VariantType, keepNA = TRUE) == 0, "Missense_Mutation", VariantType
    )
)

# Valid version of df_clinical
df_clinical <- subset(df_clinical_invalid, !duplicated(Samples))

test_that("ggoncoplot throws appropriate errors when sample or gene data is missing from mutation datasets", {

  # Create versions of df_mutations with only one type of illegal input (so order of assertions doesn't matter
  df_mutations_na_samples <- subset(df_mutations_invalid,
    subset =
    # !is.na(Samples) &
      nchar(Samples, keepNA = FALSE) > 0 &
      !is.na(Genes) &
      nchar(Genes, keepNA = FALSE) > 0
  )
  df_mutations_emptystring_samples <- subset(df_mutations_invalid,
    subset =
      !is.na(Samples) &
      !is.na(Genes) &
      nchar(Genes, keepNA = FALSE) > 0
  )

  df_mutations_na_genes <- subset(df_mutations_invalid,
    subset =
      !is.na(Samples) &
      nchar(Samples, keepNA = FALSE) > 0 &
      nchar(Genes, keepNA = FALSE) > 0
  )

  df_mutations_emptystring_genes <- subset(df_mutations_invalid,
    subset =
      !is.na(Samples) &
      nchar(Samples, keepNA = FALSE) > 0 &
      !is.na(Genes)
  )

  # Throw error due to sampleIDs = NA
  expect_snapshot(
    ggoncoplot(
      df_mutations_na_samples,
      col_genes = "Genes",
      col_samples = "Samples"
  ), error = TRUE)

  #Throw error due to sampleIDs = ""
  expect_snapshot(
    ggoncoplot(
      df_mutations_emptystring_samples,
      col_genes = "Genes",
      col_samples = "Samples"
    ), error = TRUE)


  # Throw error due to gene names = NA
  expect_snapshot(
    ggoncoplot(
      df_mutations_na_genes,
      col_genes = "Genes",
      col_samples = "Samples"
    ), error = TRUE)

  # Throw error due to gene names = ""
  expect_snapshot(
    ggoncoplot(
      df_mutations_emptystring_genes,
      col_genes = "Genes",
      col_samples = "Samples"
    ), error = TRUE)
})


test_that("ggoncoplot throws appropriate errors when clinical metadata has duplicate rows", {

  # Throw error due to duplicate sampleID (sample B)
  expect_snapshot(
    ggoncoplot(
      df_mutations,
      col_genes = "Genes",
      col_samples = "Samples",
      metadata = df_clinical_invalid
    ), error = TRUE)
})

test_that("ggoncoplot throws error if metadata isn't a dataframe", {

  # Throw error because metadata isn't a dataframe
  expect_snapshot(
    ggoncoplot(
      df_mutations,
      col_genes = "Genes",
      col_samples = "Samples",
      metadata = 1:10
    ), error = TRUE)

})



test_that("ggoncoplot works with valid data", {

  # Throw no errors with valid data
  expect_error(
    suppressMessages(ggoncoplot(
      df_mutations,
      col_genes = "Genes",
      col_samples = "Samples",
      metadata = df_clinical
    ), NA
  ))

})

test_that("ggoncoplot throws eror if Mutation Type column is empty or NA", {
  expect_snapshot(
    ggoncoplot(
      df_mutations_valid_sample_genes,
      col_genes = "Genes",
      col_samples = "Samples",
      col_mutation_type = "VariantType",
      topn = Inf
    ), error = TRUE
  )
})

test_that("ggoncoplot doesn't drop genes if all variant_type = NA", {

  #NFE2L2 & EMPTY should both be in final plot

  gg <- suppressMessages(ggoncoplot(
    df_mutations,
    col_genes = "Genes",
    col_samples = "Samples",
    col_mutation_type = "VariantType", interactive = FALSE, topn = Inf
  ))

  genes_plot <- ggplot2::layer_scales(gg)[['y']][["range"]][['range']]

  expect_true(all(c("NFE2L2", "EMPTY") %in% c(genes_plot)))
})




