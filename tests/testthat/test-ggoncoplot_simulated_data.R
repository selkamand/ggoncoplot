# Illegal versions (missing data / for clinical data: duplicated samples)
df_mutations_invalid <- read.csv(system.file(package = "ggoncoplot", "testdata/simulated_mutations.tsv"), sep = "\t", header = TRUE, check.names = FALSE)
df_clinical_invalid <- read.csv(system.file(package = "ggoncoplot", "testdata/simulated_mutations.clinical.tsv"), sep = "\t", header = TRUE, check.names = FALSE)

df_mutations <- subset(
  x = df_mutations_invalid,
  subset = !is.na(Samples) &
    Samples != "" &
    !is.na(Genes) &
    Genes != ""
)

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


  #vdiffr::expect_doppelganger(title = "Basic Simulated Oncoplot")
})
