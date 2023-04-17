# Illegal versions (missing data / for clinical data: duplicated samples)
df_mutations_invalid <- read.csv(system.file(package = "ggoncoplot", "testdata/simulated_mutations.tsv"), sep = "\t", header = TRUE, check.names = FALSE)
df_clinical_invalid <- read.csv(system.file(package = "ggoncoplot", "testdata/simulated_mutations.clinical.tsv"), sep = "\t", header = TRUE, check.names = FALSE)
df_pathway_invalid <- read.csv(system.file(package = "ggoncoplot", "testdata/simulated_mutations.pathways.tsv"), sep = "\t", header = TRUE, check.names = FALSE)


# Create Valid Mutations Data --------------------------------------------------
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


# Create Valid Clinical --------------------------------------------------
# Valid version of df_clinical
df_clinical <- subset(df_clinical_invalid, !duplicated(Samples))

# Create Valid Pathway Data --------------------------------------------------
df_pathway <- df_pathway_invalid
# Already valid so do nothing


# Write output ------------------------------------------------------------
readr::write_tsv(df_mutations, here::here("inst/testdata/simulated_mutations.valid.tsv"))
readr::write_tsv(df_clinical, here::here("inst/testdata/simulated_mutations.clinical.valid.tsv"))
readr::write_tsv(df_pathway, here::here("inst/testdata/simulated_mutations.pathways.valid.tsv"))


# Test oncoplot
devtools::load_all()

df_mutations_valid <- read.csv(system.file(package = "ggoncoplot", "testdata/simulated_mutations.valid.tsv"), sep = "\t", header = TRUE, check.names = FALSE)
df_clinical_valid <- read.csv(system.file(package = "ggoncoplot", "testdata/simulated_mutations.clinical.valid.tsv"), sep = "\t", header = TRUE, check.names = FALSE)
df_pathway_valid <- read.csv(system.file(package = "ggoncoplot", "testdata/simulated_mutations.pathways.valid.tsv"), sep = "\t", header = TRUE, check.names = FALSE)

ggoncoplot(
  .data = df_mutations_valid, col_genes = "Genes", col_samples = "Samples", col_mutation_type = "VariantType",
  metadata = df_clinical_valid, pathway = df_pathway_valid,
  plotsize_metadata_rel_height = 40, interactive_svg_height = 10
  )
