# ggoncoplot throws appropriate errors when sample or gene data is missing from mutation datasets

    Code
      ggoncoplot(df_mutations_na_samples, col_genes = "Genes", col_samples = "Samples")
    Condition
      Error in `ggoncoplot()`:
      ! 'data[[col_samples]]' must have no missing values! Found 1

---

    Code
      ggoncoplot(df_mutations_emptystring_samples, col_genes = "Genes", col_samples = "Samples")
    Condition
      Error in `ggoncoplot()`:
      ! Sample column cannot contain zero-length strings

---

    Code
      ggoncoplot(df_mutations_na_genes, col_genes = "Genes", col_samples = "Samples")
    Condition
      Error in `ggoncoplot()`:
      ! 'data[[col_genes]]' must have no missing values! Found 1

---

    Code
      ggoncoplot(df_mutations_emptystring_genes, col_genes = "Genes", col_samples = "Samples")
    Condition
      Error in `ggoncoplot()`:
      ! Gene column cannot contain zero-length strings

# ggoncoplot throws appropriate errors when clinical metadata has duplicate rows

    Code
      ggoncoplot(df_mutations, col_genes = "Genes", col_samples = "Samples",
        metadata = df_clinical_invalid)
    Condition
      Error in `ggoncoplot()`:
      ! 'Metadata Sample Column' must have no duplicates! Found 1 duplicated value: SampleB

# ggoncoplot throws error if metadata isn't a dataframe

    Code
      ggoncoplot(df_mutations, col_genes = "Genes", col_samples = "Samples",
        metadata = 1:10)
    Condition
      Error in `ggoncoplot()`:
      ! 'metadata' must be a data.frame, not a integer

# ggoncoplot throws eror if Mutation Type column is empty or NA

    Code
      ggoncoplot(df_mutations_valid_sample_genes, col_genes = "Genes", col_samples = "Samples",
        col_mutation_type = "VariantType", topn = Inf)
    Condition
      Error in `ggoncoplot()`:
      ! 'Mutation Type Column: VariantType' must have no missing values! Found 3

