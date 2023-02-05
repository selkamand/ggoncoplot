# ggoncoplot throws appropriate errors when sample or gene data is missing from mutation datasets

    Code
      ggoncoplot(df_mutations_na_samples, col_genes = "Genes", col_samples = "Samples")
    Error <rlang_error>
      '.data[[col_samples]]' must have no missing values! Found 1

---

    Code
      ggoncoplot(df_mutations_emptystring_samples, col_genes = "Genes", col_samples = "Samples")
    Error <rlang_error>
      Sample column cannot contain zero-length strings

---

    Code
      ggoncoplot(df_mutations_na_genes, col_genes = "Genes", col_samples = "Samples")
    Error <rlang_error>
      '.data[[col_genes]]' must have no missing values! Found 1

---

    Code
      ggoncoplot(df_mutations_emptystring_genes, col_genes = "Genes", col_samples = "Samples")
    Error <rlang_error>
      Gene column cannot contain zero-length strings

