test_that("score_based_on_gene_rank works", {

  # Code runs without error
  expect_error(
    score_based_on_gene_rank(c("TERT", "EGFR", "PTEN", "BRCA2"), c("EGFR", "BRCA2"), gene_rank = 1:2),
    regexp = NA
    )


   # First set of genes has a high rank since both BRCA2 and EGFR are mutated
  expect_equal(
    score_based_on_gene_rank(c("TERT", "EGFR", "PTEN", "BRCA2"), c("EGFR", "BRCA2"), gene_rank = 1:2),
    3
    )

  # If EGFR is mutated without BRCA2, we get a lower score
  expect_equal(
    score_based_on_gene_rank(
     mutated_genes = c("TERT", "EGFR", "PTEN", "IDH1"),
     genes_informing_score = c("EGFR", "BRCA2"),
     gene_rank = 1:2
     ),
    expected = 1
  )

  # If BRCA2 is mutated without EGFR,
  # we get a score lower than BRCA2+EGFR but higher than EGFR alone due to higher gene_rank of BRCA2
  expect_equal(
    score_based_on_gene_rank(
      mutated_genes = c("TERT", "IDH1", "PTEN", "BRCA2"),
      genes_informing_score = c("EGFR", "BRCA2"),
      gene_rank = 1:2
    ),
    expected = 2
  )

})
