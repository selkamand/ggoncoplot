test_that("ggoncoplot_options works", {

  # No Error
  expect_error(ggoncoplot_options(), NA)
  expect_s3_class(ggoncoplot_options(), "ggoncoplot_options")
})
