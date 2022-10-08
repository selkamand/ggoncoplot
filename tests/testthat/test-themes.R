test_that('oncoplot themes work', {

  # Runs without error
  expect_error(
    theme_oncoplot_default(),
    NA
  )

  # returns theme object
  expect_s3_class(
    theme_oncoplot_default(),
    "theme"
  )
})
