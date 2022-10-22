test_that("check_valid_dataframe_column works", {

  # runs without error when columns are in data
  expect_error(
    check_valid_dataframe_column(mtcars, c('cyl', 'mpg')),
    NA
  )

  # runs without error when columns are in data
  expect_error(
    check_valid_dataframe_column(mtcars, 'mpg'),
    NA
  )


  # Throws expected error when columns are missing
  expect_error(
    check_valid_dataframe_column(mtcars, 'asdasd'),
    'asdasd'
  )

  # Throws expected error when columns are missing
  expect_error(
    check_valid_dataframe_column(mtcars, c('cyl','asdasd')),
    'asdasd'
  )

  # Throws expected error when columns are missing
  expect_error(
    check_valid_dataframe_column(mtcars, c('asdasd', 'cyl', 'bob')),
    'asdasd'
  )
})
