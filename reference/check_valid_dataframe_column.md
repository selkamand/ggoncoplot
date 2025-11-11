# data.frame has colnames

Assert that data.frame contains a set of user defined column names.

## Usage

``` r
check_valid_dataframe_column(data, colnames, error_call = rlang::caller_env())
```

## Arguments

- data:

  dataframe that you want to assert contain specific columns
  (data.frame)

- colnames:

  Name (character)

- error_call:

  error call environment (do not change)

## Value

invisibly returns TRUE. If data is missing columns, will throw error

## Details

data.frame may have any additional colnames. It just has to have AT
LEAST the columns specified by `colnames`

Informs user about the missing columns one at a time. This may change in
future

## Examples

``` r
# Check mtcars has columns 'mpg' and 'cyl'
ggoncoplot:::check_valid_dataframe_column(mtcars, c("mpg", "cyl"))
```
