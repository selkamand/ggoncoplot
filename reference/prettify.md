# Make strings prettier for printing

Takes an input string and 'prettify' by converting underscores to
spaces, capitalising each word, etc.

## Usage

``` r
prettify(string, space_after_apostrophe = TRUE, autodetect_units = TRUE)
```

## Arguments

- string:

  input string

- space_after_apostrophe:

  add a space after any apostrophe so long as its after an alphanumeric
  character and followed by anything but a space (flag)

- autodetect_units:

  automatically detect units (e.g. mm, kg, etc) and wrap in brackets.

## Value

string
