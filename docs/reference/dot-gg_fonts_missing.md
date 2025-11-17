# Check if any fonts listed are missing

Check if any fonts listed are missing

## Usage

``` r
.gg_fonts_missing(family)
```

## Arguments

- family:

  the font family

## Value

`TRUE` if missing fonts detected

## Unit tests


    .gg_fonts_missing("Arial")
    .gg_fonts_missing(c("Roboto","Kings","ASDASDAS"))
