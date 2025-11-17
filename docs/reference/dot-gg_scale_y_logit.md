# A logit y scale

A logit y scale

## Usage

``` r
.gg_scale_y_logit(..., n = 5, sf = 2)
```

## Arguments

- ...:

  Other arguments passed on to `scale_(x|y)_continuous()`

- n:

  the number of major breaks

- sf:

  significant figures

## Value

a ggplot scale

## Unit tests


    dplyr::tibble(pvalue = c(0.001, 0.05, 0.1), fold_change = 1:3) 
     ggplot2::ggplot(ggplot2::aes(fold_change , pvalue)) +
     ggplot2::geom_point() +
     .gg_scale_y_logit(n=8)
