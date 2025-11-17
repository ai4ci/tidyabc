# logit scale

Perform logit scaling with right axis formatting. To not be used
directly but with ggplot (e.g. ggplot2::scale_y_continuous(trans =
"logit")

## Usage

``` r
.gg_transform_logit(n = 5, ...)
```

## Arguments

- n:

  number of breas

- ...:

  not used

## Value

A scales object

## Unit tests


    dplyr::tibble(pvalue = c(0.001, 0.05, 0.1), fold_change = 1:3) 
     ggplot2::ggplot(ggplot2::aes(fold_change , pvalue)) +
     ggplot2::geom_point() +
     ggplot2::scale_y_continuous(transform= .gg_transform_logit())
