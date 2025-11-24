# Scale probabilities by odds ratios

\$\$ O = \frac{p}{1-p} \\ p = \frac{O}{1+O} \\ p_k =
\frac{\frac{kp}{1-p}}{( 1 + \frac{kp}{1-p})} \\ p_k = \frac{kp}{1 - p +
kp} \$\$

## Usage

``` r
scale_probability(p, odds_ratio)
```

## Arguments

- p:

  A vector of probabilities

- odds_ratio:

  a vector of odds ratios

## Value

an adjusted vector of probabilities

## Unit tests


    scale_probability(0.5, c(0.25,0.5,1,2,4))
