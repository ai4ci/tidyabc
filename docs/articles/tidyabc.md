# Getting started

``` r
library(tidyabc)
#> 
#> Attaching package: 'tidyabc'
#> The following objects are masked from 'package:base':
#> 
#>     transform, truncate
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
```

``` r

# example simulation
# Well be trying to recover norm and gamma parameters
# We'll use this function for both example generation and fitting
test_simulation_fn = function(norm_mean, norm_sd, gamma_mean, gamma_sd) {
  
  A = rnorm(1000, norm_mean, norm_sd)
  B = rgamma2(1000, gamma_mean, gamma_sd)
  C = rgamma2(1000, gamma_mean, gamma_sd)

  return(
    list(
      data1 = A + B - C,
      data2 = B * C
    )
  )
  
}
```

``` r

test_scorer_fn = function(simdata, obsdata) {
  return(list(
    data1 = calculate_wasserstein(obsdata$data1, simdata$data1),
    data2 = calculate_wasserstein(obsdata$data2, simdata$data2)
  ))
}
```

``` r

tmp = test_simulation(
  test_simulation_fn,
  test_scorer_fn,
  norm_mean = 4, norm_sd=2, gamma_mean=6, gamma_sd=2,
  seed = 123
)

truth = tmp$truth
test_obsdata = tmp$obsdata
```

``` r

ggplot(
  tibble(data1 = test_obsdata$data1), aes(x=data1))+
  geom_histogram(,binwidth = 1)+
  xlab("A + B - C")
```

![](tidyabc_files/figure-html/unnamed-chunk-5-1.png)

``` r

ggplot(tibble(data2 = test_obsdata$data2), aes(x=data2))+
  geom_histogram(,binwidth = 1)+
  xlab("B x C")
```

![](tidyabc_files/figure-html/unnamed-chunk-5-2.png)

``` r
test_priors = priors(
  norm_mean ~ unif(0, 10),
  norm_sd ~ unif(0, 10),
  gamma_mean ~ unif(0, 10),
  gamma_sd ~ unif(0, 10),
  # enforces convex gamma:
  ~ gamma_mean > gamma_sd 
)

test_priors
#> Parameters: 
#> * norm_mean: unif(min = 0, max = 10)
#> * norm_sd: unif(min = 0, max = 10)
#> * gamma_mean: unif(min = 0, max = 10)
#> * gamma_sd: unif(min = 0, max = 10)
#> Constraints:
#> * gamma_mean > gamma_sd
```

``` r

rejection_fit = abc_rejection(
  obsdata = test_obsdata,
  priors_list = test_priors,
  sim_fn = test_simulation_fn,
  scorer_fn = test_scorer_fn,
  n_sims = 10000,
  acceptance_rate = 0.01,
  parallel= FALSE
)
#> ABC rejection, 1 wave.

summary(rejection_fit)
#> ABC rejection fit: single wave
#> Parameter estimates:
#> # A tibble: 4 × 4
#> # Groups:   param [4]
#>   param      mean_sd       median_95_CrI           ESS
#>   <chr>      <chr>         <chr>                 <dbl>
#> 1 gamma_mean 5.984 ± 0.319 6.004 [5.164 — 6.866]  70.9
#> 2 gamma_sd   2.426 ± 0.525 2.449 [1.314 — 3.617]  70.9
#> 3 norm_mean  4.046 ± 0.786 4.087 [2.122 — 5.856]  70.9
#> 4 norm_sd    2.165 ± 1.574 1.639 [0.127 — 6.273]  70.9
```

``` r
plot(rejection_fit, truth=truth)
```

![](tidyabc_files/figure-html/unnamed-chunk-8-1.png)
