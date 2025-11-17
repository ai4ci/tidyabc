# ABC Adaptive

``` r
# library(tidyabc)

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
devtools::load_all()
#> ℹ Loading tidyabc
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

![](abc-adaptive_files/figure-html/unnamed-chunk-5-1.png)

``` r

ggplot(tibble(data2 = test_obsdata$data2), aes(x=data2))+
  geom_histogram(,binwidth = 1)+
  xlab("B x C")
```

![](abc-adaptive_files/figure-html/unnamed-chunk-5-2.png)

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

adaptive_fit = abc_adaptive(
  obsdata = test_obsdata,
  priors_list = test_priors,
  sim_fn = test_simulation_fn,
  scorer_fn = test_scorer_fn,
  n_sims = 1000,
  acceptance_rate = 0.5,
  parallel= FALSE,
  bw=0.05
)
#> ABC-Adaptive
#> Adaptive waves:  ■                                  0% | wave 1 ETA:  7m
#> Adaptive waves:  ■                                  1% | wave 5 ETA:  5m
#> Adaptive waves:  ■■                                 2% | wave 9 ETA:  5m
#> Adaptive waves:  ■■                                 3% | wave 13 ETA:  5m
#> Converged on wave: 15
#> Adaptive waves:  ■■                                 4% | wave 14 ETA:  5m

summary(adaptive_fit)
#> ABC adaptive fit: 15 waves - (converged)
#> Parameter estimates:
#> # A tibble: 4 × 4
#> # Groups:   param [4]
#>   param      mean_sd       median_95_CrI           ESS
#>   <chr>      <chr>         <chr>                 <dbl>
#> 1 gamma_mean 5.982 ± 0.055 5.967 [5.810 — 6.147] 1677.
#> 2 gamma_sd   2.056 ± 0.097 2.059 [1.742 — 2.295] 1677.
#> 3 norm_mean  3.978 ± 0.142 3.959 [3.613 — 4.424] 1677.
#> 4 norm_sd    1.990 ± 0.390 2.034 [0.651 — 2.849] 1677.
```

``` r
plot(adaptive_fit, truth=truth)
```

![](abc-adaptive_files/figure-html/unnamed-chunk-8-1.png)

``` r
plot_convergence(adaptive_fit)
```

![](abc-adaptive_files/figure-html/unnamed-chunk-9-1.png)

``` r
plot_evolution(adaptive_fit,truth)
```

![](abc-adaptive_files/figure-html/unnamed-chunk-10-1.png)

``` r
plot_correlations(adaptive_fit$posteriors, truth)
```

![](abc-adaptive_files/figure-html/unnamed-chunk-11-1.png)
