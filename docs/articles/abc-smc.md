# ABC Sequential Monte-Carlo

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
    data1 = calculate_wasserstein(simdata$data1, obsdata$data1),
    data2 = calculate_wasserstein(simdata$data2, obsdata$data2)
  ))
}

test_scorewt = c(data1 = 1,data2 = 2)
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

![](abc-smc_files/figure-html/unnamed-chunk-5-1.png)

``` r

ggplot(tibble(data2 = test_obsdata$data2), aes(x=data2))+
  geom_histogram(,binwidth = 1)+
  xlab("B x C")
```

![](abc-smc_files/figure-html/unnamed-chunk-5-2.png)

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

smc_fit = abc_smc(
  obsdata = test_obsdata,
  priors_list = test_priors,
  sim_fn = test_simulation_fn,
  scorer_fn = test_scorer_fn,
  n_sims = 1000,
  acceptance_rate = 0.5,
  distance_method = "euclidean",
  parallel= FALSE,
  scoreweights = test_scorewt
)
#> ABC-SMC
#> SMC waves:  ■                                  0% | wave 1 ETA:  6m
#> SMC waves:  ■                                  1% | wave 2 ETA:  5m
#> SMC waves:  ■■                                 2% | wave 5 ETA:  5m
#> SMC waves:  ■■                                 3% | wave 8 ETA:  5m
#> SMC waves:  ■■                                 4% | wave 11 ETA:  5m
#> SMC waves:  ■■                                 4% | wave 14 ETA:  5m
#> Converged on wave: 17
#> SMC waves:  ■■■                                5% | wave 16 ETA:  5m

summary(smc_fit)
#> ABC SMC fit: 17 waves - (converged)
#> Parameter estimates:
#> # A tibble: 4 × 4
#> # Groups:   param [4]
#>   param      mean_sd       median_95_CrI           ESS
#>   <chr>      <chr>         <chr>                 <dbl>
#> 1 gamma_mean 5.971 ± 0.049 5.969 [5.844 — 6.097]  797.
#> 2 gamma_sd   1.972 ± 0.073 1.970 [1.792 — 2.167]  797.
#> 3 norm_mean  3.991 ± 0.171 3.990 [3.553 — 4.417]  797.
#> 4 norm_sd    1.965 ± 0.399 1.989 [0.864 — 2.845]  797.
```

``` r
plot(smc_fit, truth=truth)
```

![](abc-smc_files/figure-html/unnamed-chunk-8-1.png)

``` r
plot_convergence(smc_fit)
```

![](abc-smc_files/figure-html/unnamed-chunk-9-1.png)

``` r
plot_evolution(smc_fit,truth)
```

![](abc-smc_files/figure-html/unnamed-chunk-10-1.png)

``` r
plot_correlations(smc_fit$posteriors, truth)
```

![](abc-smc_files/figure-html/unnamed-chunk-11-1.png)
