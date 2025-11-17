# Example generators

These are a set of internally cached functions to support examples. They
are exported as internal functions so that the examples can run
correctly and cache their output to prevent excessive repetition of the
code examples.

## Usage

``` r
example_sim_fn(mean, sd1, sd2)

example_scorer_fn(simdata, obsdata)

example_obs()

example_truth()

example_obsdata()

example_priors_list()

example_rejection_fit()

example_smc_fit()

example_adaptive_fit()
```

## Arguments

- mean:

  the example mean

- sd1:

  the example normal sd (A)

- sd2:

  the example gamma sd (B)

- simdata:

  output of example_sim_fn()

- obsdata:

  e.g. example_obsdata()

## Value

a list of A and B of samples from a normal and gamma

a list of scores

output of
[`test_simulation()`](https://ai4ci.github.io/tidyabc/reference/test_simulation.md)

example `truth` parameter

example `obsdata` parameter

example `priors_list` parameter

example `abc_rejection` output

example `abc_smc` output

example `abc_adaptive` output

## Functions

- `example_sim_fn()`: Example simulation function

- `example_scorer_fn()`: Example simulation function

- `example_obs()`: Example output of
  [`test_simulation()`](https://ai4ci.github.io/tidyabc/reference/test_simulation.md)

- `example_truth()`: Example for `truth` parameter

- `example_obsdata()`: Example for `obsdata` parameter

- `example_priors_list()`: Example for `priors_list` parameter

- `example_rejection_fit()`: Example `abc_rejection` output

- `example_smc_fit()`: Example `abc_smc` output

- `example_adaptive_fit()`: Example `abc_adaptive` output

## Examples

``` r
example_sim_fn(3,2,1) %>% lapply(head,10)
#> $A
#>  [1]  6.395490  1.142441  5.519913  1.466012  6.805626  2.421382 -2.136159
#>  [8]  2.435590  1.426870 -1.350074
#> 
#> $B
#>  [1] 4.507524 2.467182 3.138320 2.659531 1.670072 5.975053 2.879220 3.255144
#>  [9] 2.768422 3.060035
#> 
example_scorer_fn(
  example_obsdata(),
  example_sim_fn(3,2,1)
)
#> $A
#> [1] 1.27649
#> 
#> $B
#> [1] 2.489633
#> 
example_obs() %>% rapply(head,n=10, how="replace")
#> $obsdata
#> $obsdata$A
#>  [1] 3.879049 4.539645 8.117417 5.141017 5.258575 8.430130 5.921832 2.469878
#>  [9] 3.626294 4.108676
#> 
#> $obsdata$B
#>  [1] 3.963790 4.935891 5.197293 3.450608 5.984255 5.150334 7.583826 5.601786
#>  [9] 4.467521 4.746891
#> 
#> 
#> $truth
#> mean  sd1  sd2 
#>    5    2    1 
#> 
#> $sim_fn
#> [1] "function (mean, sd1, sd2) "                                                    
#> [2] "{"                                                                             
#> [3] "    return(list(A = stats::rnorm(1000, mean, sd1), B = tidyabc::rgamma2(1000, "
#> [4] "        mean, sd2)))"                                                          
#> [5] "}"                                                                             
#> 
#> $obsscores
#> $obsscores$A
#> [1] 0
#> 
#> $obsscores$B
#> [1] 0
#> 
#> 
#> $scorer_fn
#> [1] "function (simdata, obsdata) "                                      
#> [2] "{"                                                                 
#> [3] "    return(list(A = tidyabc::calculate_wasserstein(obsdata$A, "    
#> [4] "        simdata$A), B = tidyabc::calculate_wasserstein(obsdata$B, "
#> [5] "        simdata$B)))"                                              
#> [6] "}"                                                                 
#> 
example_truth()
#> mean  sd1  sd2 
#>    5    2    1 
example_obsdata() %>% lapply(head,10)
#> $A
#>  [1] 3.879049 4.539645 8.117417 5.141017 5.258575 8.430130 5.921832 2.469878
#>  [9] 3.626294 4.108676
#> 
#> $B
#>  [1] 3.963790 4.935891 5.197293 3.450608 5.984255 5.150334 7.583826 5.601786
#>  [9] 4.467521 4.746891
#> 
example_priors_list()
#> Parameters: 
#> * mean: unif(min = 0, max = 10)
#> * sd1: unif(min = 0, max = 5)
#> * sd2: unif(min = 0, max = 5)
#> Constraints:
#> * mean > sd2
example_rejection_fit()
#> ABC rejection, 1 wave.
#> ABC rejection fit: single wave
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  5.014 ± 0.179 5.034 [4.591 — 5.429]  75.3
#> 2 sd1   2.425 ± 0.556 2.322 [1.329 — 3.934]  75.3
#> 3 sd2   1.254 ± 0.307 1.212 [0.727 — 2.024]  75.3
example_smc_fit()
#> ABC SMC fit: 8 waves - (converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.974 ± 0.037 4.975 [4.882 — 5.074]  241.
#> 2 sd1   2.001 ± 0.078 2.001 [1.812 — 2.241]  241.
#> 3 sd2   0.992 ± 0.037 0.991 [0.909 — 1.102]  241.
example_adaptive_fit()
#> ABC-Adaptive
#> Adaptive waves:  ■■■■■■■                           22% | wave 2 ETA:  4s
#> Adaptive waves:  ■■■■■■■■■■■■■■                    44% | wave 4 ETA:  3s
#> Converged on wave: 8
#> Adaptive waves:  ■■■■■■■■■■■■■■■■■■■■■■■■■         79% | wave 7 ETA:  1s
#> ABC adaptive fit: 8 waves - (converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.979 ± 0.026 4.978 [4.891 — 5.075]  379.
#> 2 sd1   2.053 ± 0.064 2.072 [1.850 — 2.211]  379.
#> 3 sd2   1.006 ± 0.036 1.008 [0.893 — 1.099]  379.
```
