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
#>  [1]  1.540629  5.824354  2.998694  6.382562  3.991957  3.551254 -0.242528
#>  [8]  4.163835  5.073847  5.011748
#> 
#> $B
#>  [1] 3.690281 4.677996 2.640020 3.331717 1.540111 1.555870 2.184533 1.562332
#>  [9] 2.817359 2.469767
#> 
example_scorer_fn(
  example_obsdata(),
  example_sim_fn(3,2,1)
)
#> $A
#> [1] 0.9761397
#> 
#> $B
#> [1] 1.981469
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
#> 1 mean  5.002 ± 0.185 5.004 [4.556 — 5.461]  78.1
#> 2 sd1   2.288 ± 0.537 2.224 [1.339 — 3.983]  78.1
#> 3 sd2   1.189 ± 0.312 1.150 [0.660 — 2.095]  78.1
example_smc_fit()
#> ABC SMC fit: 7 waves - (not yet converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.978 ± 0.036 4.979 [4.886 — 5.063]  260.
#> 2 sd1   1.998 ± 0.089 1.999 [1.786 — 2.217]  260.
#> 3 sd2   0.995 ± 0.045 0.993 [0.894 — 1.128]  260.
example_adaptive_fit()
#> ABC-Adaptive
#> Adaptive waves:  ■■■■■■■■■                         25% | wave 2 ETA:  4s
#> Adaptive waves:  ■■■■■■■■■■■■■                     39% | wave 3 ETA:  3s
#> Adaptive waves:  ■■■■■■■■■■■■■■■■■                 52% | wave 4 ETA:  2s
#> Adaptive waves:  ■■■■■■■■■■■■■■■■■■■■■             66% | wave 5 ETA:  2s
#> Adaptive waves:  ■■■■■■■■■■■■■■■■■■■■■■■■■         79% | wave 6 ETA:  1s
#> Adaptive waves:  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     92% | wave 7 ETA:  0s
#> ABC adaptive fit: 8 waves - (not yet converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.979 ± 0.020 4.982 [4.917 — 5.041]  349.
#> 2 sd1   2.034 ± 0.045 2.030 [1.854 — 2.174]  349.
#> 3 sd2   1.002 ± 0.028 1.000 [0.899 — 1.077]  349.
```
