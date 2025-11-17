# Pegasus R0 / Generation time - EXERCISE

## Simulation

## Setup siumulation —-

``` r
sim_params = list(
  # A short generation time
  mean_gt = 4,
  sd_gt = 2,
  R0 = 2,
  I0 = 10,
  # Add a longish and very variable delay to symptoms
  p_symptomatic = 0.3,
  mean_onset = 7,
  sd_onset = 5,
  # and a slightly shorter delay to observation:
  # Only symptomatic cases are observed
  p_detected_given_symptoms = 0.7,
  mean_obs = 5,
  sd_obs = 3,
  # Observation cutoff:
  T_obs = 40
)

# Run simulation ----

sim_ip = ggoutbreak::make_gamma_ip(
  median_of_mean = sim_params$mean_gt,
  median_of_sd = sim_params$sd_gt
)

sim_params$r0 = ggoutbreak::inv_wallinga_lipsitch(sim_params$R0, sim_ip)

truth = ggoutbreak::sim_branching_process(
  fn_Rt = ~ sim_params$R0,
  fn_ip = ~sim_ip,
  fn_imports = \(t) ifelse(t == 1, sim_params$I0, 0),
  max_time = 40
)
```

    ## ....................complete
    ## interfacer: development mode active (local function).

``` r
delayed = truth %>%
  ggoutbreak::sim_delay(
    p_fn = ~ sim_params$p_symptomatic,
    delay_fn = ~ ggoutbreak::rgamma2(
      .x,
      sim_params$mean_onset,
      sim_params$sd_onset
    ),
    input = "time",
    output = "symptom"
  ) %>%
  ggoutbreak::sim_delay(
    p_fn = \(t, symptom) {
      ifelse(symptom, sim_params$p_detected_given_symptoms, 0)
    },
    delay_fn = ~ ggoutbreak::rgamma2(
      .x,
      sim_params$mean_obs,
      sim_params$sd_obs
    ),
    input = "symptom_time",
    output = "observation"
  )


observed = delayed %>% dplyr::filter(observation_time < sim_params$T)

traced_contacts = observed %>%
  dplyr::semi_join(observed, by = c("infector" = "id"))
```

### Index case onset

``` r
index_case_onset = observed %>% dplyr::transmute(onset_time = floor(symptom_time))

ggplot(index_case_onset, aes(x=onset_time))+geom_histogram(binwidth = 1)
```

![](sim_example_files/figure-html/unnamed-chunk-2-1.png)

### Delay to observation

``` r
# Data

delay_distribution = observed %>% dplyr::transmute(
  obs_delay = floor(observation_time) - floor(symptom_time)
)

ggplot(delay_distribution, aes(x = obs_delay))+geom_histogram(binwidth = 1)+
  xlab("symptom to observation")
```

![](sim_example_files/figure-html/unnamed-chunk-3-1.png)

### Observed serial interval

``` r
serial_pairs = observed %>%
  inner_join(
    traced_contacts,
    by = c("id" = "infector"),
    suffix = c(".1", ".2")
  ) %>%
  transmute(
    serial_interval = floor(symptom_time.2) - floor(symptom_time.1) #order known
    # serial_interval = abs(floor(symptom_time.2) - floor(symptom_time.1)) #order uncertain
  ) 



ggplot(serial_pairs) +
  geom_histogram(aes(x = serial_interval), binwidth = 1)+
  xlab("symptom serial interval (given observed)")
```

![](sim_example_files/figure-html/unnamed-chunk-4-1.png)

``` r
obsdata = list(
  onset = as.numeric(index_case_onset$onset_time),
  diff = as.numeric(delay_distribution$obs_delay),
  si = as.numeric(serial_pairs$serial_interval)
)
```

## Model

Aim is to fit a model to all 3 aspects of the data simultaneously. The
model makes the following hard assumptions:

- constant exponential growth in infection times (hidden t_inf)
- delayed onset to symptoms - gamma distributed (t_onset = t_inf +
  onset_delay)
- delayed observation of symptoms - gamma distributed (t_obs = onset +
  obs_delay)
- secondary case delayed by generation time (hidden t_inf2 = t_inf +
  gt_delay)
- secondary case symptom onset and observation as above (t_onset_2 /
  t_obs_2)
- cases only observed if t_obs / t_obs_2 within time window (0-T)

``` math

\begin{align}
t_{max} - T_{inf} &\sim Exp(r_0) \\
\Delta T_{inf \rightarrow onset} &\sim Gamma(\mu_{onset},\sigma_{onset}) \\
\Delta T_{onset \rightarrow obs} &\sim Gamma(\mu_{obs},\sigma_{obs}) \\
\Delta T_{gt} &\sim Gamma(\mu_{gt},\sigma_{gt}) \\
\end{align}
```

``` math

\begin{align}
T_{onset} &= T_{inf} + \Delta T_{inf \rightarrow onset} \\
T_{obs} &= T_{onset} + \Delta T_{onset \rightarrow obs}\\
T_{inf_2} &= T_{inf_1} + \Delta T_{gt} \\
\Delta T_{onset_1 \rightarrow onset_2} &=  \Delta T_{gt} + \Delta T_{inf_2 \rightarrow onset_2} -  \Delta T_{inf_1 \rightarrow onset_1} \\
\end{align}
```

``` math

\begin{align}
O_1 &= I(t_0 \le T_{onset_1}, T_{obs_1} \le t_{max}) \\
O_{1,2} &= I(O_1, t_0 \le T_{onset_2}, T_{obs_2} \le t_{max})\\
T_{onset_1}|O_1 &\Rightarrow \text{primary case times}\\
\Delta T_{onset_1 \rightarrow obs_1}|O_1 &\Rightarrow \text{onset to interview delay}\\
\Delta T_{onset_1 \rightarrow onset_2}|O_{1,2} &\Rightarrow \text{onset to onset serial interval}\\
\end{align}
```

``` r
n = nrow(observed)

sim1_fn = carrier::crate(
  function(r0, mean_onset, sd_onset, mean_obs, sd_obs, mean_gt, sd_gt, R0, ...) {
    
    # Primary case infection time
    # exponentially distributed in time. Need to make sure we have enough samples 
    # before t0 observation cutoff to account for early observed cases.
    # if (r0>0) {
    #   t_inf_1 = T_obs - stats::rexp(n, r0)
    # } else {
    #   # Decling epidemic (quite unlikely)
    #   t_early = - stats::qgamma(0.99,mean_onset,sd_onset) # t starts at 0
    #   if (r0==0) {
    #     t_inf_1 = stats::runif(n,t_early,T_obs)
    #   } else {
    #     t_inf_1 = t_early+stats::rexp(n, -r0)
    #   } 
    # }
    
    t_early = - stats::qgamma(0.99,mean_onset,sd_onset) # t starts at 0
    t_inf_1 = tidyabc::rexpgrowth(n, r0, T_obs, t_early)
    
    onset_delay = ggoutbreak::rgamma2(n, mean_onset, sd_onset, convex = FALSE)
    obs_delay = ggoutbreak::rgamma2(n, mean_obs, sd_obs, convex = FALSE)
    
    t_onset_1 = t_inf_1 + onset_delay
    t_obs_1 = t_onset_1 + obs_delay
    
    # Primary case observations:
    # Onset after t0 and observed before T
    obs_1 = t_obs_1 < T_obs & t_onset_1 > 0
    
    t_inf_1 = t_inf_1[obs_1]
    t_onset_1 = t_onset_1[obs_1]
    t_obs_1 = t_obs_1[obs_1]
    
    n1 = length(t_inf_1)
    
    # Secondary case. Numbers of secondary cases are poission(R0). Could add 
    # dispersion here and fit it also
    # Only observed primary will be observed secondary so we can restrict to 
    # observed subset
    # browser()
    case_2ary = stats::rpois(n1,R0)
    index_1ary = rep(seq_along(case_2ary), case_2ary)
    n2 = length(index_1ary)
    
    gt_delay = ggoutbreak::rgamma2(n2, mean_gt, sd_gt, convex = FALSE)
    onset_delay_2 = ggoutbreak::rgamma2(n2, mean_onset, sd_onset, convex = FALSE)
    obs_delay_2 = ggoutbreak::rgamma2(n2, mean_obs, sd_obs, convex = FALSE)
    
    t_inf_2 = t_inf_1[index_1ary] + gt_delay
    t_onset_2 = t_inf_2 + onset_delay_2
    t_obs_2 = t_onset_2 + obs_delay_2
    
    # order dependent
    serial_interval = floor(t_onset_2) - floor(t_onset_1[index_1ary])
    # order independent
    # serial_interval = abs(floor(t_onset_2) - floor(t_onset_1[index_1ary]))
    
    # Secondary case observations
    obs_2 = t_obs_2 < T_obs & t_onset_2 > 0
    
    serial_interval = serial_interval[obs_2]
    t_onset_2 = t_onset_2[obs_2]
    
    return(list(
      onset = floor(t_onset_1),
      diff = floor(t_obs_1) - floor(t_onset_1),
      si = serial_interval
    ))
  },
  T_obs = sim_params$T_obs,
  n=n
)
```

``` r
scorer1_fn = function(simdata, obsdata) {
  
  onset = calculate_wasserstein(obsdata$onset, simdata$onset)
  diff = calculate_wasserstein(obsdata$diff, simdata$diff)
  si = calculate_wasserstein(obsdata$si, simdata$si)
  mad_si = abs(mean(obsdata$si) - mean(simdata$si))
  
  return(list(
    sim_onset = onset,
    sim_diff = diff,
    sim_si=si,
    sim_mad_si = mad_si
  ))
}

# make the serial interval fitting much more important:
scoreweights1 = c(
    sim_onset = 2,
    sim_diff = 0.5,
    sim_si = 3,
    sim_mad_si = 3
  )
```

``` r
test = tidyabc::test_simulation(
  sim_fn = sim1_fn, 
  scorer_fn = scorer1_fn,
  params = sim_params,
  obsdata = obsdata
  # debug=TRUE
)

# .gg_hist(test$obsdata$onset)
```

### Scoring

We are going to fit the model in a custom ABC with rejection framework

I’m assessing model fits to the data using an earth movers distance at
the level of individual’s observed times versus simulations observed
times. Given that simulations cutoff different numbers of people I have
to match the size of the simulation and observed data

### Priors

Use a set of uninformative priors. RO is derived from Wallinga-Lipsitch,
Gamma distributions hyperparameters constrained so that they are convex:

``` r
priors = priors(
  r0 ~ unif(-0.1, 0.7),
  mean_onset ~ unif(0, 12),
  sd_onset ~ unif(0, 8),
  mean_obs ~ unif(0, 12),
  sd_obs ~ unif(0, 8),
  mean_gt ~ unif(0, 12),
  sd_gt ~ unif(0, 8),
  R0 ~ (1+r0*sd_gt^2/mean_gt) ^ (mean_gt^2 / sd_gt^2),
  ~ is.finite(R0) & R0 > 0,
  ~ mean_onset > sd_onset,
  ~ mean_obs > sd_obs,
  ~ mean_gt > sd_gt
)

priors
```

    ## Parameters: 
    ## * r0: unif(min = -0.1, max = 0.7)
    ## * mean_onset: unif(min = 0, max = 12)
    ## * sd_onset: unif(min = 0, max = 8)
    ## * mean_obs: unif(min = 0, max = 12)
    ## * sd_obs: unif(min = 0, max = 8)
    ## * mean_gt: unif(min = 0, max = 12)
    ## * sd_gt: unif(min = 0, max = 8)
    ## Constraints:
    ## * is.finite(R0) & R0 > 0
    ## * mean_onset > sd_onset
    ## * mean_obs > sd_obs
    ## * mean_gt > sd_gt
    ## Derived values:
    ## * R0 = (1 + r0 * sd_gt^2/mean_gt)^(mean_gt^2/sd_gt^2)

``` r
abc_fit = abc_rejection(
  obsdata = obsdata,
  priors_list = priors,
  sim_fn = sim1_fn,
  scorer_fn = scorer1_fn,
  n_sims = 1000000,
  acceptance_rate = 0.001,
  parallel = TRUE,
  scoreweights = scoreweights1
)
```

    ## ABC rejection, 1 wave.

    ## Warning in stats::qnorm(q, -0.205013900337245, 0.0663277312263509): NaNs
    ## produced

``` r
summary(abc_fit)
```

    ## ABC rejection fit: single wave
    ## Parameter estimates:
    ## # A tibble: 8 × 4
    ## # Groups:   param [8]
    ##   param      mean_sd       median_95_CrI            ESS
    ##   <chr>      <chr>         <chr>                  <dbl>
    ## 1 R0         2.227 ± 0.545 2.115 [1.467 — 3.490]   591.
    ## 2 mean_gt    5.261 ± 1.681 4.938 [2.196 — 9.861]   591.
    ## 3 mean_obs   6.230 ± 1.903 6.032 [2.837 — 10.770]  591.
    ## 4 mean_onset 8.112 ± 1.953 8.067 [3.408 — 11.686]  591.
    ## 5 r0         0.186 ± 0.030 0.184 [0.121 — 0.268]   591.
    ## 6 sd_gt      3.696 ± 2.175 3.508 [0.240 — 7.721]   591.
    ## 7 sd_obs     4.491 ± 1.863 4.412 [1.108 — 7.779]   591.
    ## 8 sd_onset   5.094 ± 1.396 5.012 [2.149 — 7.753]   591.

``` r
plot(abc_fit, truth = sim_params)
```

![](sim_example_files/figure-html/unnamed-chunk-10-1.png)

``` r
smc_fit = abc_smc(
  obsdata = obsdata,
  priors_list = priors,
  sim_fn = sim1_fn,
  scorer_fn = scorer1_fn,
  n_sims = 8000,
  acceptance_rate = 0.5,
  #debug_errors = TRUE,
  parallel = TRUE,
  scoreweights = scoreweights1
)
```

    ## ABC-SMC

    ## Warning in stats::qnorm(q, -0.18810727283479, 0.0594144584708982): NaNs
    ## produced

    ## SMC waves:  ■                                  1% | wave 1 ETA:  5m

    ## SMC waves:  ■■                                 2% | wave 2 ETA:  5m

    ## SMC waves:  ■■                                 4% | wave 3 ETA:  5m

    ## SMC waves:  ■■■                                7% | wave 4 ETA:  5m

    ## SMC waves:  ■■■■                               9% | wave 5 ETA:  5m

    ## SMC waves:  ■■■■                              12% | wave 6 ETA:  4m

    ## SMC waves:  ■■■■■                             14% | wave 7 ETA:  4m

    ## SMC waves:  ■■■■■■                            17% | wave 8 ETA:  4m

    ## SMC waves:  ■■■■■■■                           19% | wave 9 ETA:  4m

    ## SMC waves:  ■■■■■■■■                          22% | wave 10 ETA:  4m

    ## SMC waves:  ■■■■■■■■                          25% | wave 11 ETA:  4m

    ## SMC waves:  ■■■■■■■■■                         27% | wave 12 ETA:  4m

    ## SMC waves:  ■■■■■■■■■■                        30% | wave 13 ETA:  4m

    ## SMC waves:  ■■■■■■■■■■■                       32% | wave 14 ETA:  3m

    ## SMC waves:  ■■■■■■■■■■■                       35% | wave 15 ETA:  3m

    ## SMC waves:  ■■■■■■■■■■■■                      38% | wave 16 ETA:  3m

    ## SMC waves:  ■■■■■■■■■■■■■                     40% | wave 17 ETA:  3m

    ## SMC waves:  ■■■■■■■■■■■■■■                    43% | wave 18 ETA:  3m

    ## SMC waves:  ■■■■■■■■■■■■■■■                   46% | wave 19 ETA:  3m

    ## Converged on wave: 20

``` r
summary(smc_fit)
```

    ## ABC SMC fit: 20 waves - (converged)
    ## Parameter estimates:
    ## # A tibble: 8 × 4
    ## # Groups:   param [8]
    ##   param      mean_sd       median_95_CrI            ESS
    ##   <chr>      <chr>         <chr>                  <dbl>
    ## 1 R0         2.101 ± 0.289 2.062 [1.605 — 2.723]  6328.
    ## 2 mean_gt    5.168 ± 1.119 5.040 [2.992 — 8.360]  6328.
    ## 3 mean_obs   5.468 ± 1.331 5.264 [3.129 — 9.261]  6328.
    ## 4 mean_onset 7.159 ± 1.574 7.035 [3.888 — 11.027] 6328.
    ## 5 r0         0.181 ± 0.020 0.180 [0.133 — 0.238]  6328.
    ## 6 sd_gt      3.947 ± 1.666 3.972 [0.552 — 7.404]  6328.
    ## 7 sd_obs     3.840 ± 1.479 3.649 [1.324 — 7.300]  6328.
    ## 8 sd_onset   4.662 ± 1.089 4.558 [2.456 — 7.329]  6328.

``` r
plot(smc_fit,truth = sim_params)
```

![](sim_example_files/figure-html/unnamed-chunk-12-1.png)

``` r
plot_evolution(smc_fit,truth = sim_params)
```

![](sim_example_files/figure-html/unnamed-chunk-12-2.png)

``` r
adaptive_fit = abc_adaptive(
  obsdata = obsdata,
  priors_list = priors,
  sim_fn = sim1_fn,
  scorer_fn = scorer1_fn,
  n_sims = 4000,
  acceptance_rate = 0.5,
  # debug_errors = TRUE,
  parallel = TRUE,
  scoreweights = scoreweights1
)
```

    ## ABC-Adaptive

    ## Warning in stats::qnorm(q, -0.206791663418369, 0.0654413385144222): NaNs
    ## produced

    ## Adaptive waves:  ■                                  0% | wave 1 ETA:  6m

    ## Adaptive waves:  ■                                  1% | wave 2 ETA:  5m

    ## Adaptive waves:  ■■                                 2% | wave 4 ETA:  5m

    ## Adaptive waves:  ■■                                 3% | wave 6 ETA:  5m

    ## Adaptive waves:  ■■                                 4% | wave 7 ETA:  5m

    ## Adaptive waves:  ■■■                                5% | wave 9 ETA:  5m

    ## Adaptive waves:  ■■■                                6% | wave 10 ETA:  5m

    ## Adaptive waves:  ■■■                                7% | wave 12 ETA:  5m

    ## Converged on wave: 13

``` r
summary(adaptive_fit)
```

    ## ABC adaptive fit: 13 waves - (converged)
    ## Parameter estimates:
    ## # A tibble: 8 × 4
    ## # Groups:   param [8]
    ##   param      mean_sd       median_95_CrI             ESS
    ##   <chr>      <chr>         <chr>                   <dbl>
    ## 1 R0         2.253 ± 0.674 2.111 [1.391 — 3.586]  15023.
    ## 2 mean_gt    4.707 ± 1.076 4.454 [2.492 — 9.342]  15023.
    ## 3 mean_obs   6.791 ± 1.766 6.630 [2.416 — 10.930] 15023.
    ## 4 mean_onset 8.215 ± 1.730 8.003 [3.500 — 11.609] 15023.
    ## 5 r0         0.195 ± 0.044 0.192 [0.069 — 0.323]  15023.
    ## 6 sd_gt      2.879 ± 1.420 2.858 [0.306 — 6.634]  15023.
    ## 7 sd_obs     4.060 ± 1.696 3.991 [0.723 — 7.516]  15023.
    ## 8 sd_onset   4.767 ± 1.402 4.668 [1.489 — 7.601]  15023.

``` r
plot(adaptive_fit,truth = sim_params)
```

![](sim_example_files/figure-html/unnamed-chunk-14-1.png)

``` r
plot_evolution(adaptive_fit,truth = sim_params)
```

![](sim_example_files/figure-html/unnamed-chunk-14-2.png)

``` r
plot_correlations(adaptive_fit,truth = sim_params) & ggplot2::theme(
   axis.title.y = ggplot2::element_text(angle=70,vjust=0.1)
)
```

![](sim_example_files/figure-html/unnamed-chunk-14-3.png)

``` r
plot_convergence(adaptive_fit)
```

![](sim_example_files/figure-html/unnamed-chunk-14-4.png)

``` r
plot_simulations(obsdata, adaptive_fit, sim_fn = sim1_fn)
```

![](sim_example_files/figure-html/unnamed-chunk-14-5.png)

``` r
priors2 = priors(
  r0 ~ unif(-0.1, 0.7),
  mean_onset ~ lnorm2(7, 2),
  sd_onset ~ lnorm2(5, 1),
  mean_obs ~ lnorm2(5, 1),
  sd_obs ~ lnorm2(3, 1),
  mean_gt ~ unif(0, 12),
  sd_gt ~ unif(0, 8),
  R0 ~ (1+r0*sd_gt^2/mean_gt) ^ (mean_gt^2 / sd_gt^2),
  ~ is.finite(R0) & R0 > 0,
  ~ mean_onset > sd_onset,
  ~ mean_obs > sd_obs,
  ~ mean_gt > sd_gt
)

priors2
```

    ## Parameters: 
    ## * r0: unif(min = -0.1, max = 0.7)
    ## * mean_onset: lnorm2(mean = 7, sd = 2)
    ## * sd_onset: lnorm2(mean = 5, sd = 1)
    ## * mean_obs: lnorm2(mean = 5, sd = 1)
    ## * sd_obs: lnorm2(mean = 3, sd = 1)
    ## * mean_gt: unif(min = 0, max = 12)
    ## * sd_gt: unif(min = 0, max = 8)
    ## Constraints:
    ## * is.finite(R0) & R0 > 0
    ## * mean_onset > sd_onset
    ## * mean_obs > sd_obs
    ## * mean_gt > sd_gt
    ## Derived values:
    ## * R0 = (1 + r0 * sd_gt^2/mean_gt)^(mean_gt^2/sd_gt^2)

``` r
adaptive_fit2 = abc_adaptive(
  obsdata = obsdata,
  priors_list = priors2,
  sim_fn = sim1_fn,
  scorer_fn = scorer1_fn,
  n_sims = 4000,
  acceptance_rate = 0.5,
  # debug_errors = TRUE,
  parallel = TRUE,
  scoreweights = scoreweights1
)
```

    ## ABC-Adaptive

    ## Warning in stats::qnorm(q, -0.200590853219569, 0.0666229038704867): NaNs
    ## produced

    ## Adaptive waves:  ■                                  0% | wave 1 ETA:  6m

    ## Adaptive waves:  ■                                  2% | wave 3 ETA:  5m

    ## Adaptive waves:  ■■                                 3% | wave 5 ETA:  5m

    ## Adaptive waves:  ■■                                 3% | wave 6 ETA:  5m

    ## Adaptive waves:  ■■                                 5% | wave 8 ETA:  5m

    ## Adaptive waves:  ■■■                                5% | wave 9 ETA:  5m

    ## Adaptive waves:  ■■■                                7% | wave 11 ETA:  5m

    ## Adaptive waves:  ■■■                                8% | wave 12 ETA:  5m

    ## Adaptive waves:  ■■■■                               8% | wave 13 ETA:  5m

    ## Adaptive waves:  ■■■■                               9% | wave 14 ETA:  5m

    ## Converged on wave: 15

``` r
summary(adaptive_fit2)
```

    ## ABC adaptive fit: 15 waves - (converged)
    ## Parameter estimates:
    ## # A tibble: 8 × 4
    ## # Groups:   param [8]
    ##   param      mean_sd       median_95_CrI             ESS
    ##   <chr>      <chr>         <chr>                   <dbl>
    ## 1 R0         2.158 ± 0.471 2.068 [1.474 — 3.161]  18741.
    ## 2 mean_gt    4.724 ± 0.976 4.430 [2.736 — 8.976]  18741.
    ## 3 mean_obs   4.936 ± 0.781 4.860 [3.509 — 6.764]  18741.
    ## 4 mean_onset 7.319 ± 1.483 7.055 [4.756 — 11.063] 18741.
    ## 5 r0         0.186 ± 0.035 0.184 [0.088 — 0.293]  18741.
    ## 6 sd_gt      2.909 ± 1.453 2.812 [0.354 — 6.784]  18741.
    ## 7 sd_obs     2.774 ± 0.712 2.735 [1.553 — 4.570]  18741.
    ## 8 sd_onset   4.817 ± 0.805 4.798 [3.384 — 6.608]  18741.

``` r
plot(adaptive_fit2,truth = sim_params)
```

![](sim_example_files/figure-html/unnamed-chunk-17-1.png)

``` r
plot_evolution(adaptive_fit2,truth = sim_params)
```

![](sim_example_files/figure-html/unnamed-chunk-17-2.png)

``` r
plot_correlations(adaptive_fit2,truth = sim_params)
```

![](sim_example_files/figure-html/unnamed-chunk-17-3.png)

``` r
plot_convergence(adaptive_fit2)
```

![](sim_example_files/figure-html/unnamed-chunk-17-4.png)

``` r
plot_simulations(obsdata, adaptive_fit2, sim_fn = sim1_fn)
```

![](sim_example_files/figure-html/unnamed-chunk-17-5.png)
