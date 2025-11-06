## Workflow wrappers ----

#' Perfom simple ABC rejection algorithm
#'
#' This function will execute a simulation for a random selection of parameters
#' and identify the best matching `acceptance_rate` percent, as defined by the
#' summary distance metric. A large number of simulations and a low acceptance
#' rate are best here.
#'
#' @inheritParams tidyabc_common
#'
#' @returns a list containing 3 items. A flag assessing convergence, a list of
#'   summary metrics for each wave (`waves`), and `posteriors` - the filtered
#'   parameter posteriors including component and summary score columns, plus
#'   weights. Optionally this will also include a simulation data list column
#'   containing the raw simulated data.
#' @export
#' @concept workflow
#' @examples
#'
#' fit = abc_rejection(
#'   example_obsdata(),
#'   example_priors_list(),
#'   example_sim_fn,
#'   example_scorer_fn,
#'   n_sims = 10000,
#'   acceptance_rate = 0.01
#' )
#'
#' summary(fit)
#'
abc_rejection = function(
  obsdata,
  priors_list,
  sim_fn,
  scorer_fn,
  n_sims,
  acceptance_rate,
  converged_fn = default_termination_fn(),
  obsscores = NULL,
  distance_method = "euclidean",
  keep_simulations = FALSE,
  seed = NULL,
  parallel = FALSE
) {
  message("ABC rejection, 1 wave.")

  prev_sim_df = .sample_priors(priors_list, n_sims)
  summ_df = posterior_summarise(prev_sim_df, priors_list) %>%
    dplyr::mutate(wave = 0)

  sim_df = .abc_do_one(
    obsdata = obsdata,
    sim_df = prev_sim_df,
    sim_fn = sim_fn,
    scorer_fn = scorer_fn,
    obsscores = obsscores,
    distance_method = distance_method,
    keep_simulations = keep_simulations,
    seed = seed,
    parallel = parallel
  )

  cutoff = stats::quantile(sim_df$abc_summary_distance, acceptance_rate)
  epsilon = stats::quantile(sim_df$abc_summary_distance, acceptance_rate / 2)

  # Hard SMC
  sim_df = sim_df %>%
    dplyr::filter(abc_summary_distance < cutoff)

  sim_df = sim_df %>%
    .calculate_weights_wave_one(
      epsilon = epsilon
    )

  # browser()
  # ggplot(sim_df)+geom_point(aes(x=abc_summary_distance,y = abc_weight))

  metric = .compare_waves(
    sim_df = sim_df,
    prev_sim_df = prev_sim_df,
    priors_list = priors_list
  ) %>%
    dplyr::mutate(wave = 1)

  converged = converged_fn(
    summary = metric$summary[[1]],
    per_param = metric$per_param[[1]]
  )

  summ_df = dplyr::bind_rows(
    summ_df,
    posterior_summarise(sim_df, priors_list) %>% dplyr::mutate(wave = 1)
  )

  return(new_abc_fit(
    type = "rejection",
    iterations = 1,
    converged = converged,
    priors_list = priors_list,
    wave_df = metric,
    summ_df = summ_df,
    sim_df = sim_df
  ))
}

#' Perform ABC sequential Monte Carlo fitting
#'
#' This function will execute a simulation for a random selection of parameters.
#' Based on the `acceptance_rate` it will reject a proportion of the results.
#' The remaining results are weighted (using a kernel with a tolerance
#' equivalent to half the acceptance rate). Weighted parameter particles
#' generate proposals for further waves but a particle perturbation. Waves
#' are executed until a maximum is reached or the results converge sufficiently
#' that the changes between waves are small. A relatively small number of
#' simulations may be attempted with a high acceptance rate, over multiple waves.
#'
#' @inheritParams tidyabc_common
#' @param max_waves the maximum number of ABC iterations to do before admitting
#'   defeat. The number of simulations executed will be `n_sims*max_waves`.
#'   Fewer waves may be performed if the algorithm converges early.
#'
#' @returns a list containing 3 items. A flag assessing convergence, a list
#'   of summary metrics for each wave (`waves`), and `posteriors` - the last
#'   wave filtered parameter posteriors including component and summary score
#'   columns, plus weights.
#' @export
#' @concept workflow
#' @examples
#'
#' fit = abc_smc(
#'   obsdata = example_obsdata(),
#'   priors_list = example_priors_list(),
#'   sim_fn = example_sim_fn,
#'   scorer_fn = example_scorer_fn,
#'   n_sims = 1000,
#'   acceptance_rate = 0.25,
#'   max_waves = 5,
#'   parallel = FALSE,
#'   allow_continue = FALSE
#' )
#'
#' summary(fit)
#'
abc_smc = function(
  obsdata,
  priors_list,
  sim_fn,
  scorer_fn,
  n_sims,
  acceptance_rate,
  kernel_t = 0.2,
  max_waves = ceiling(log(0.001) / log(acceptance_rate)),
  converged_fn = default_termination_fn(),
  obsscores = NULL,
  distance_method = "euclidean",
  seed = NULL,
  parallel = FALSE,
  allow_continue = interactive()
) {
  if (!is.null(seed)) {
    seed = set.seed(seed)
    on.exit(set.seed(seed), add = TRUE)
  }

  # Generate a set of initial proposals
  sim_df = .sample_priors(priors_list, n_sims)
  # browser()
  message("ABC-SMC, waves: ", max_waves)
  bar_id = cli::cli_progress_bar(name = "SMC waves: ", total = max_waves)
  prev_sim_df = NULL
  wave_df = NULL
  wave1_cov = NULL
  summ_df = posterior_summarise(sim_df, priors_list) %>% dplyr::mutate(wave = 0)
  converged = FALSE
  i = 0
  remaining_waves = max_waves
  while (remaining_waves > 0) {
    remaining_waves = remaining_waves - 1
    i = i + 1

    # Execute simulation
    sim_df = .abc_do_one(
      obsdata = obsdata,
      sim_df,
      sim_fn = sim_fn,
      scorer_fn = scorer_fn,
      obsscores = obsscores,
      distance_method = distance_method,
      keep_simulations = FALSE,
      seed = NULL,
      parallel = parallel,
      wave1_cov = wave1_cov
    )

    wave1_cov = attr(sim_df$abc_summary_distance, "cov")

    cutoff = stats::quantile(sim_df$abc_summary_distance, acceptance_rate)
    epsilon = stats::quantile(sim_df$abc_summary_distance, acceptance_rate / 2)

    # Hard SMC
    sim_df = sim_df %>%
      dplyr::filter(abc_summary_distance <= cutoff)

    sim_df = sim_df %>%
      .calculate_weights_smc(
        epsilon = epsilon,
        prev_sim_df = prev_sim_df,
        kernel_t
      )

    metric = .compare_waves(
      sim_df = sim_df,
      prev_sim_df = prev_sim_df,
      priors_list = priors_list
    ) %>%
      dplyr::mutate(wave = i)

    wave_df = dplyr::bind_rows(wave_df, metric)

    summ_df = dplyr::bind_rows(
      summ_df,
      posterior_summarise(sim_df, priors_list) %>% dplyr::mutate(wave = i)
    )

    prev_sim_df = sim_df

    if (
      converged_fn(
        summary = metric$summary[[1]],
        per_param = metric$per_param[[1]]
      )
    ) {
      converged = TRUE
      message("Converged on wave: ", i)
      break
    }

    if (remaining_waves == 0) {
      if (allow_continue) {
        message("Current estimates: ")
        cat(.format_summ(summ_df), sep = "\n")
        if (
          utils::askYesNo(
            sprintf(
              "Still converging. Do you want another %d waves?",
              max_waves
            )
          )
        ) {
          remaining_waves = max_waves
        } else {
          break
        }
      } else {
        break
      }
    }

    sim_df = .sample_using_perturbation(
      n_sims,
      prev_sim_df,
      priors_list,
      kernel_t
    )

    cli::cli_progress_update(id = bar_id, total = i + remaining_waves)
  }

  cli::cli_progress_done(id = bar_id)

  return(new_abc_fit(
    type = "SMC",
    iterations = i,
    converged = converged,
    priors_list = priors_list,
    wave_df = wave_df,
    summ_df = summ_df,
    sim_df = sim_df
  ))
}


#' Perform ABC sequential adaptive fitting
#'
#' This function will execute a simulation for a random selection of parameters.
#' Based on the `acceptance_rate` it will reject a proportion of the results.
#' The remaining results are weighted (using a kernel with a tolerance
#' equivalent to half the acceptance rate). Empirical distributions are fitted
#' to weighted parameter particles and from these proposals are generated for
#' further waves by fresh sampling. Waves are executed until a maximum is
#' reached or the results converge sufficiently that the changes between waves
#' are small. A relatively small number of simulations may be attempted with a
#' high acceptance rate, over multiple waves.
#'
#' @inheritParams tidyabc_common
#' @param max_waves the maximum number of ABC iterations to do before admitting
#'   defeat. The number of simulations executed will be `n_sims*max_waves`.
#'   Fewer waves may be performed if the algorithm converges early.
#' @inheritParams empirical
#'
#' @returns a list containing 3 items. A flag assessing convergence, a list
#'   of summary metrics for each wave (`waves`), and `posteriors` - the last
#'   wave filtered parameter posteriors including component and summary score
#'   columns, plus weights. Optionally this will also include a simulation data
#'   list column containing the raw simulated data.
#' @export
#' @concept workflow
#' @examples
#'
#' fit = abc_adaptive(
#'   obsdata = example_obsdata(),
#'   priors_list = example_priors_list(),
#'   sim_fn = example_sim_fn,
#'   scorer_fn = example_scorer_fn,
#'   n_sims = 1000,
#'   acceptance_rate = 0.25,
#'   max_waves = 5,
#'   parallel = FALSE,
#'   allow_continue = FALSE
#' )
#'
#' summary(fit)
#'
abc_adaptive = function(
  obsdata,
  priors_list,
  sim_fn,
  scorer_fn,
  n_sims,
  acceptance_rate,
  kernel_t = 0.2,
  max_waves = ceiling(log(0.001) / log(acceptance_rate)),
  converged_fn = default_termination_fn(),
  obsscores = NULL,
  distance_method = "euclidean",
  seed = NULL,
  knots = NULL,
  parallel = FALSE,
  max_recover = 3,
  allow_continue = interactive()
) {
  if (!is.null(seed)) {
    seed = set.seed(seed)
    on.exit(set.seed(seed), add = TRUE)
  }

  # browser()
  message("ABC-Adaptive, waves: ", max_waves)
  bar_id = cli::cli_progress_bar(name = "Adaptive waves: ", total = max_waves)

  # Setup first wave

  wave_df = NULL
  wave1_cov = NULL
  # N.b. could be refactored inside the loop?
  # Generate a set of initial proposals
  sim_df = .sample_priors(priors_list, n_sims)
  summ_df = posterior_summarise(sim_df, priors_list) %>% dplyr::mutate(wave = 0)
  converged = FALSE
  i = 0
  remaining_waves = max_waves
  proposal_list = priors_list
  prev_sim_df = sim_df
  recover = 0

  while (remaining_waves > 0) {
    remaining_waves = remaining_waves - 1
    i = i + 1

    while (recover <= max_recover) {
      # Execute simulation
      sim_df = .abc_do_one(
        obsdata = obsdata,
        sim_df,
        sim_fn = sim_fn,
        scorer_fn = scorer_fn,
        obsscores = obsscores,
        distance_method = distance_method,
        keep_simulations = FALSE,
        seed = NULL,
        parallel = parallel,
        wave1_cov = wave1_cov
      )

      cutoff = stats::quantile(sim_df$abc_summary_distance, acceptance_rate)
      epsilon = stats::quantile(
        sim_df$abc_summary_distance,
        acceptance_rate / 2
      )

      # Hard SMC
      sim_df = sim_df %>%
        dplyr::filter(abc_summary_distance <= cutoff)

      sim_df = sim_df %>%
        .calculate_weights_adaptive(
          priors_list = priors_list,
          epsilon = epsilon,
          proposal_list = proposal_list
        )

      ESS = 1 / sum(sim_df$abc_weight^2)
      # Exit the retry loop if ESS is big enough
      if (ESS > 200) {
        break
      }

      # ESS is too small. we will try and recover by
      # rerunning the wave with a larger sample size.
      sim_df = prev_sim_df
      recover = recover + 1
      message("Effective sample size has reduced below 200.")
      if (recover > max_recover) {
        warning("No more retries. Aborting.")
        return(new_abc_fit(
          type = "adaptive",
          iterations = i,
          converged = converged,
          priors_list = priors_list,
          wave_df = wave_df,
          summ_df = summ_df,
          sim_df = prev_sim_df
        ))
      } else {
        message("Attempting recovery with larger sample size.")
        n_sims = n_sims * 2
        sim_df = .sample_priors(proposal_list, n_sims)
      }

      # re-enter recovery loop
    }

    # End of recovery loop:
    # Extract the covariance if it has been calculated
    # This stays the same throughout all waves
    wave1_cov = attr(sim_df$abc_summary_distance, "cov")

    metric = .compare_waves(
      sim_df = sim_df,
      prev_sim_df = prev_sim_df,
      priors_list = priors_list
    ) %>%
      dplyr::mutate(
        wave = i,
        # keep the proposal distribution
        proposal_distribution = list(proposal_list)
      )

    wave_df = dplyr::bind_rows(wave_df, metric)

    summ_df = dplyr::bind_rows(
      summ_df,
      posterior_summarise(sim_df, priors_list) %>% dplyr::mutate(wave = i)
    )

    # Store current wave as previous:
    prev_sim_df = sim_df

    # Check whether current wave meets convergence criteria
    if (
      converged_fn(
        summary = metric$summary[[1]],
        per_param = metric$per_param[[1]]
      )
    ) {
      converged = TRUE
      message("Converged on wave: ", i)
      break
    }

    if (remaining_waves == 0) {
      if (allow_continue) {
        message("Current estimates: ")
        cat(.format_summ(summ_df), sep = "\n")
        if (
          utils::askYesNo(
            sprintf(
              "Still converging. Do you want another %d waves?",
              max_waves
            )
          )
        ) {
          remaining_waves = max_waves
        } else {
          break
        }
      } else {
        break
      }
    }

    proposal_list = posterior_fit_empirical(
      sim_df,
      priors_list = proposal_list, # could be proposal_list...?
      knots = knots
    )

    sim_df = .sample_priors(proposal_list, n_sims)
    # browser()

    cli::cli_progress_update(id = bar_id, total = i + remaining_waves)
  }

  cli::cli_progress_done(id = bar_id)

  return(new_abc_fit(
    type = "adaptive",
    iterations = i,
    converged = converged,
    priors_list = priors_list,
    wave_df = wave_df,
    summ_df = summ_df,
    sim_df = prev_sim_df
  ))
}

## Configuration helpers ----

#' Set up default convergence criteria for SMC and adaptive ABC
#'
#' Convergence is assessed on firstly whether the central estimate of the parameters
#' being assessed is stable, and not changing from one wave to the next, and
#' secondly if the 95 percent credible interval is stable between waves. If the
#' parameter central estimate is stationary but the credible intervals are still
#' dropping then continuing simulation may get better estimates of confidence.
#'
#' @param stability how close do sequential estimates need to be before
#'   declaring them as a good set of parameter estimates. This is in the units
#'   of the parameter. If this is given as a single number it applies to all
#'   parameters equally, alternatively a named vector can be used to set parameter
#'   specific cutoffs.
#' @param confidence how stable do the 95% confidence intervals need to be
#'   before we are happy with the parameter estimates. This is in the scale of
#'   the parameters, but represents a change in IQR from wave to wave. If this
#'   is given as a single number it applies to all parameters equally,
#'   alternatively a named vector can be used to set parameter specific cutoffs.
#'
#' @returns a function that specifies the convergence.
#' @export
#' @concept workflow
#'
#' @examples
#'
#' # A more permissive definition of convergence has
#' # less strict stability criteria (sequential estimates varying by less than 5%)
#' # and confidence intervals not changing by more than 1 unit between waves.)
#' check = default_termination_fn(0.05, 1.0)
#'
#'
#' fit = example_smc_fit()
#'
#' # This is performed as an integral part of the SMC and adaptive
#' # fitting and is here only for example
#' last_wave_metrics = tail(fit$waves,1)
#' converged = check(
#'   last_wave_metrics$summary[[1]],
#'   last_wave_metrics$per_param[[1]]
#' )
#'
#' if (isTRUE(converged)) print("Converged (permissive definition)")
#'
default_termination_fn = function(stability = 0.01, confidence = 0.1) {
  function(summary, per_param) {
    if (length(stability) != nrow(per_param)) {
      stability = rep(stability, nrow(per_param))
    }
    if (is.null(names(stability))) {
      names(stability) = per_param$param
    }
    if (length(confidence) != nrow(per_param)) {
      confidence = rep(confidence, nrow(per_param))
    }
    if (is.null(names(confidence))) {
      names(confidence) = per_param$param
    }

    isTRUE(
      all(per_param$rel_mean_change < stability[per_param$param]) &&
        all(per_param$IQR_95_redn < confidence[per_param$param])
    )
  }
}

#' Run the simulation for one set of parameters
#'
#' @inheritParams tidyabc_common
#' @param ... simulation parameters, must be named
#'
#' @returns a list containing the parameters as `truth` and an
#'   instance of the simulation as `obsdata`
#' @export
#' @concept workflow
test_simulation = function(sim_fn, scorer_fn, ..., seed = NULL) {
  if (!is.null(seed)) {
    seed = set.seed(seed)
    on.exit(set.seed(seed), add = TRUE)
  }

  sim_fn = .crate_sim_fn(sim_fn)
  stopifnot(all(nzchar(names(rlang::list2(...)))))
  test = sim_fn(...)
  params = rlang::list2(...)

  scorer_fn = .crate_scorer_fn(scorer_fn, obsdata = test)
  scores = scorer_fn(test)

  return(list(
    obsdata = test,
    truth = unlist(params),
    sim_fn = format(unclass(rlang::fn_env(sim_fn)$sim_fn)),
    obsscores = scores,
    scorer_fn = format(unclass(rlang::fn_env(scorer_fn)$scorer_fn))
  ))
}


## Functions for working with posteriors ----

#' Generate a set of samples from selected posteriors
#'
#' Once an ABC model fitting is complete the simulation data is generally
#' only one possible realisation of the parameters, which has been selected for
#' closeness. To properly compare the output with the observed data we need a
#' set of posterior re-samples, which are selected from posteriors according to
#' importance.
#'
#' @inheritParams tidyabc_common
#' @param n_resamples the number of resamples for each parameter combination.
#' @param max_samples the maximum total number of resamples to pick.
#'
#' @returns a dataframe of the posteriors with an `abc_sim` list column
#'   containing the output of `sim_fn` called with the parameters in that row.
#' @export
#' @concept workflow
#'
#' @examples
#'
#' fit = example_adaptive_fit()
#'
#' sample_df = posterior_resample(
#'   fit$posteriors,
#'   sim_fn = example_sim_fn
#' )
#'
#' # the fitted simulations are in the `abc_sim` column
#' sim1 = sample_df$abc_sim[[1]]
#'
#' sim1 %>% lapply(head, 10)
#'
posterior_resample = function(
  posteriors_df,
  sim_fn,
  n_resamples = 1,
  seed = NULL,
  parallel = FALSE,
  max_samples = 200
) {
  if (!is.null(seed)) {
    seed = set.seed(seed)
    on.exit(set.seed(seed), add = TRUE)
  }

  if (nrow(posteriors_df) * n_resamples > max_samples) {
    posteriors_df = posteriors_df %>%
      dplyr::slice_sample(
        weight_by = abc_weight,
        n = max_samples %/% n_resamples
      )
  }

  sim_crate = .crate_sim_fn(sim_fn)

  sim_df = .abc_do_simulation_and_scoring(
    sim_df = posteriors_df,
    sim_crate = sim_crate,
    n_resamples = n_resamples,
    keep_simulations = TRUE,
    scorer_crate = NULL,
    parallel = parallel
  )

  # exclude the mvn columns
  sim_df %>% dplyr::select(-dplyr::starts_with("abc_mvn_"))
}


#' Calculate a basket of summaries from a weighted list of posterior samples
#'
#' @inheritParams tidyabc_common
#' @returns a dataframe indexed by parameter with useful summary metrics.
#' @export
#' @concept workflow
#'
#' @examples
#'
#' fit = example_adaptive_fit()
#' summ = posterior_summarise(fit$posteriors)
#'
#' summ %>% dplyr::glimpse()
#'
posterior_summarise = function(
  posteriors_df,
  priors_list = NULL,
  p = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
) {
  weights = suppressWarnings(posteriors_df$abc_weight)
  if (is.null(weights)) {
    weights = rep(1 / nrow(posteriors_df), nrow(posteriors_df))
  }
  df = posteriors_df %>%
    dplyr::select(-dplyr::starts_with("abc"))

  summ_df = df %>%
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "param",
      values_to = "value"
    ) %>%
    dplyr::group_by(param) %>%
    dplyr::summarise(
      mean = stats::weighted.mean(value, w = weights),
      sd = sqrt(Hmisc::wtd.var(value, weights = weights, normwt = TRUE)),
      quantiles = list(
        tibble::tibble(
          p = p,
          x = unname(
            Hmisc::wtd.quantile(
              value,
              weights = weights,
              probs = p,
              normwt = TRUE
            )
          )
        )
      ),
      q.0.025 = Hmisc::wtd.quantile(
        value,
        weights = weights,
        probs = 0.025,
        normwt = TRUE
      ),
      q.0.5 = Hmisc::wtd.quantile(
        value,
        weights = weights,
        probs = 0.5,
        normwt = TRUE
      ),
      q.0.975 = Hmisc::wtd.quantile(
        value,
        weights = weights,
        probs = 0.975,
        normwt = TRUE
      ),
      # IQR_95 = q.0.975 - q.0.025,
      density = c(empirical(
        value,
        weights = weights,
        name = dplyr::cur_group()$param,
        link = if (is.null(priors_list)) {
          "ident"
        } else {
          priors_list[[dplyr::cur_group()$param]]
        }
      )),
      ESS = 1 / sum(weights^2)
    )

  return(summ_df)
}


#' Fit empirical distribution to posterior samples for generating more waves
#'
#' This function allows "updating" of the prior with an empirical posterior
#' distribution which will retain the bounds of the prior, and is effectively a
#' spline based transform of the prior distribution, based on the density of
#' data in prior space. This gives a clean density when data is close to a
#' prior distribution limit and work better than a standard density for
#'
#' @inheritParams tidyabc_common
#'
#' @returns a list of `dist_fns` approximating the distribution of the posterior
#'   samples, adhering to prior support if supplied.
#' @export
#' @concept workflow
#' @examples
#'
#' fit = example_smc_fit()
#' proposals = posterior_fit_empirical(fit$posteriors)
#'
#' proposals
#'
posterior_fit_empirical = function(
  posteriors_df,
  priors_list = list(),
  knots = NULL
) {
  weights = suppressWarnings(posteriors_df$abc_weight) #maybe null
  nms = colnames(posteriors_df)
  nms = nms[!grepl("^abc_", nms)]
  nms = nms[nms != ""]

  # Get particles in proposal space
  # posteriors_df = posteriors_df %>%
  #   dplyr::select(dplyr::all_of(nms))

  # TODO: generate appropriate link function
  # currently we are using the prior as a link function. This is OK but does
  # mean the transformations will get nested after multiple waves. This
  # enforces constraints on the priors at the cost of complexity. It's also
  # possible that the empirical transformation will be unfavourable in prior
  # space. There is a potential performance bottleneck here. A simpler approach
  # to identify a standard link from the distribution support would be simpler.
  tmp = lapply(nms, function(nm) {
    col = posteriors_df[[nm]]
    prior = priors_list[[nm]]
    if (is.null(prior)) {
      prior = "ident"
    }
    # limits = prior$q(c(0, 1))
    empirical(col, knots = knots, weights = weights, link = prior, name = nm)
  })
  names(tmp) = nms

  # Get particles in MVN space
  theta = posteriors_df %>%
    dplyr::select(dplyr::starts_with("abc_mvn_")) %>%
    as.matrix()
  # weighted correlation:
  cov = stats::cov.wt(theta, wt = weights, method = "ML")$cov
  # browser()

  priors_list = utils::modifyList(priors_list, tmp)
  return(
    structure(
      as.dist_fns_list(tmp),
      cor = stats::cov2cor(cov)
    )
  )
}

# # Fit truncated normal distribution to posterior samples for generating more waves
# #
# # This function will match mean and SD of untruncated distribution, and allows updating of the prior which will retain derived value
# # definitions and constraint checks from the prior.
# #
# # @inheritParams tidyabc_common
# #
# # @returns a list of `dist_fns` approximating the distribution of the posterior
# #   samples, plus constraints and derived values if they were supplied.
# # @export
# posterior_fit_tnorm = function(
#   posteriors_df,
#   priors_list = list(),
#   knots = NULL
# ) {
#   weights = suppressWarnings(posteriors_df$abc_weights) #maybe null
#
#   nms = names(priors_list)
#   nms = nms[nms != ""]
#   posteriors_df = posteriors_df %>%
#     dplyr::select(dplyr::all_of(nms))
#
#   tmp = lapply(colnames(posteriors_df), function(nm) {
#     col = posteriors_df[[nm]]
#     fit = fitdistrplus::fitdist(col, "norm")
#     dfn = as.dist_fns(fit)
#     prior = priors_list[[nm]]
#     if (!is.null(prior)) {
#       low = prior$q(0)
#       high = prior$q(1)
#       if (!is.finite(low)) {
#         low = NA
#       }
#       if (!is.finite(high)) {
#         high = NA
#       }
#       dfn = truncate(dfn, low, high)
#     }
#     return(dfn)
#   })
#   names(tmp) = colnames(posteriors_df)
#   priors_list = modifyList(priors_list, tmp)
#   return(priors_list)
# }
