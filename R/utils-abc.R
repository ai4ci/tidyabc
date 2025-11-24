## Common ----

#' Common internal parameters
#'
#' @param sim_df the output of a wave of simulation including a
#'   `abc_summary_distance` column
#' @param priors_list the list of priors as a named list of `dist_fn`s, plus one
#'   sided formulae as constraints and 2 sided formulae as derived values.
#'   A correlation matrix may also be present as an attribute `cor`.
#' @param epsilon epsilon is a tolerance threshold that controls how closely
#'   simulated summaries must match the observed ones to be considered
#'   plausible. This is in the unit of `abc_summary_distance`. Initially the 0.5
#'   quantile of distances, in subsequent waves this might be decreased.
#'   It is the scale parameter of the kernel function. $K_h(|u|)$
#' @param prev_sim_df the output of a previous ABC wave including a
#'   `abc_weight` column
#' @param proposal_list a list of empirical probability distributions that map
#'   MVN space to proposal space, and are the "prior" for each adaptive wave.
#'   This is already used to generate the proposals and their mapping in `sim_df`
#'
#' @name common_internal
#' @keywords internal
NULL


#' Combine simulation scores to a single distance metric.
#'
#'
#' This function calculates the overall score for an individual simulation,
#' based on combining the components output by `score_list` comparing them to
#' the target component scores of the original observed data. The Mahalanobis
#' distance is calculated relative to the distribution of the first wave (i.e.
#' when the priors were directly sampled).
#'
#' @param component_scores - a list column of scores. each entry is itself a list
#' @inheritParams tidyabc_common
#'
#' @returns a vector of distances
#'
#' @keywords internal
.summary_distance = function(
  component_scores,
  obsscores = NULL,
  distance_method = c("euclidean", "manhattan", "mahalanobis"),
  wave1_metrics = NULL,
  scoreweights = NULL
) {
  method = match.arg(distance_method)
  result = rep(NA, length(component_scores))

  nonnulls = which(
    !sapply(component_scores, is.null) &
      sapply(component_scores, function(l) all(sapply(l, is.finite)))
  )

  metrics = posterior_distance_metrics(
    component_scores,
    obsscores,
    keep_data = TRUE
  )

  # simulation scores as matrix:
  simscores = metrics$simscores

  # reference scores
  obsscores = metrics$obsscores
  obsscores = matrix(
    unlist(obsscores),
    nrow = nrow(simscores),
    ncol = ncol(simscores),
    byrow = TRUE
  )

  # score weights
  if (is.null(scoreweights)) {
    scoreweights = rep(1, ncol(simscores))
  }
  if (is.null(names(scoreweights))) {
    names(scoreweights) = colnames(simscores)
  }
  if (!identical(sort(names(scoreweights)), sort(colnames(simscores)))) {
    stop(
      "`scoreweights` names do not match `scorer_fn` outputs. Should be a vector with names: ",
      paste0(colnames(simscores), collapse = ",")
    )
  }
  scoreweights = unlist(scoreweights)
  scoreweights = scoreweights[colnames(simscores)]

  if (is.null(wave1_metrics)) {
    wave1_scoreweights = metrics$scoreweights
    wave1_cov = metrics$cov

    if (method == "normalised") {
      scoreweights = scoreweights * wave1_scoreweights
    }
  }

  delta = metrics$deltascores
  scoreweights = matrix(
    unlist(scoreweights),
    nrow = nrow(simscores),
    ncol = ncol(simscores),
    byrow = TRUE
  )

  if (method == "euclidean" || method == "normalised") {
    delta = delta * scoreweights
    result[nonnulls] = sqrt(rowSums(delta^2))
  }

  if (method == "manhattan") {
    delta = delta * scoreweights
    result[nonnulls] = rowSums(abs(delta))
  }

  if (method == "mahalanobis") {
    # I want columns that are high weighted to appear further away from the
    result[nonnulls] = suppressWarnings(stats::mahalanobis(
      simscores * scoreweights,
      obsscores * scoreweights,
      wave1_cov
    ))
  }

  return(result)
}

#' Covariance from component scores
#'
#' Supports calculation of Mahalanobis distance.
#'
#' @param component_scores a list of lists of scores, one per simulation
#'
#' @returns a covariance matrix
#' @keywords internal
.distance_covariance = function(
  component_scores
) {
  nonnulls = which(
    !sapply(component_scores, is.null) &
      sapply(component_scores, function(l) all(sapply(l, is.finite)))
  )

  # unnest scores from list of lists into a matrix:
  # this relies on the naming order being consistent, which it should be as output
  # from purrr...
  simscores = matrix(
    unname(unlist(component_scores[nonnulls])),
    nrow = length(nonnulls),
    byrow = TRUE
  )
  wave1_cov = stats::cov(simscores)
  return(wave1_cov)
}

# x = list(list(a=1,b=2,c=3),list(a=1,b=2,c=3),list(a=1,b=2,c=3))

#' Generate comparison metrics for two sequential waves
#'
#' @inheritParams common_internal
#'
#' @returns a nested tibble with 2 columns `summary` and `per_parameter`
#'   with stats in each. The `summary` stats are
#' @keywords internal
.compare_waves = function(sim_df, prev_sim_df = NULL, priors_list, wave) {
  nms = priors_list@params

  sim_weight = suppressWarnings(sim_df$abc_weight)

  if (is.null(sim_weight)) {
    # This really shouldn't happen...
    stop("sim_df has no weight?")
    sim_weight = rep(1 / nrow(sim_df), nrow(sim_df))
  }

  if (is.null(prev_sim_df)) {
    return(dplyr::tibble(
      summary = list(dplyr::tibble(
        abs_distance_redn = NA,
        rel_distance_redn = NA,
        ESS = 1 / sum(sim_weight^2)
      )),
      per_param = list(dplyr::tibble(
        param = nms,
        IQR_95_redn = NA,
        abs_variance_redn = NA,
        rel_mean_change = NA
      ))
    ))
  }

  # In wave 1 the prev_sim_df is the pre-simulated version of the
  # priors, and will have no distances or weights.
  # The distances we cannot make up but the weights we can assume to be uniform
  if (!"abc_weight" %in% colnames(prev_sim_df)) {
    prev_sim_df$abc_weight = rep(1 / nrow(prev_sim_df), nrow(prev_sim_df))
  }

  if (!"abc_summary_distance" %in% colnames(prev_sim_df)) {
    prev_sim_df$abc_summary_distance = rep(NA, nrow(prev_sim_df))
  }

  prev_weight = prev_sim_df$abc_weight

  cov_new = sim_df %>%
    dplyr::select(dplyr::all_of(nms)) %>%
    as.matrix() %>%
    stats::cov.wt(wt = sim_weight)

  cov_old = prev_sim_df %>%
    dplyr::select(dplyr::all_of(nms)) %>%
    as.matrix() %>%
    stats::cov.wt(wt = prev_weight)

  epsilon_new = stats::quantile(sim_df$abc_summary_distance, 0.5, na.rm = TRUE)
  epsilon_old = stats::quantile(
    prev_sim_df$abc_summary_distance,
    0.5,
    na.rm = TRUE
  )

  # The difference between the 95% IQR pre and post update. If this is low then
  # things are not improving.
  quantile_range_redn = lapply(nms, function(nm) {
    if (all(is.na(prev_weight))) {
      return(NA)
    }
    unname(
      diff(wquantile(
        p = c(0.025, 0.975),
        x = prev_sim_df[[nm]],
        w = prev_weight,
        link = priors_list[[nm]]
      )) -
        diff(wquantile(
          p = c(0.025, 0.975),
          x = sim_df[[nm]],
          w = sim_weight,
          link = priors_list[[nm]]
        ))
    )
  })
  names(quantile_range_redn) = nms

  return(dplyr::tibble(
    wave = wave,
    summary = list(dplyr::tibble(
      abs_distance = epsilon_new,
      abs_distance_redn = (epsilon_old - epsilon_new),
      rel_distance_redn = (epsilon_old - epsilon_new) / epsilon_old,
      ESS = 1 / sum(sim_weight^2)
    )),
    per_param = list(dplyr::tibble(
      param = names(cov_new$center),
      IQR_95_redn = unlist(quantile_range_redn[names(cov_new$center)]),
      abs_variance_redn = (diag(cov_old$cov) - diag(cov_new$cov)),
      rel_mean_change = abs((cov_new$center - cov_old$center) / cov_old$center)
    ))
  ))
}


# Crate score function with observed data for transfer to different threads
# crate includes .p progressr option.
.crate_scorer_fn = function(scorer_fn, obsdata, debug = FALSE) {
  if (is.function(scorer_fn)) {
    if (!identical(names(formals(scorer_fn))[1:2], c("simdata", "obsdata"))) {
      stop(
        "`scorer_fn` must have 2 parameters, named `simdata` (or `.x`) then `obsdata` (or `.y`)"
      )
    }
  } else {
    scorer_fn = rlang::as_function(scorer_fn)
  }
  # Crate function for parallelisation.
  scorer_crate = carrier::crate(
    function(simdata, .p = NULL) {
      if (is.null(simdata)) {
        return(NULL)
      }
      if (!is.null(.p)) {
        .p()
      }
      tryCatch(
        scorer_fn(simdata = simdata, obsdata = obsdata),
        error = function(e) {
          warning(e$message)
          if (!!debug) {
            debug(scorer_fn)
            on.exit(undebug(scorer_fn), add = TRUE)
            scorer_fn(simdata = simdata, obsdata = obsdata)
          }
          return(NULL)
        }
      )
    },
    scorer_fn = .autocrate_fn(scorer_fn),
    obsdata = obsdata
  )
  return(scorer_crate)
}

# Crate sim function with .p progressr option
.crate_sim_fn = function(sim_fn, debug = FALSE) {
  # Crate function for parallelisation.

  sim_crate = carrier::crate(
    function(..., .p = NULL) {
      if (!is.null(.p)) {
        .p()
      }
      e = NULL
      args = rlang::list2(...)
      args = args[names(args) %in% names(formals(sim_fn))]
      tryCatch(
        {
          do.call(sim_fn, args)
        },
        error = function(e) {
          warning(e$message)
          if (!!debug) {
            debug(sim_fn)
            on.exit(undebug(sim_fn), add = TRUE)
            do.call(sim_fn, args)
          }
          return(NULL)
        }
      )
    },
    sim_fn = .autocrate_fn(sim_fn)
  )
}


## Internal workflow functions ----

# Do simulation and scoring with settings for SMC/Adaptive
# and not for resampling.
.abc_do_one = function(
  obsdata,
  sim_df,
  sim_fn,
  scorer_fn,
  obsscores = NULL,
  distance_method = "euclidean",
  keep_simulations = FALSE,
  seed = NULL,
  parallel = FALSE,
  wave1_metrics = NULL,
  n_resamples = 1,
  debug = FALSE,
  scoreweights = NULL,
  wave = 0
) {
  if (!is.null(seed)) {
    seed = set.seed(seed)
    on.exit(set.seed(seed), add = TRUE)
  }

  sim_crate = .crate_sim_fn(sim_fn, debug = debug)

  sim_df = .abc_do_simulation_and_scoring(
    sim_df = sim_df,
    sim_crate = sim_crate,
    n_resamples = 1,
    keep_simulations = keep_simulations,
    scorer_crate = .crate_scorer_fn(scorer_fn, obsdata, debug = debug),
    obsscores = obsscores,
    distance_method = distance_method,
    parallel = parallel,
    wave1_metrics = wave1_metrics,
    scoreweights = scoreweights,
    wave = wave
  )

  return(sim_df)
}


# Handles running a simulation in parallel and scoring it
# May keep the simulation output or not depending on options
.abc_do_simulation_and_scoring = function(
  sim_df,
  sim_crate,
  n_resamples,
  keep_simulations,
  scorer_crate = NULL,
  obsscores = NULL,
  distance_method = "euclidean",
  parallel = FALSE,
  wave1_metrics = NULL,
  scoreweights = NULL,
  wave = 0
) {
  weights = suppressWarnings(sim_df$abc_weight)

  if (!is.null(weights)) {
    # weighted resampling
    sim_df = sim_df %>%
      dplyr::slice_sample(
        n = nrow(sim_df) * n_resamples,
        weight_by = weights,
        replace = TRUE
      )
  } else {
    # no weights - just repeat the data n_resamples times
    sim_df = dplyr::bind_rows(
      rep(list(sim_df), n_resamples)
    )
  }

  # TODO: Switch to progressr progress bars
  # Issue URL: https://github.com/ai4ci/tidyabc/issues/12
  # I think this would not support wave ID for progressbar?
  # .wave = sprintf("Wave %d simulations:", wave)
  # p = progressr::progressor(steps = nrow(sim_df))

  .wave = interactive() &&
    !is.null(getOption("knitr.in.progress")) &&
    !identical(Sys.getenv("IN_PKGDOWN"), "true")

  if (parallel) {
    fn_pmap = function(.l, .f, ..., .wave = TRUE) {
      furrr::future_pmap(
        .l,
        .f,
        ...,
        .options = furrr::furrr_options(seed = TRUE),
        .progress = .wave
      )
    }
    fn_map = function(.x, .f, ..., .wave = TRUE) {
      furrr::future_map(
        .x,
        .f,
        ...,
        .options = furrr::furrr_options(seed = TRUE),
        .progress = .wave
      )
    }
  } else {
    fn_pmap = function(.l, .f, ..., .wave = TRUE) {
      # had to switch a standalone purrr for performance reasons
      pmap(.l, .f, ..., .progress = .wave)
    }
    fn_map = function(.x, .f, ..., .wave = TRUE) {
      map(.x, .f, ..., .progress = .wave)
    }
  }

  p = NULL

  if (is.null(scorer_crate)) {
    # Score fn not present so we keep simulations.
    sim_df = sim_df %>%
      dplyr::mutate(abc_sim = fn_pmap(., sim_crate, .p = p, .wave = .wave))
  } else {
    if (keep_simulations) {
      # keep both sim_df and scores
      sim_df = sim_df %>%
        dplyr::mutate(
          abc_sim = fn_pmap(., sim_crate, .p = p, .wave = .wave)
        ) %>%
        dplyr::mutate(
          abc_component_score = fn_map(
            abc_sim,
            scorer_crate
            # no progress needed
          )
        )
    } else {
      # just calculate scores - e.g. from within the .abc_do_one
      sim_df = sim_df %>%
        dplyr::mutate(
          abc_component_score = fn_pmap(
            .,
            function(...) {
              sim_df = sim_crate(...)
              return(scorer_crate(sim_df, .p = p))
            },
            .wave = .wave
          )
        )
    }

    # given the components calculate the overall distance.
    sim_df = sim_df %>%
      dplyr::mutate(
        abc_summary_distance = .summary_distance(
          abc_component_score,
          obsscores = obsscores,
          distance_method = distance_method,
          wave1_metrics = wave1_metrics,
          scoreweights = scoreweights
        )
      )
  }
  return(sim_df)
}

## Adaptive ----

#' Calculate weights for particles in a new wave
#'
#' The ABC weights need to be calculated for sampling from the proposal
#' distribution. They depend on the priors, an acceptance tolerance and the
#' proposal probability distribution.
#'
#' @inheritParams common_internal
#' @inheritParams tidyabc_common
#' @inheritParams .log_kernel
#'
#' @returns the `sim_df` with an `abc_weight` column
#' @keywords internal
.calculate_weights_adaptive = function(
  sim_df,
  priors_list,
  acceptance_rate,
  proposal_list,
  kernel,
  use_proposal_correlation,
  ess_limit,
  max_recover
) {
  # N.B. this is not currently being used
  # as tends to produce worse results than a kernel distance only weight.

  params = names(priors_list)
  params = params[params != ""]

  if (identical(priors_list, proposal_list)) {
    # In the first wave the effect of priors and proposal distrbituion
    # cancels out:
    log_q = rep(0, nrow(sim_df))
    log_prior = rep(0, nrow(sim_df))
  } else {
    # In the adaptive approach the MVN space for every wave is always centred at 0
    # regardless of the wave, however correlation structure will change after each
    # one. After each wave the mapping between MVN and proposal space (copula) is
    # held in the proposal_list distribution functions It is these mappings from
    # MVN -> proposal space that get updated at each round.

    # From proposal mappings
    cor = NULL
    if (!is.null(proposal_list)) {
      cor = attr(proposal_list, "cor")
    }
    if (is.null(cor) || !use_proposal_correlation) {
      cor = diag(length(params))
    }

    # We have to map back from real to proposal MVN space as particles will
    # be a mix of particles from multiple waves.
    theta_new = sapply(params, function(nm) {
      prior = proposal_list[[nm]]
      # This is the current set of proposals in proposal space:
      theta_star = sim_df[[nm]]
      # need to map this back to MVN space using prior copula
      # rather than proposal
      # convert target to uniform (prior$p) & map to MVN (qnorm):
      stats::qnorm(prior$p(theta_star))
    })

    # In adaptive we cannot do this as meaning of MVN space changes over waves:
    # theta_new = sim_df %>%
    #   dplyr::select(dplyr::starts_with("abc_mvn_")) %>%
    #   as.matrix()

    # In proposal MVN space:
    # This should be the same copula that was used to generate the proposals:
    # Correlation is accounted for unless switched off
    log_q = mvtnorm::dmvnorm(theta_new, sigma = cor, log = TRUE)
    # Constant terms:
    log_q = log_q - max(log_q, na.rm = TRUE)

    # The prior probability of this particle is defined in the MVN space using the
    # prior mappings. Everything is independent in the prior so no correlation.

    theta_prior = sapply(params, function(nm) {
      prior = priors_list[[nm]]
      # This is the current set of proposals in proposal space:
      theta_star = sim_df[[nm]]
      # need to map this back to MVN space using prior copula
      # rather than proposal
      # convert target to uniform (prior$p) & map to MVN (qnorm):
      stats::qnorm(prior$p(theta_star))
    })

    # In prior MVN space:
    log_prior = mvtnorm::dmvnorm(theta_prior, log = TRUE)
    # Constant terms:
    log_prior = log_prior - max(log_prior, na.rm = TRUE)
    # Fix dmvnorn NaNs if not finite inputs:
    # -Inf because log(P=0)
    log_prior[!apply(is.finite(theta_prior), MARGIN = 1, all)] = -Inf

    # P_prior(theta) / P_proposal(theta)
    # term can dominate, and create proposals that are outside of the
    # gradually decreasing kernel radius, leading to low ESS and
    # lack of convergence. Its the log_q term that causes the issue here but
    # adjustment would be non linear in normal space.
    # The rationale here is we want to prioritise surprising results that are
    # close to the result we want. This is rationale for the expit in the
    # conversion
  }

  distances = sim_df$abc_summary_distance

  ess = 0
  recover = 0
  tmp_sim_df = sim_df

  while (recover < max_recover) {
    recover = recover + 1
    epsilon = stats::quantile(
      distances,
      probs = acceptance_rate
    )

    # calculate a weight based on kernel:
    log_abc_kernel = .log_kernel(distances, epsilon, kernel) -
      .log_kernel(0, epsilon = epsilon, kernel)

    log_weight = log_abc_kernel + log_prior - log_q

    # convert log weight to importance
    # log_weight = pmin(log_weight,0)
    # abc_weight = exp(log_weight - max(log_weight))
    # Convert log_weight to importance as probability:
    # abc_weight = .expit(log_weight)
    tmp_sim_df = sim_df %>%
      dplyr::mutate(abc_weight = .expit(log_weight)) %>%
      dplyr::filter(abc_weight > sqrt(.Machine$double.eps)) %>%
      dplyr::mutate(abc_weight = abc_weight / sum(abc_weight))

    ess = 1 / sum(tmp_sim_df$abc_weight^2)

    # browser()

    if (ess < ess_limit[1] && acceptance_rate < 0.99) {
      # increase acceptance rate if ESS is low and acceptance not already at limit
      acceptance_rate = scale_probability(acceptance_rate, 1.5)
    } else if (ess > ess_limit[2] && acceptance_rate > 0.05) {
      # decrease acceptance rate if ESS is high up to small 5% acceptance limit
      # This would be 5 retries if it started at 50%
      acceptance_rate = scale_probability(acceptance_rate, 1 / 1.5)
    } else {
      break
    }
  }

  return(list(
    acceptance_rate = acceptance_rate,
    sim_df = tmp_sim_df,
    ess = ess
  ))
}


## Basic ----

#' Calculate weights for particles in a new wave
#'
#' The ABC weights need to be calculated for sampling from the proposal
#' distribution. In the situation where we are sampling from the prior.
#' this is just their distance by a kernel function.
#'
#' @inheritParams common_internal
#' @inheritParams tidyabc_common
#' @inheritParams .log_kernel
#'
#' @returns the `sim_df` with an `abc_weight` column
#' @keywords internal
.calculate_weights_wave_one = function(
  sim_df,
  epsilon,
  kernel
) {
  distances = sim_df$abc_summary_distance
  log_weight = .log_kernel(distances, epsilon, kernel) -
    .log_kernel(0, epsilon, kernel)

  sim_df %>%
    dplyr::mutate(
      # Compute unnormalized ABC weights (Gaussian kernel)
      # abc_weight = exp(log_weight - max(log_weight))
      abc_weight = .expit(log_weight)
    ) %>%
    dplyr::mutate(
      abc_weight = abc_weight / sum(abc_weight)
    )
}

## SMC ----

#' Calculate weights for particles in a new wave
#'
#' The ABC weights need to be calculated for sampling from the proposal
#' distribution. They depend on the priors, an acceptance tolerance and the
#' weighted particles from the previous wave, which define the proposal
#' probability distribution.
#'
#' @inheritParams common_internal
#' @inheritParams tidyabc_common
#' @inheritParams .log_kernel
#'
#' @returns the `sim_df` with an `abc_weight` column
#' @keywords internal
.calculate_weights_smc = function(
  sim_df,
  # priors_list,
  epsilon,
  prev_sim_df = NULL,
  kernel
) {
  # priors = names(priors_list)
  # priors = priors[priors != ""]

  theta_new = sim_df %>%
    dplyr::select(dplyr::starts_with("abc_mvn_")) %>%
    as.matrix()

  distances = sim_df$abc_summary_distance

  # log_prior = sapply(priors, function(nm) {
  #   prior = priors_list[[nm]]
  #   col = sim_df[[nm]]
  #   return(prior$d(col, log = TRUE))
  # })

  # The prior probability of this particle is defined in the MVN space.
  # everything is assumed independent at this stage
  # log_prior = dnorm(theta_new, log = TRUE)
  # log_prior = apply(log_prior, MARGIN = 1, sum)

  log_prior = mvtnorm::dmvnorm(theta_new, log = TRUE)

  log_abc_kernel = .log_kernel(distances, epsilon, kernel)
  # log_abc_kernel = log(as.integer(distances < epsilon))

  if (!is.null(prev_sim_df)) {
    theta_prev = prev_sim_df %>%
      dplyr::select(dplyr::starts_with("abc_mvn_")) %>%
      as.matrix()
    w_prev = suppressWarnings(prev_sim_df$abc_weight)
    Sigma_score = .posterior_covariance(theta_prev, w_prev)
    log_q = .log_q_proposal_vectorized(
      theta_new,
      theta_prev,
      w_prev,
      Sigma_score
    )
  } else {
    log_q = log_prior
  }

  # Constant terms:
  log_M = .log_kernel(0, epsilon = epsilon, kernel) +
    max(log_prior) -
    max(log_q)
  log_weight = log_prior + log_abc_kernel - log_q - log_M

  sim_df %>%
    dplyr::mutate(
      # convert log weight to importance
      # log_weight = pmin(log_weight,0)
      # abc_weight = exp(log_weight - max(log_weight))
      # Convert log_weight to importance as probability:
      abc_weight = .expit(log_weight)
    ) %>%
    dplyr::mutate(
      abc_weight = abc_weight / sum(abc_weight)
    )
}

#' Generate a new set of particles from a previous wave using perturbation
#'
#' This set of particles will be close to existing ones, depending on the
#' kernel_t parameter.
#'
#' @inheritParams common_internal
#' @inheritParams tidyabc_common
#'
#' @returns a `sim_df`
#' @keywords internal
.sample_using_perturbation = function(
  priors_list,
  n_sims,
  prev_sim_df
) {
  # kernel_t
  # kernel_t A kernel bandwidth parameter for proposals. This controls the
  #   amount of noise that particles are perturbed by (and hence the spread of
  #   the PDF of the proposal distribution), and 1 is approximately 34% of the
  #   proposal distribution at the centre. Smaller values (default is 0.2) give a
  #   smaller step size in generating new proposals, and proposals will be closer
  #   to currently accepted particles.

  theta_prev = prev_sim_df %>%
    dplyr::select(dplyr::starts_with("abc_mvn_")) %>%
    as.matrix()

  w_prev = suppressWarnings(prev_sim_df$abc_weight)
  if (is.null(w_prev)) {
    stop("Weights must be calculated first.")
  }
  Sigma_pert = .posterior_covariance(theta_prev, w_prev) #, kernel_t)

  # Step 1: resample indices with probability w_prev
  idx = sample(
    seq_len(nrow(theta_prev)),
    size = n_sims,
    replace = TRUE,
    prob = w_prev
  )

  # Step 2: perturb each selected particle in MVN space
  z = theta_prev[idx, , drop = FALSE] +
    mvtnorm::rmvnorm(
      n_sims,
      mean = rep(0, ncol(theta_prev)),
      sigma = Sigma_pert
    )

  sim_df = z %>% dplyr::as_tibble()

  nms = priors_list@params

  for (i in seq_along(nms)) {
    nm = nms[[i]]
    mvn_nm = sprintf("abc_mvn_%s", nm)
    prior_dist = priors_list[[nm]]
    # Map to proposal space
    sim_df[[nm]] = prior_dist$q(stats::pnorm(z[, mvn_nm]))
  }

  return(sim_df)
}


#' Covariance matrix from the posteriors
#'
#' @param theta a set of particles as a matrix
#' @param weights the weights of the previous particles assumed calculated using
#'   the same kernel
#' @param kernel_t A kernel bandwidth parameter for proposals. This controls the
#'   amount of noise that particles are perturbed by (and hence the spread of
#'   the PDF of the proposal distribution), and 1 is approximately 34% of the
#'   proposal distribution at the centre. Smaller values (default is 0.2) give a
#'   smaller step size in generating new proposals, and proposals will be closer
#'   to currently accepted particles.
#' @inheritParams tidyabc_common
#'
#' @returns a covariance matrix, optionally scaled
#' @keywords internal
.posterior_covariance = function(theta, weights, kernel_t = 1) {
  if (is.null(weights)) {
    weights = rep(1 / nrow(theta), nrow(theta))
  }

  # Weighted covariance
  cov_w = stats::cov.wt(theta, wt = weights, method = "ML")$cov

  # Scale it (standard trick for random walk)
  d = ncol(theta) # dimension = 4
  cov_w = kernel_t^2 / d * cov_w

  # Ensure positive definite (add jitter if needed)
  Sigma_pert = as.matrix(Matrix::nearPD(cov_w, keepDiag = TRUE)$mat)
  return(Sigma_pert)
}


#' Calculate the probability of new proposals based on previous
#'
#' @param theta_new a set of new particles as a matrix
#' @param theta_prev a set of previous particles as a matrix
#' @param w_prev the weights of the previous particles assumed calculated using
#'   the same kernel
#' @param Sigma_pert the covariance of the kernel.
#'
#' @returns a vector of probabilities for each row of theta_new
#' @keywords internal
#'
#' @unit
#' # Previous particles: 100 samples from N(0, 1)
#' N_prev <- 100
#' theta_prev <- matrix(rnorm(N_prev, 0, 1), ncol = 1)
#' w_prev <- rep(1/N_prev, N_prev)  # uniform weights
#' Sigma_pert <- matrix(0.5^2, nrow = 1, ncol = 1)  # scalar variance
#'
#' # New point
#' theta_new <- matrix(c(0.5), ncol = 1)
#'
#' # Compute via vectorized function
#' q_val_fast <- .log_q_proposal_vectorized(theta_new, theta_prev, w_prev, Sigma_pert)
#'
#' # Compute via slow loop (ground truth)
#' q_val_slow <- 0
#' for (j in 1:N_prev) {
#'   q_val_slow <- q_val_slow + w_prev[j] * dnorm(theta_new, mean = theta_prev[j], sd = sqrt(Sigma_pert))
#' }
#' expect_equal(exp(q_val_fast), as.vector(q_val_slow))
.log_q_proposal_vectorized <- function(
  theta_new,
  theta_prev,
  w_prev,
  Sigma_pert
) {
  # theta_new: M x d
  # theta_prev: N x d
  # w_prev: N-vector (sums to 1)
  # Returns: M-vector of proposal densities

  M <- nrow(theta_new)
  N <- nrow(theta_prev)
  d <- ncol(theta_new)

  # 1. Whiten the space: solve L %*% t(L) = Sigma_pert
  L <- chol(Sigma_pert) # Sigma = L'L, so L is upper triangular
  # Solve L' %*% X = theta  => X = solve(t(L), theta)
  theta_new_white <- t(solve(t(L), t(theta_new))) # M x d
  theta_prev_white <- t(solve(t(L), t(theta_prev))) # N x d

  # 2. Compute pairwise squared Euclidean distances: M x N matrix
  # Use proxy::dist with method = "Euclidean", but we need cross-distance
  # Manual: ||a_i - b_j||^2 = ||a_i||^2 + ||b_j||^2 - 2 a_i b_j^T
  a2 <- rowSums(theta_new_white^2) # M-vector
  b2 <- rowSums(theta_prev_white^2) # N-vector
  ab <- tcrossprod(theta_new_white, theta_prev_white) # M x N

  D2 <- outer(a2, b2, "+") - 2 * ab # M x N matrix of squared Mahalanobis distances

  # 3. Compute log-density constant
  log_det <- determinant(Sigma_pert, logarithm = TRUE)$modulus
  log_const <- -0.5 * (d * log(2 * pi) + log_det)

  # 4. Compute densities: M x N matrix
  log_dens_matrix <- log_const - 0.5 * D2

  # 5. Weighted sum over N (columns) for each M (rows)
  log_w_prev <- log(w_prev + .Machine$double.eps)

  # For each row i: log(sum_j exp(log_dens[i,j] + log_w_prev[j]))
  log_q_vals <- apply(
    log_dens_matrix + matrix(log_w_prev, nrow = M, ncol = N, byrow = TRUE),
    1,
    function(row) {
      max_val <- max(row)
      max_val + log(sum(exp(row - max_val)))
    }
  )

  as.vector(log_q_vals)
}
