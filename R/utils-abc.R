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
#'   quantile of distances, in subsequent waves this might be decreased
#' @param prev_sim_df the output of a previous ABC wave including a
#'   `abc_weight` column
#' @param proposal_list a list of empirical probability distributions that map
#'   MVN space to proposal space, and are the "prior" for each adaptive wave.
#'   This is already used to generate the proposals and their mapping in `sim_df`
#' @name common_internal
#' @keywords internal
NULL

#' Sample for the prior distribution
#'
#' This uses a multivariate normal and a copula to generate correlated
#' structure. The correlation is held as an attribute in the proposals.
#' N.B. this used to construct derived values and apply constraints but this
#' was complex and has been deferred. Constraints change the probability
#' distribution also and will affect validity of proposal distribution.
#'
#' @inheritParams common_internal
#' @inheritParams tidyabc_common
#'
#' @returns a data frame of samples in MVN (prefixed `abc_mvn_`) and parameter
#' space.
#' @keywords internal
.sample_priors = function(proposal_list, n_sims) {
  sim_df = dplyr::tibble(.rows = n_sims)

  nms = names(proposal_list)
  nms = nms[nms != ""]

  cor = attr(proposal_list, "cor")
  if (is.null(cor)) {
    cor = diag(length(nms))
  }

  # generate the correlated measures as normals:
  z = mvtnorm::rmvnorm(n_sims, mean = rep(0, length(nms)), sigma = cor)

  # convert to uniform (pnorm) & map to targets (prior$q:
  for (i in seq_along(nms)) {
    nm = nms[[i]]
    mvn_nm = sprintf("abc_mvn_%s", nm)
    prior_dist = proposal_list[[nm]]
    if (!is.dist_fns(prior_dist)) {
      stop("Prior `", nm, "` is not an object of type `dist_fns`.")
    }
    sim_df[[mvn_nm]] = z[, i]
    sim_df[[nm]] = prior_dist$q(stats::pnorm(z[, i]))
  }

  # derived = unname(proposal_list[nms == ""])
  # 2 sided formulae are derived / 1 sided constraints
  # tmp = sapply(sapply(derived, rlang::f_lhs), is.null)
  # constraints = derived[as.logical(tmp)]
  # derived = derived[!as.logical(tmp)]

  # derived values:
  # for (form in derived) {
  #   tgt = rlang::f_lhs(form)
  #   expr = rlang::f_rhs(form)
  #   sim_df = sim %>% dplyr::mutate(!!tgt := !!expr)
  # }

  # sim_df = .apply_constraints_and_derived(sim_df, proposal_list)

  return(sim_df)
}


#' Generate comparison metrics for two sequential waves
#'
#' @inheritParams common_internal
#'
#' @returns a nested tibble with 2 columns `summary` and `per_parameter`
#'   with stats in each. The `summary` stats are
#' @keywords internal
.compare_waves = function(sim_df, prev_sim_df = NULL, priors_list) {
  nms = names(priors_list)
  nms = nms[nms != ""]

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
      Hmisc::wtd.quantile(
        x = prev_sim_df[[nm]],
        weights = prev_weight,
        probs = 0.975,
        normwt = TRUE,
        na.rm = TRUE
      ) -
        Hmisc::wtd.quantile(
          x = prev_sim_df[[nm]],
          weights = prev_weight,
          probs = 0.025,
          normwt = TRUE,
          na.rm = TRUE
        ) -
        Hmisc::wtd.quantile(
          x = sim_df[[nm]],
          weights = sim_weight,
          probs = 0.975,
          normwt = TRUE,
          na.rm = TRUE
        ) +
        Hmisc::wtd.quantile(
          x = sim_df[[nm]],
          weights = sim_weight,
          probs = 0.025,
          normwt = TRUE,
          na.rm = TRUE
        )
    )
  })
  names(quantile_range_redn) = nms

  return(dplyr::tibble(
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
.crate_scorer_fn = function(scorer_fn, obsdata) {
  if (is.function(scorer_fn)) {
    if (!all(names(formals(scorer_fn))[1:2] != c("simdata", "obsdata"))) {
      stop(
        "`scorer_fn` must have 2 parameters, named `simdata` (or `.x`) then `obsdata` (or `.y`)"
      )
    }
  } else {
    scorer_fn = rlang::as_function(scorer_fn)
  }

  require(stats)
  require(utils)
  # Crate function for parallelisation.
  scorer_crate = .autocrate(
    function(simdata, .p = NULL) {
      if (!is.null(.p)) {
        .p()
      }
      scorer_fn(simdata, obsdata)
    },
    scorer_fn = scorer_fn,
    obsdata = obsdata
  )
}

# Crate sim function with .p progressr option
.crate_sim_fn = function(sim_fn) {
  # Crate function for parallelisation.
  require(stats)
  require(utils)
  sim_crate = .autocrate(
    function(..., .p = NULL) {
      if (!is.null(.p)) {
        .p()
      }
      args = rlang::list2(...)
      args = args[names(args) %in% names(formals(sim_fn))]
      do.call(sim_fn, args)
    },
    sim_fn = sim_fn
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
  wave1_cov = NULL,
  n_resamples = 1
) {
  if (!is.null(seed)) {
    seed = set.seed(seed)
    on.exit(set.seed(seed), add = TRUE)
  }

  sim_crate = .crate_sim_fn(sim_fn)

  sim_df = .abc_do_simulation_and_scoring(
    sim_df = sim_df,
    sim_crate = sim_crate,
    n_resamples = 1,
    keep_simulations = keep_simulations,
    scorer_crate = .crate_scorer_fn(scorer_fn, obsdata),
    obsscores = obsscores,
    distance_method = distance_method,
    parallel = parallel,
    wave1_cov = wave1_cov
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
  wave1_cov = NULL
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

  #TODO: wave ID for progressbar?

  if (parallel) {
    fn_pmap = function(...) furrr::future_pmap(..., .progress = TRUE)
    fn_map = function(...) furrr::future_map(..., .progress = TRUE)
  } else {
    fn_pmap = function(...) purrr::pmap(..., .progress = TRUE)
    fn_map = function(...) purrr::map(..., .progress = TRUE)
  }
  # TODO:
  # p = progressr::progressor(steps = nrow(sim_df))
  p = NULL

  if (is.null(scorer_crate)) {
    # Score fn not present so we keep simulations.
    sim_df = sim_df %>%
      dplyr::mutate(abc_sim = fn_pmap(., sim_crate, .p = p))
  } else {
    if (keep_simulations) {
      # keep both sim_df and scores
      sim_df = sim_df %>%
        dplyr::mutate(
          abc_sim = fn_pmap(., sim_crate, .p = p)
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
            }
          )
        )
    }

    # given the components calculate the overall distance.
    sim_df = sim_df %>%
      dplyr::mutate(
        abc_summary_distance = summary_distance(
          abc_component_score,
          obsscores = obsscores,
          distance_method = distance_method,
          wave1_cov = wave1_cov
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
#'
#' @returns the `sim_df` with an `abc_weight` column
#' @keywords internal
.calculate_weights_adaptive = function(
  sim_df,
  priors_list,
  epsilon,
  proposal_list
) {
  params = names(priors_list)
  params = params[params != ""]

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
  if (is.null(cor)) {
    cor = diag(length(params))
  }

  theta_new = sim_df %>%
    dplyr::select(dplyr::starts_with("abc_mvn_")) %>%
    as.matrix()

  # In proposal MVN space:
  log_q = mvtnorm::dmvnorm(theta_new, sigma = cor, log = TRUE)

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
  # Fix dmvnorn NaNs if not finite inputs:
  # -Inf because log(P=0)
  log_prior[!apply(is.finite(theta_prior), MARGIN = 1, all)] = -Inf

  distances = sim_df$abc_summary_distance
  log_abc_kernel = -0.5 * (distances / epsilon)^2

  log_weight = log_prior + log_abc_kernel - log_q

  if (any(is.na(log_weight))) {
    browser()
  }

  sim_df %>%
    dplyr::mutate(
      # Compute unnormalized ABC weights (Gaussian kernel)
      abc_weight = exp(log_weight - max(log_weight))
    ) %>%
    dplyr::mutate(
      abc_weight = abc_weight / sum(abc_weight)
    )
}

# # Filter proposals to match constraints and boundaries imposed by priors
# #
# # @inheritParams common_internal
# # @returns the `sim_df` with invalid combinations removed, and replaced with
# # copies of valid proposals randomly.
# # @keywords internal
# .apply_constraints_and_derived = function(sim_df, priors_list) {
#   n = nrow(sim_df)
#   # derived
#   nms = names(priors_list)
#   param_nms = nms[nms != ""]
#   derived = unname(priors_list[nms == ""])
#   # 2 sided formulae are derived / 1 sided constraints
#   tmp = sapply(sapply(derived, rlang::f_lhs), is.null)
#   constraints = derived[as.logical(tmp)]
#   derived = derived[!as.logical(tmp)]
#
#   # derived values:
#   for (form in derived) {
#     tgt = rlang::f_lhs(form)
#     expr = rlang::f_rhs(form)
#     sim_df = sim %>% dplyr::mutate(!!tgt := !!expr)
#   }
#
#   max_it = getOption("tidyabc.max_oob", 1000)
#   while (max_it > 1) {
#     # check validity
#     sim_df = sim_df %>% dplyr::mutate(.valid = TRUE)
#     for (nm in param_nms) {
#       prior = priors_list[[nm]]
#       col = sim_df[[nm]]
#       keep = col >= prior$q(0) & col <= prior$q(1)
#       sim_df = sim_df %>% dplyr::mutate(.valid = .valid & keep)
#     }
#
#     for (constraint in constraints) {
#       expr = rlang::f_rhs(constraint)
#       sim_df = sim_df %>% dplyr::mutate(.valid = .valid & !!expr)
#     }
#
#     if (any(!sim_df$.valid)) {
#       # browser()
#       # fix invalid
#       valid = sim_df %>% dplyr::filter(.valid) %>% dplyr::select(-.valid)
#       invalid = sim_df %>% dplyr::filter(!.valid) %>% dplyr::select(-.valid)
#       alt = valid %>% dplyr::slice_sample(n = nrow(invalid), replace = TRUE)
#       # invalid replaced with midpoint of invalid and valid particles
#       rpl = dplyr::bind_rows(
#         invalid %>% dplyr::mutate(.id = row_number()),
#         alt %>% dplyr::mutate(.id = row_number())
#       ) %>%
#         dplyr::group_by(.id) %>%
#         dplyr::summarise(dplyr::across(dplyr::everything(), mean)) %>%
#         dplyr::ungroup() %>%
#         dplyr::select(-.id)
#
#       sim_df = dplyr::bind_rows(
#         valid,
#         rpl
#       )
#
#       max_it = max_it - 1
#     } else {
#       # exit loop and function
#       return(sim_df %>% dplyr::select(-.valid))
#     }
#   }
#   warning("Invalid proposals generated that could not be fixed.")
#   return(sim_df %>% dplyr::filter(.valid) %>% dplyr::select(-.valid))
# }

## Basic ----

#' Calculate weights for particles in a new wave
#'
#' The ABC weights need to be calculated for sampling from the proposal
#' distribution. In the situation where we are sampling from the prior.
#' this is just their distance by a kernel function.
#'
#' @inheritParams common_internal
#' @inheritParams tidyabc_common
#'
#' @returns the `sim_df` with an `abc_weight` column
#' @keywords internal
.calculate_weights_wave_one = function(
  sim_df,
  epsilon
) {
  distances = sim_df$abc_summary_distance
  log_weight = -0.5 * (distances / epsilon)^2

  sim_df %>%
    dplyr::mutate(
      # Compute unnormalized ABC weights (Gaussian kernel)
      abc_weight = exp(log_weight - max(log_weight))
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
#'
#' @returns the `sim_df` with an `abc_weight` column
#' @keywords internal
.calculate_weights_smc = function(
  sim_df,
  # priors_list,
  epsilon,
  prev_sim_df = NULL,
  kernel_t
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

  log_abc_kernel = -0.5 * (distances / epsilon)^2
  # log_abc_kernel = log(as.integer(distances < epsilon))

  if (!is.null(prev_sim_df)) {
    theta_prev = prev_sim_df %>%
      dplyr::select(dplyr::starts_with("abc_mvn_")) %>%
      as.matrix()
    w_prev = suppressWarnings(prev_sim_df$abc_weight)
    Sigma_score = .posterior_covariance(theta_prev, w_prev, kernel_t)

    log_q = .log_q_proposal_vectorized(
      theta_new,
      theta_prev,
      w_prev,
      Sigma_score
    )
  } else {
    log_q = log_prior
  }

  log_weight = log_prior + log_abc_kernel - log_q
  # log_weight = log_abc_kernel

  sim_df %>%
    dplyr::mutate(
      # Compute unnormalized ABC weights (Gaussian kernel)
      abc_weight = exp(log_weight - max(log_weight))
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
  n_sims,
  prev_sim_df,
  priors_list,
  kernel_t
) {
  theta_prev = prev_sim_df %>%
    dplyr::select(dplyr::starts_with("abc_mvn_")) %>%
    as.matrix()

  w_prev = suppressWarnings(prev_sim_df$abc_weight)
  if (is.null(w_prev)) {
    stop("Weights must be calculated first.")
  }
  Sigma_pert = .posterior_covariance(theta_prev, w_prev, kernel_t)

  # Step 1: resample indices with probability w_prev
  idx = sample(
    seq_len(nrow(theta_prev)),
    size = n_sims,
    replace = TRUE,
    prob = w_prev
  )

  # Step 2: perturb each selected particle in z space
  z = theta_prev[idx, , drop = FALSE] +
    mvtnorm::rmvnorm(
      n_sims,
      mean = rep(0, ncol(theta_prev)),
      sigma = Sigma_pert
    )

  sim_df = z %>% dplyr::as_tibble()

  nms = names(priors_list)
  nms = nms[nms != ""]
  for (i in seq_along(nms)) {
    nm = nms[[i]]
    mvn_nm = sprintf("abc_mvn_%s", nm)
    prior_dist = priors_list[[nm]]
    sim_df[[nm]] = prior_dist$q(stats::pnorm(z[, mvn_nm]))
  }

  return(sim_df)
}


#' Covariance matrix from the posteriors
#'
#' @param theta a set of particles as a matrix
#' @param weights the weights of the previous particles assumed calculated using
#'   the same kernel
#' @inheritParams tidyabc_common
#'
#' @returns a covariance matrix, optionally scaled
#' @keywords internal
.posterior_covariance = function(theta, weights, kernel_t) {
  if (is.null(weights)) {
    weights = rep(1 / nrow(theta), nrow(theta))
  }

  # Weighted covariance
  cov_w = stats::cov.wt(theta, wt = weights, method = "ML")$cov

  # Scale it (standard trick for random walk)
  d = ncol(theta) # dimension = 4
  scaling_factor = 2.38^2 / d # optimal for Gaussian targets
  cov_w = kernel_t * scaling_factor * cov_w

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
