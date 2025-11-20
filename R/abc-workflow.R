## Workflow wrappers ----

#' Perfom simple ABC rejection algorithm
#'
#' This function will execute a simulation for a random selection of parameters
#' and identify the best matching `acceptance_rate` percent, as defined by the
#' summary distance metric. A large number of simulations and a low acceptance
#' rate are best here.
#'
#' Parameters \eqn{\theta^{(i)}} are sampled independently from the prior
#' distribution \eqn{P(\theta)} for \eqn{i = 1, \dots, n_{\text{sims}}}.
#' For each \eqn{\theta^{(i)}}, simulated data \eqn{D_s^{(i)} = M(\theta^{(i)})}
#' is generated via the simulator function `sim_fn`, and a vector of summary
#' statistics \eqn{S_s^{(i)} = \texttt{scorer\_fn}(D_s^{(i)})}
#' is computed and compared to the observed summary statistics
#' \eqn{S_o = \texttt{scorer\_fn}(D_o)}.
#'
#' A distance metric \eqn{d^{(i)} = d(S_s^{(i)}, S_o)} is computed. By default,
#' this is the Euclidean distance:
#' \deqn{
#'   d^{(i)} = \left\| W \circ (S_s^{(i)} - S_o) \right\|_2,
#' }
#' where \eqn{W} is a vector of optional summary statistic weights
#' (`scoreweights`), and \eqn{\circ} denotes element-wise multiplication.
#' Other supported metrics include Manhattan (\eqn{\ell_1}) and Mahalanobis
#' distance (using the empirical covariance from the first wave).
#'
#' The tolerance threshold \eqn{\epsilon} is set to the
#' \eqn{\alpha = \texttt{acceptance\_rate}} quantile of the distances
#' \eqn{\{d^{(i)}\}_{i=1}^{n_{\text{sims}}}}:
#' \deqn{
#'   \epsilon = \text{quantile}\big(\{d^{(i)}\}, \alpha\big).
#' }
#'
#' Unnormalized ABC weights are then assigned using a kernel function
#' \eqn{K_\epsilon(\cdot)}:
#' \deqn{
#'   \tilde{w}^{(i)} = K_\epsilon\big(d^{(i)}\big),
#' }
#' where \eqn{K_\epsilon(d)} is one of the kernels defined in `kernels.R`
#' (e.g., Epanechnikov: \eqn{K_\epsilon(d) = \frac{3}{4\epsilon}
#' (1 - d^2\!/\!\epsilon^2)\, \mathbb{I}(d \leq \epsilon)}).
#' These weights are then transformed via a logistic ("expit") function and
#' normalized to sum to one:
#' \deqn{
#'   w^{(i)} = \frac{\text{expit}(\log \tilde{w}^{(i)})}{
#'     \sum_j \text{expit}(\log \tilde{w}^{(j)})}
#'   = \frac{ \tilde{w}^{(i)} / (1 + \tilde{w}^{(i)}) }{
#'     \sum_j \tilde{w}^{(j)} / (1 + \tilde{w}^{(j)}) }.
#' }
#' The resulting weighted sample \eqn{\{(\theta^{(i)}, w^{(i)})\}} approximates
#' the ABC posterior distribution \eqn{P_\epsilon(\theta \mid D_o)}.
#'
#' @inheritParams tidyabc_common
#' @param ... must be empty
#'
#' @inherit new_abc_fit return
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
  ...,
  converged_fn = default_termination_fn(),
  obsscores = NULL,
  distance_method = "euclidean",
  keep_simulations = FALSE,
  seed = NULL,
  parallel = FALSE,
  debug_errors = FALSE,
  kernel = "epanechnikov",
  scoreweights = NULL
) {
  rlang::check_dots_empty()
  message("ABC rejection, 1 wave.")

  priors_list = as.abc_prior(priors_list)

  prev_sim_df = .sample_constrained(
    priors_list,
    n_sims,
    sampler_fn = .sample_priors
  )
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
    parallel = parallel,
    debug = debug_errors,
    scoreweights = scoreweights
  )

  if (all(is.na(sim_df$abc_summary_distance))) {
    stop(
      "No non null scores from simulation. Rerun with debug flag to test `sim_fn` an `scorer_fn`."
    )
  }

  sim_df = sim_df %>%
    dplyr::filter(
      !is.na(abc_summary_distance) &
        is.finite(abc_summary_distance)
    )

  epsilon = stats::quantile(
    sim_df$abc_summary_distance,
    probs = acceptance_rate
  )

  sim_df = sim_df %>%
    .calculate_weights_wave_one(
      epsilon = epsilon,
      kernel = kernel
    )

  sim_df = sim_df %>%
    dplyr::filter(
      abc_weight > 0
    )

  # browser()
  # ggplot(sim_df)+geom_point(aes(x=abc_summary_distance,y = abc_weight))

  metric = .compare_waves(
    sim_df = sim_df,
    prev_sim_df = prev_sim_df,
    priors_list = priors_list,
    wave = 1
  )

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
#' Performs the ABC Sequential Monte Carlo (SMC) algorithm.
#' This iterative method refines parameter estimates across multiple waves.
#'
#' 1. **Initialization (Wave 1):**
#'    Parameters \eqn{\theta^{(i)}} are sampled from the prior \eqn{P(\theta)}.
#'    Simulations \eqn{D_s^{(i)} = M(\theta^{(i)})} are run, summaries \eqn{S_s^{(i)}}
#'    are computed, and distances \eqn{d^{(i)} = d(S_s^{(i)}, S_o)} are calculated.
#'    A tolerance threshold \eqn{\epsilon_1} is set as the
#'    \eqn{\alpha = \texttt{acceptance\_rate}} quantile of these initial distances.
#'    Unnormalized weights \eqn{\tilde{w}^{(i)}_1} are calculated using a kernel
#'    \eqn{K_{\epsilon_1}(d^{(i)})}.
#'
#' 2. **Subsequent Waves (\eqn{t > 1}):**
#'    - **Proposal Generation:** A particle \eqn{\theta^{(j)}_{t-1}} is selected
#'      from the previous wave's accepted particles with probability proportional
#'      to its weight \eqn{w^{(j)}_{t-1}}. The particle is then perturbed in a
#'      transformed MVN space using a multivariate normal kernel with covariance
#'      \eqn{\Sigma_t = \frac{\kappa_t^2}{d} \text{Cov}_{w_{t-1}}(\theta_{t-1})},
#'      where \eqn{\text{Cov}_{w_{t-1}}} is the weighted covariance from wave \eqn{t-1},
#'      \eqn{d} is the parameter dimension, and \eqn{\kappa_t} is the `kernel_t` parameter.
#'      The new proposal \eqn{\theta^{(i)}_t} is generated as:
#'      \deqn{
#'        \theta^{(i)}_t = \theta^{(j)}_{t-1} + \zeta, \quad \zeta \sim \mathcal{N}(0, \Sigma_t)
#'      }
#'      This proposal is mapped back to the original parameter space using the
#'      prior's copula transformation (from MVN space defined by prior CDFs).
#'
#'    - **Simulation and Weighting:** Simulations \eqn{D_s^{(i)} = M(\theta^{(i)}_t)}
#'      are run for the new proposals. Distances \eqn{d^{(i)}_t} are computed.
#'      The tolerance \eqn{\epsilon_t} is set as the \eqn{\alpha}-quantile of
#'      distances from the *current* wave's simulations. The unnormalized weight
#'      for particle \eqn{i} in wave \eqn{t} is calculated as:
#'      \deqn{
#'        \tilde{w}^{(i)}_t = \frac{P(\theta^{(i)}_t) K_{\epsilon_t}(d^{(i)}_t)}{q_t(\theta^{(i)}_t)}
#'      }
#'      where \eqn{P(\theta^{(i)}_t)} is the prior density, \eqn{K_{\epsilon_t}}
#'      is the ABC kernel, and \eqn{q_t(\theta^{(i)}_t)} is the proposal density
#'      from the previous wave's weighted particles (calculated using the
#'      perturbation kernel). This proposal density is computed as a weighted sum:
#'      \deqn{
#'        q_t(\theta^{(i)}_t) = \sum_j w^{(j)}_{t-1} \phi(\theta^{(i)}_t; \theta^{(j)}_{t-1}, \Sigma_t)
#'      }
#'      where \eqn{\phi(\cdot; \mu, \Sigma)} is the PDF of a multivariate normal
#'      with mean \eqn{\mu} and covariance \eqn{\Sigma}.
#'
#' 3. **Normalization:** Weights \eqn{w^{(i)}_t} are normalized to sum to one.
#'    Particles with negligible weights are typically filtered out.
#'
#' 4. **Termination:** The process repeats until a maximum number of waves or
#'    time is reached, or convergence criteria are met based on changes in
#'    parameter estimates or effective sample size (ESS).
#'
#' @inheritParams tidyabc_common
#' @param ... must be empty
#' @inherit new_abc_fit return
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
#'   max_time = 5, # 5 seconds to fit within examples limit
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
  ...,
  max_time = 5 * 60,
  converged_fn = default_termination_fn(),
  obsscores = NULL,
  distance_method = "euclidean",
  seed = NULL,
  parallel = FALSE,
  allow_continue = interactive(),
  debug_errors = FALSE,
  kernel = "epanechnikov",
  scoreweights = NULL
) {
  rlang::check_dots_empty()
  if (!is.null(seed)) {
    seed = set.seed(seed)
    on.exit(set.seed(seed), add = TRUE)
  }

  priors_list = as.abc_prior(priors_list)

  # Generate a set of initial proposals
  sim_df = .sample_constrained(priors_list, n_sims, sampler_fn = .sample_priors)
  # browser()
  message("ABC-SMC")
  bar_id = cli::cli_progress_bar(name = "SMC waves: ", total = max_time)
  prev_sim_df = NULL
  wave_df = NULL
  wave1_metrics = NULL
  summ_df = posterior_summarise(sim_df, priors_list) %>% dplyr::mutate(wave = 0)
  converged = FALSE
  i = 0
  start_time = Sys.time()
  finish_time = start_time + max_time
  while (Sys.time() < finish_time) {
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
      wave1_metrics = wave1_metrics,
      debug = debug_errors,
      scoreweights = scoreweights,
      wave = i
    )

    sim_df = sim_df %>%
      dplyr::filter(
        !is.na(abc_summary_distance) &
          is.finite(abc_summary_distance)
      )

    # Extract the covariance if it has been calculated
    # This stays the same throughout all waves
    if (is.null(wave1_metrics)) {
      wave1_metrics = posterior_distance_metrics(
        sim_df$abc_component_score,
        obsscores = obsscores,
        keep_data = FALSE
      )
    }

    if (i > 1) {
      sim_df = dplyr::bind_rows(
        sim_df,
        prev_sim_df %>% dplyr::select(dplyr::all_of(colnames(sim_df)))
      )
    }

    epsilon = stats::quantile(
      sim_df$abc_summary_distance,
      probs = acceptance_rate
    )

    # A major problem here is that the weights are calculated as K(distance) *
    # P_prior(Particle) / P_proposal(Particle) This can have issues because the
    # 1/proposal part is tends to push away from convergence. If the K(distance)
    # function is not steep enough (because we are a long way away for the
    # correct answer) the resulting weight distribution is U shaped away from
    # where the current proposal is. For uniform kernels This is probably less
    # of an issue because it is effectively the same as the Hard SMC cutoff,
    # similarly for kernels with zero weight in the tails. The gaussian kernel
    # however can be quite flat if the distances are large. This can lead
    # to degeneracy collapse. We ahve mitigated this somewhat by calculating the
    # weight as the logit (the prior/proposal term is effectively an odds) rather
    # than the plain exponent. This tends to keep the distribution within range

    sim_df = sim_df %>%
      .calculate_weights_smc(
        epsilon = epsilon,
        prev_sim_df = prev_sim_df,
        kernel = kernel
      )

    sim_df = sim_df %>%
      dplyr::filter(
        abc_weight > sqrt(.Machine$double.eps)
      )

    metric = .compare_waves(
      sim_df = sim_df,
      prev_sim_df = prev_sim_df,
      priors_list = priors_list,
      wave = i
    )

    wave_df = dplyr::bind_rows(wave_df, metric)

    summ_df = dplyr::bind_rows(
      summ_df,
      posterior_summarise(sim_df, priors_list) %>% dplyr::mutate(wave = i)
    )

    prev_sim_df = sim_df

    if (
      converged_fn(
        wave = i,
        summary = metric$summary[[1]],
        per_param = metric$per_param[[1]]
      )
    ) {
      converged = TRUE
      message("Converged on wave: ", i)
      break
    }

    if (Sys.time() > finish_time) {
      if (allow_continue) {
        message("Current estimates: ")
        cat(.format_summ(summ_df), sep = "\n")
        if (
          utils::askYesNo(
            "Still converging. Do you want to continue?",
          )
        ) {
          finish_time = Sys.time() + max_time
        } else {
          break
        }
      } else {
        break
      }
    }

    sim_df = .sample_constrained(
      priors_list,
      n_sims = n_sims,
      sampler_fn = .sample_using_perturbation,
      prev_sim_df = prev_sim_df
    )

    # browser()
    cli::cli_progress_update(
      id = bar_id,
      status = sprintf("wave %d", i),
      set = as.double(
        min(c(Sys.time(), finish_time)) - start_time,
        units = "secs"
      ),
      total = as.double(finish_time - start_time, units = "secs")
    )
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
#' Performs the ABC Adaptive algorithm.
#' This iterative method refines parameter estimates across waves by fitting
#' empirical proposal distributions to the weighted posterior samples from the
#' previous wave. Unlike ABC-SMC, which uses a fixed perturbation kernel,
#' `abc_adaptive` constructs a new proposal distribution \eqn{Q_t(\theta)} at
#' each wave \eqn{t}.
#'
#' 1. **Initialization (Wave 1):**
#'    Parameters \eqn{\theta^{(i)}} are sampled from the prior \eqn{P(\theta)}.
#'    Simulations are run, summary statistics \eqn{S_s^{(i)}} are computed,
#'    and distances \eqn{d^{(i)} = d(S_s^{(i)}, S_o)} are calculated.
#'    A tolerance threshold \eqn{\epsilon_1} is set as the
#'    \eqn{\alpha = \texttt{acceptance\_rate}} quantile of these distances.
#'    Unnormalized weights \eqn{\tilde{w}^{(i)}_1} are calculated using a kernel
#'    \eqn{K_{\epsilon_1}(d^{(i)})}.
#'
#' 2. **Subsequent Waves (\eqn{t > 1}):**
#'    - **Proposal Generation:** An empirical joint proposal distribution
#'      \eqn{Q_t(\theta)} is constructed from the weighted posterior sample
#'      \eqn{\{(\theta^{(i)}_{t-1}, w^{(i)}_{t-1})\}} of the previous wave.
#'      This is done by fitting marginal empirical distributions
#'      \eqn{Q_{t,j}(\theta_j)} to each parameter \eqn{\theta_j}, using the
#'      `empirical()` function with the prior \eqn{P_j(\theta_j)} as a link to
#'      enforce support constraints. These marginals are assumed independent,
#'      but their weighted correlation matrix \eqn{R_t} is retained as an
#'      attribute and used to induce dependence in the MVN sampling space.
#'      New proposals \eqn{\theta^{(i)}_t} are generated by:
#'      \itemize{
#'        \item Sampling a vector \eqn{Z \sim \mathcal{N}(0, R_t)} in a
#'              correlated standard normal space.
#'        \item Mapping each component \eqn{Z_j} to uniform space:
#'              \eqn{U_j = \Phi(Z_j)}.
#'        \item Mapping to the parameter space using the empirical quantile
#'              functions: \eqn{\theta^{(i)}_{t,j} = Q_{t,j}^{-1}(U_j)}.
#'      }
#'
#'    - **Simulation and Weighting:** Simulations are run for the new proposals.
#'      Distances \eqn{d^{(i)}_t} are computed and the tolerance \eqn{\epsilon_t}
#'      is set as the \eqn{\alpha}-quantile of the current wave's distances.
#'      The unnormalized weight for particle \eqn{i} in wave \eqn{t} is calculated as:
#'      \deqn{
#'        \tilde{w}^{(i)}_t = \frac{P(\theta^{(i)}_t) K_{\epsilon_t}(d^{(i)}_t)}{Q_t(\theta^{(i)}_t)}
#'      }
#'      where \eqn{P(\theta^{(i)}_t) = \prod_j P_j(\theta^{(i)}_{t,j})} is the prior
#'      density (assuming independence), \eqn{K_{\epsilon_t}} is the ABC kernel,
#'      and \eqn{Q_t(\theta^{(i)}_t) = \prod_j Q_{t,j}(\theta^{(i)}_{t,j})} is the
#'      empirical proposal density (also assuming independence for density
#'      calculation, consistent with the marginal fitting). The correlation
#'      structure is handled in the sampling process, and optionally in the density
#'      evaluation.
#'
#' 3. **Normalization:** Weights \eqn{w^{(i)}_t} are normalized to sum to one.
#'    The algorithm includes a recovery mechanism: if the Effective Sample Size
#'    (ESS) falls below a threshold (e.g., 200), the acceptance rate is
#'    increased (i.e., \eqn{\epsilon_t} is made larger) to accept more particles
#'    and improve the ESS.
#'
#' 4. **Termination:** The process repeats until a maximum time is reached or
#'    convergence criteria based on parameter stability and credible interval
#'    contraction are met.
#'
#' @inheritParams tidyabc_common
#' @inheritParams empirical
#' @param ... must be empty

#'
#' @inherit new_abc_fit return
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
#'   max_time = 5, # 5 seconds to fit within examples limit
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
  ...,
  max_time = 5 * 60,
  converged_fn = default_termination_fn(),
  obsscores = NULL,
  distance_method = "euclidean",
  seed = NULL,
  knots = NULL,
  parallel = FALSE,
  max_recover = 3,
  allow_continue = interactive(),
  debug_errors = FALSE,
  kernel = "epanechnikov",
  bw = 0.1,
  widen_by = 1.05,
  scoreweights = NULL,
  use_proposal_correlation = TRUE
) {
  rlang::check_dots_empty()
  if (!is.null(seed)) {
    seed = set.seed(seed)
    on.exit(set.seed(seed), add = TRUE)
  }

  priors_list = as.abc_prior(priors_list)

  # browser()
  message("ABC-Adaptive")
  bar_id = cli::cli_progress_bar(name = "Adaptive waves: ", total = max_time)

  # Setup first wave

  wave_df = NULL
  wave1_metrics = NULL
  # N.b. could be refactored inside the loop?
  # Generate a set of initial proposals
  sim_df = .sample_constrained(priors_list, n_sims, sampler_fn = .sample_priors)
  summ_df = posterior_summarise(sim_df, priors_list) %>% dplyr::mutate(wave = 0)
  converged = FALSE
  i = 0
  start_time = Sys.time()
  finish_time = start_time + max_time
  proposal_list = priors_list
  prev_sim_df = sim_df
  recover = 0

  while (Sys.time() < finish_time) {
    i = i + 1

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
      wave1_metrics = wave1_metrics,
      debug = debug_errors,
      scoreweights = scoreweights,
      wave = i
    )

    # Exclude failed simulations and infinite distances:
    sim_df = sim_df %>%
      dplyr::filter(
        !is.na(abc_summary_distance) &
          is.finite(abc_summary_distance)
      )

    while (recover <= max_recover) {
      # combine with previous particles.
      # before refining kernel
      if (i > 1) {
        tmp_sim_df = dplyr::bind_rows(
          sim_df,
          prev_sim_df %>% dplyr::select(dplyr::all_of(colnames(sim_df)))
        )
      } else {
        tmp_sim_df = sim_df
      }

      epsilon = stats::quantile(
        sim_df$abc_summary_distance,
        probs = acceptance_rate
      )

      tmp_sim_df = tmp_sim_df %>%
        .calculate_weights_adaptive(
          priors_list = priors_list,
          epsilon = epsilon,
          proposal_list = proposal_list,
          kernel = kernel,
          use_proposal_correlation = use_proposal_correlation
        )

      tmp_sim_df = tmp_sim_df %>%
        dplyr::filter(
          abc_weight > sqrt(.Machine$double.eps)
        )

      # ggplot(sim_df, aes(x = abc_weight, y = abc_summary_distance)) +
      #   geom_point()

      ESS = 1 / sum(tmp_sim_df$abc_weight^2)
      # Exit the retry loop if ESS is big enough
      if (ESS > 200) {
        sim_df = tmp_sim_df
        break
      }

      # ESS is too small. we will try and recover by
      # increasing epsilon
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
        acceptance_rate = 1 - (1 - acceptance_rate) / 2
        message(
          "Attempting recovery with larger acceptance rate: ",
          sprintf("%1.3f%%", acceptance_rate * 100)
        )
      }
    }

    # End of recovery loop:
    # Extract the covariance if it has been calculated
    # This stays the same throughout all waves
    if (is.null(wave1_metrics)) {
      wave1_metrics = posterior_distance_metrics(
        sim_df$abc_component_score,
        obsscores,
        keep_data = FALSE
      )
    }

    metric = .compare_waves(
      sim_df = sim_df,
      prev_sim_df = prev_sim_df,
      priors_list = priors_list,
      wave = i
    ) %>%
      dplyr::mutate(
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
        wave = i,
        summary = metric$summary[[1]],
        per_param = metric$per_param[[1]]
      )
    ) {
      converged = TRUE
      message("Converged on wave: ", i)
      break
    }

    if (Sys.time() > finish_time) {
      if (allow_continue) {
        message("Current estimates: ")
        cat(.format_summ(summ_df), sep = "\n")
        if (
          utils::askYesNo(
            "Still converging. Do you want to continue?",
          )
        ) {
          finish_time = Sys.time() + max_time
        } else {
          break
        }
      } else {
        break
      }
    }

    proposal_list = posterior_fit_empirical(
      sim_df,
      priors_list = priors_list,
      # priors_list = proposal_list,
      # might think this could be proposal_list but either seems to work...
      knots = knots,
      bw = bw,
      widen_by = widen_by
    )

    sim_df = .sample_constrained(
      proposal_list,
      n_sims,
      sampler_fn = .sample_priors
    )

    cli::cli_progress_update(
      id = bar_id,
      status = sprintf("wave %d", i),
      set = as.double(
        min(c(Sys.time(), finish_time)) - start_time,
        units = "secs"
      ),
      total = as.double(finish_time - start_time, units = "secs")
    )
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
#' last_wave_metrics = utils::tail(fit$waves,1)
#' converged = check(
#'   wave = 0,
#'   summary = last_wave_metrics$summary[[1]],
#'   per_param = last_wave_metrics$per_param[[1]]
#' )
#'
#' if (isTRUE(converged)) print("Converged (permissive definition)")
#'
default_termination_fn = function(stability = 0.01, confidence = 0.1) {
  function(wave, summary, per_param) {
    stability = .recycle_and_name(stability, per_param$param)
    confidence = .recycle_and_name(confidence, per_param$param)

    isTRUE(
      all(per_param$rel_mean_change < stability) &&
        all(per_param$IQR_95_redn < confidence)
    )
  }
}


#' Run the SMC or adaptive algorithm for a set number of waves
#'
#' @param max_wave the number of waves to run
#'
#' @returns A function that will test for completion
#' @export
#'
#' @examples
#' # Declare converged after 3 waves:
#' converged_fn = fixed_wave_termination_fn(3)
fixed_wave_termination_fn = function(max_wave) {
  # TODO: The termination function should have mulitple possible statuses
  # Issue URL: https://github.com/ai4ci/tidyabc/issues/29
  function(wave, summary, per_param) {
    if (wave < max_wave) {
      return(FALSE)
    }
    return(TRUE) # sprintf("Finished at wave %d",wave))
  }
}

#' Run the simulation for one set of parameters
#'
#' @inheritParams tidyabc_common
#' @param ... simulation parameters, must be named
#' @param params a named list of simulation parameters to test (as an alternative
#'   to including them in `...`)
#' @param debug start the simulation function in debug mode. This will step
#'   through both the `sim_fn` and the `scorer_fn` line by line to check that
#'   the behaviour is as intended.
#'
#' @returns a list containing the parameters as `truth` and an
#'   instance of the simulation as `obsdata`, `obsscores` is the result of
#'   comparing the `obsdata` with itself and is usually going to result in
#'   zeros. Both `sim_fn` and `scorer_fn` are rewritten to make sure that all
#'   package names are fully qualified. The rewritten versions are refurned in
#'   the result list (as `sim_fn` and `scorer_fn` respectively)
#' @export
#' @concept workflow
test_simulation = function(
  sim_fn,
  scorer_fn,
  ...,
  params = NULL,
  obsdata = NULL,
  seed = NULL,
  debug = FALSE,
  debug_errors = !debug
) {
  # TODO: add in check for score weights
  # Issue URL: https://github.com/ai4ci/tidyabc/issues/1

  if (!is.null(seed)) {
    seed = set.seed(seed)
    on.exit(set.seed(seed), add = TRUE)
  }

  sim_fn = .crate_sim_fn(sim_fn, debug = debug_errors)

  if (is.null(params)) {
    params = rlang::list2(...)
    stopifnot(all(nzchar(names(params))))
  }

  int_sim_fn = rlang::fn_env(sim_fn)$sim_fn
  if (debug) {
    debug(int_sim_fn)
    on.exit(undebug(int_sim_fn), add = TRUE)
  }

  test = do.call(sim_fn, params)
  if (is.null(obsdata)) {
    obsdata = test
  }

  # debug= = !debug here bacuse we are setting up an interactive session
  # anyway but we always want to
  scorer_fn = .crate_scorer_fn(
    scorer_fn,
    obsdata = obsdata,
    debug = debug_errors
  )

  int_scorer_fn = rlang::fn_env(scorer_fn)$scorer_fn
  if (debug) {
    debug(int_scorer_fn)
    on.exit(undebug(int_scorer_fn), add = TRUE)
  }

  scores = scorer_fn(test)

  return(list(
    obsdata = test,
    truth = unlist(params),
    sim_fn = format(unclass(int_sim_fn)),
    obsscores = scores,
    scorer_fn = format(unclass(int_scorer_fn))
  ))
}


#' Generate a set of metrics from component scores
#'
#' The component scores are summary statistics output by the user supplied
#' `scorer_fn` as a named list. These can be variable in scale and location and
#' various options exist for combining them. They may need to be weighted by
#' scale as well as importance to get a model that works well. Such weights can
#' be input into the ABC algortihms using the `scoreweights` parameter. This
#' function helps provide diagnostics for calibrating the `scoreweights`
#' parameter.
#'
#' Given a list of component scores \eqn{S_s^{(i)} = (s_1^{(i)}, \dots, s_m^{(i)})}
#' from \eqn{n} simulations and observed scores \eqn{S_o = (s_1^{(o)}, \dots, s_m^{(o)})},
#' this function calculates:
#'
#' - The mean \eqn{\mu_j} and standard deviation \eqn{\sigma_j} of each score component \eqn{j}:
#'   \deqn{
#'     \mu_j = \frac{1}{n}\sum_{i=1}^n s_j^{(i)}, \quad
#'     \sigma_j = \sqrt{\frac{1}{n-1}\sum_{i=1}^n (s_j^{(i)} - \mu_j)^2}
#'   }
#'
#' - The covariance matrix \eqn{\Sigma} of the score components.
#'
#' - The root mean squared difference (RMSD) between simulated and observed scores for each component:
#'   \deqn{
#'     \text{RMSD}_j = \sqrt{\frac{1}{n}\sum_{i=1}^n (s_j^{(i)} - s_j^{(o)})^2}
#'   }
#'
#' - A recommended vector of `scoreweights` \eqn{w_j}, calculated as the ratio of the component's
#'   standard deviation to its RMSD, normalized to sum to 1:
#'   \deqn{
#'     w_j = \frac{\sigma_j / \text{RMSD}_j}{\sum_{k=1}^m (\sigma_k / \text{RMSD}_k)}
#'   }
#'   These weights help balance the influence of different summary statistics in
#'   the overall distance metric, especially when using Euclidean distance.
#'
#' @inheritParams  tidyabc_common
#' @param keep_data mainly for internal use this flag gives you the component
#'  scores as a matrix
#'
#' @returns a list containing the following items:
#' - `obsscores`: the input reference scores for each component
#' - `means`, `sds`: the means and sds of each score component
#' - `cov`: the covariance matrix for the scores
#' - `mad`: the mean absolute differences between the `simscores` and the `obsscores`
#' - `rmsd`: the root mean squared differences between the `simscores` and the `obsscores`
#' - `scoreweights`: a the `sds` divided by the `rmsd`. This weight should
#' mean that the weights of the  individual summary scores have similar influence
#' in the overall `abc_summary_distance` output once combined during SMC and
#' adaptive waves, especially if euclidean distances are involved.
#' - `simscores`: (if `keep_data`) a matrix of all the scores from these input
#' simulation posteriors
#' - `deltascores`: (if `keep_data`)  a matrix of the differences between `simscores` and `obsscores`.
#'
#' @export
#' @concept workflow
#' @examples
#'
#' fit = example_rejection_fit()
#' metrics = posterior_distance_metrics(fit$posteriors)
#' # other elements available:
#' metrics$scoreweights
posterior_distance_metrics = function(
  posteriors_df,
  obsscores = NULL,
  keep_data = FALSE
) {
  if (inherits(posteriors_df, "abc_fit")) {
    posteriors_df = posteriors_df$posteriors
  }

  if (!is.data.frame(posteriors_df)) {
    # when we use this internally we pass a list of scores direct to it.
    component_scores = posteriors_df
  } else {
    component_scores = posteriors_df$abc_component_score
  }

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
  colnames(simscores) = names(component_scores[[1]])

  # Get obscores and check in correct format:
  if (is.null(obsscores)) {
    obsscores = stats::setNames(rep(0, ncol(simscores)), colnames(simscores))
  }
  if (is.null(names(obsscores))) {
    names(obsscores) = colnames(simscores)
  }
  if (!identical(sort(names(obsscores)), sort(colnames(simscores)))) {
    stop(
      "`obsscores` names do not match `scorer_fn` outputs. Should be a vector with names: ",
      paste0(colnames(simscores), collapse = ",")
    )
  }
  obsscores = unlist(obsscores)
  obsscores = obsscores[colnames(simscores)]

  means = apply(simscores, MARGIN = 2, FUN = mean)
  sds = apply(simscores, MARGIN = 2, FUN = stats::sd)
  cov = stats::cov(simscores)
  deltascores = simscores -
    matrix(
      rep(obsscores, nrow(simscores)),
      ncol = ncol(simscores),
      byrow = TRUE
    )
  mad = apply(abs(deltascores), MARGIN = 2, mean)
  rmsd = sqrt(apply(deltascores^2, MARGIN = 2, mean))
  scoreweights = sds / rmsd # 1 / rmsd / sum(1 / rmsd)
  scoreweights = scoreweights / sum(scoreweights)

  return(list(
    obsscores = obsscores,
    simscores = if (keep_data) simscores else NULL,
    deltascores = if (keep_data) deltascores else NULL,
    means = means,
    sds = sds,
    mad = mad,
    rmsd = rmsd,
    cov = cov,
    scoreweights = scoreweights
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
#' summ = posterior_summarise(fit$posteriors, fit$priors)
#'
#' summ %>% dplyr::glimpse()
#'
posterior_summarise = function(
  posteriors_df,
  priors_list,
  p = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
) {
  priors_list = as.abc_prior(priors_list)

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
    dplyr::group_modify(function(d, g, ...) {
      density = empirical(
        x = d$value,
        w = weights,
        name = g$param,
        link = if (is.null(priors_list[[g$param]])) {
          # Use numeric constructor?
          d$value
        } else {
          priors_list[[g$param]]
        },
        knots = 20,
        bw = 0.3
      )

      return(dplyr::tibble(
        mean = wmean(x = d$value, w = weights),
        sd = wsd(x = d$value, w = weights),
        quantiles = list(
          tibble::tibble(
            p = p,
            x = unname(density$q(p))
          )
        ),
        q.0.025 = density$q(0.025),
        q.0.5 = density$q(0.5),
        q.0.975 = density$q(0.975),
        density = c(density),
        ESS = 1 / sum(weights^2)
      ))
    })

  return(summ_df)
}


#' Fit empirical distribution to posterior samples for generating more waves
#'
#' This function allows "updating" of the prior with an empirical posterior
#' distribution which will retain the bounds of the prior, and is effectively a
#' spline based transform of the prior distribution, based on the density of
#' data in prior space. This gives a clean density when data is close to a
#' prior distribution limit and work better than a standard density.
#'
#' This takes weighted posterior samples \eqn{\{(\theta^{(i)}, w^{(i)})\}}
#' for each parameter \eqn{\theta_j} and constructs an empirical distribution
#' function \eqn{Q_j(\theta_j)} that approximates the posterior marginal for
#' that parameter.
#'
#' For each parameter \eqn{\theta_j}, the empirical distribution \eqn{Q_j} is
#' fitted using the `empirical()` function. The fitting is performed on the
#' weighted samples \eqn{(\theta_{j}^{(i)}, w^{(i)})}. A key feature is the use
#' of a link function \eqn{h_j(\cdot)} during the fitting process. This link
#' function is typically taken from the original prior distribution \eqn{P_j(\theta_j)}
#' (i.e., \eqn{h_j = P_j}). This ensures that the fitted empirical distribution
#' \eqn{Q_j} respects the support constraints defined by the prior (e.g., positive
#' values for a log-normal prior). The fitting effectively occurs in the transformed
#' space defined by the link: \eqn{h_j(\theta_j^{(i)})}.
#'
#' The resulting empirical distribution \eqn{Q_j} is a `dist_fns` object that
#' includes the support constraints from the original prior. The set of all fitted
#' marginal distributions \eqn{\{Q_j\}} forms the new proposal list for the next
#' wave.
#'
#' Additionally, the weighted covariance matrix \eqn{R} of the samples in the
#' MVN space (defined by the prior copula) is calculated:
#' \deqn{
#'   R = \text{Cov}_{w}(\Phi^{-1}(P_1(\theta_1)), \dots, \Phi^{-1}(P_d(\theta_d)))
#' }
#' where \eqn{\Phi^{-1}} is the quantile function of the standard normal.
#' This covariance matrix is stored as an attribute (`"cor"`) of the returned
#' proposal list and is used to induce correlation structure when sampling
#' new proposals from the empirical distributions in subsequent waves.
#'
#' @inheritParams tidyabc_common
#'
#' @returns an `abc_prior` S3 object approximating the distribution of the
#'   posterior samples, and adhering to the support of the provided priors.
#' @export
#' @concept workflow
#' @examples
#'
#' fit = example_smc_fit()
#' proposals = posterior_fit_empirical(fit$posteriors, fit$priors)
#'
#' proposals
#'
posterior_fit_empirical = function(
  posteriors_df,
  priors_list,
  knots = NULL,
  bw = 0.1,
  widen_by = 1
) {
  weights = suppressWarnings(posteriors_df$abc_weight) #maybe null

  priors_list = as.abc_prior(priors_list)
  nms = priors_list@params

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
    name = sprintf(
      "posterior [%1.2f \u00B1 %1.2f]",
      wmean(col, weights),
      wsd(col, weights)
    )
    # limits = prior$q(c(0, 1))
    tmp2 = empirical(
      x = col,
      w = weights,
      fit_spline = TRUE,
      knots = knots,
      link = prior,
      name = name,
      bw = bw
    )
    if (!isFALSE(widen_by)) {
      tmp2 = widen(tmp2, scale = widen_by, name = name)
    }
    return(tmp2)
  })
  names(tmp) = nms

  # Get particles in MVN space
  theta = posteriors_df %>%
    dplyr::select(dplyr::starts_with("abc_mvn_")) %>%
    as.matrix()
  # weighted correlation:
  if (is.null(weights)) {
    cov = stats::cov(theta)
  } else {
    cov = stats::cov.wt(theta, wt = weights, method = "ML")$cov
  }

  proposal_list = utils::modifyList(priors_list, tmp)
  attr(proposal_list, "cor") = stats::cov2cor(cov)

  return(proposal_list)
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
#   nms = nms = priors_list@params
#
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
