#' `abc_prior` S3 class
#'
#' @name abc_prior
#' @param x an `abc_prior` S3 object
#' @param .dists distribution functions as a named list of S3 `dist_fns` objects
#' @param .constraints a list of one sided formulae the result each of which should
#'   evaluate to a boolean when compared against the names of the priors and
#'   derived values.
#' @param .derived a list of two sided formulae. The RHS refer to the
#'   priors, and the LHS as a name to derive.
#' @param .cor (optional) a correlation matrix for the priors
#' @param ... passed on to methods
#'
#' @concept abc_prior_s3
NULL


#' Construct a set of priors
#'
#' `abc_prior` S3 objects are used to hold the specification of prior and
#' intermediate proposal distributions. They are inputs to the main `abc_...()`
#' workflow functions.
#'
#' @param ... a list of formulae. Two sided will be interpreted as distribution
#'   or derived value specifications. One sided as constraints between parameters.
#'   A distribution is specified as the name of the family of statistical distributions
#'   and their parameters: e.g.: `x ~ norm(mean=3,sd=2)`. The name will be matched
#'   to the first hit on the search path.
#'
#' @inherit new_abc_prior return
#' @concept workflow
#' @export
#'
#' @examples
#' p = priors(
#'   mean ~ tidyabc::rgamma2(4,2),
#'   sd ~ gamma2(2,1),
#'   shape ~ mean^2 / sd^2,
#'   rate ~ mean / sd^2,
#'   ~ mean > sd
#' )
#'
#' print(p)
#'
#' # Plot methods are also provided:
#' if (interactive()) plot(p)
#'
#' # constraints:
#' p@constraints
#'
priors = function(...) {
  dots = rlang::list2(...)
  if (!all(sapply(dots, rlang::is_formula))) {
    stop("priors should be specified as a list of two sided formulae")
  }
  if (any(nzchar(names(dots)))) {
    stop("priors should not be named. Did you use `=` instead of `~`?")
  }

  # Constraints are one sided formulae
  c = sapply(dots, function(f) is.null(rlang::f_lhs(f)))
  constraints = dots[c]
  others = dots[!c]

  names = sapply(others, rlang::f_lhs)
  names = sapply(names, rlang::as_label)

  fns = lapply(others, function(x) {
    x = rlang::f_rhs(x)
    setdiff(all.vars(x, functions = TRUE), all.vars(x))
  })

  # TODO: Improve handling of qualified distribution function names
  # Issue URL: https://github.com/ai4ci/tidyabc/issues/7
  # at the moment this strips off package info and ignores it.
  # it would be better if we could somehow provide the package info to
  # as.dist_fns.character.

  # browser()

  # Stat functions are identifiable
  sfns = sapply(fns, function(x) {
    if (length(x) != 1) {
      if (x[[1]] == "::") {
        x = gsub("^(p|q|r|d)", "", x[[3]])
      } else {
        x = x[[1]]
      }
    }
    sfnnm = paste0(c("p", "q", "r", "d"), x)
    # browser()
    return(all(sapply(sfnnm, exists, mode = "function")))
  })

  derived = others[!sfns]
  stat_fns = others[sfns]

  dists = lapply(stat_fns, function(stat_fn) {
    param_nm = rlang::as_label(rlang::f_lhs(stat_fn))
    stat_expr = rlang::f_rhs(stat_fn)
    dist_nm = rlang::as_label(stat_expr[[1]])
    dist_nm = gsub("^[^:]+::[pqrd]", "", dist_nm)
    hyperparams = as.list(stat_expr[-1])
    dfn = as.dist_fns(dist_nm, params = hyperparams)
    # dfn$name = param_nm
    return(dfn)
  })
  names(dists) = names[sfns]

  return(new_abc_prior(
    .dists = dists,
    .constraints = constraints,
    .derived = derived
  ))
}

#' @describeIn abc_prior Create a new prior
#' @returns an S3 object of class `abc_prior` which contains
#' - a list of `dist_fns`
#' - a `cor` attribute describing their correlation
#' - a `derived` attribute describing derive values
#' - a `constraints` attribute listing the constraints
#' - a `params` attribute listing the names of the parameters
#' @concept abc_prior_s3
#' @unit
#' p = new_abc_prior(
#'   .dists = list(
#'     mean = as.dist_fns("norm",4,2),
#'     sd = as.dist_fns("gamma",2)
#'   ),
#'   .derived = list(
#'     shape ~ mean^2 / sd^2,
#'     rate ~ mean / sd^2
#'   ),
#'   .constraints = list(
#'     ~ mean > sd
#'   )
#' )
#'
#' testthat::expect_equal(
#'   format(p),
#'   "Parameters: \n* mean: norm(mean = 4, sd = 2)\n* sd: gamma(shape = 2, rate = 1)\nConstraints:\n* mean > sd\nDerived values:\n* shape = mean^2/sd^2\n* rate = mean/sd^2"
#' )
#'
new_abc_prior = function(
  .dists,
  .constraints = list(),
  .derived = list(),
  .cor = NULL
) {
  if (!is.list(.dists)) {
    stop("`.dists` must be a named list")
  }
  if (!(all(nzchar(names(.dists))))) {
    stop("All `.dists` must be named")
  }
  if (!all(sapply(.constraints, rlang::is_formula))) {
    stop("All `.constraints` must be a formula")
  }
  if (!all(sapply(.derived, rlang::is_formula))) {
    stop("All `.derived` must be a formula")
  }
  tmp = lapply(names(.dists), function(nm) {
    return(as.dist_fns(.dists[[nm]]))
  })
  names(tmp) = names(.dists)
  if (is.null(.cor)) {
    .cor = diag(length(tmp))
  }
  return(
    structure(
      tmp,
      class = c("abc_prior", "list"),
      constraints = .constraints,
      derived = .derived,
      params = names(tmp),
      cor = .cor
    )
  )
}


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
#' @concept abc_prior_s3
#'
#' @returns a data frame of samples in MVN (prefixed `abc_mvn_`) and parameter
#' space.
#' @keywords internal
#' @concept abc_prior_s3
#' @unit
#' p = new_abc_prior(
#'   .dists = list(
#'     mean = as.dist_fns("norm",4,2),
#'     sd = as.dist_fns("gamma",2)
#'   ),
#'   .derived = list(
#'     shape ~ mean^2 / sd^2,
#'     rate ~ mean / sd^2
#'   )
#' )
#'
#' s = .sample_priors(p,10)
#' testthat::expect_equal(
#'   colnames(s),
#'   c("abc_mvn_mean", "mean", "abc_mvn_sd", "sd")
#' )
.sample_priors = function(proposal_list, n_sims) {
  if (!is.abc_prior(proposal_list)) {
    stop("proposal list must be in the format of an `abc_prior`")
  }

  sim_df = dplyr::tibble(.rows = n_sims)

  nms = proposal_list@params
  cor = proposal_list@cor

  # generate the correlated measures as normals:
  z = mvtnorm::rmvnorm(n_sims, mean = rep(0, length(nms)), sigma = cor)

  # convert to uniform (pnorm) & map to targets (prior$q:
  for (i in seq_along(nms)) {
    nm = nms[[i]]
    mvn_nm = sprintf("abc_mvn_%s", nm)
    prior_dist = proposal_list[[nm]]
    sim_df[[mvn_nm]] = z[, i]
    sim_df[[nm]] = prior_dist$q(stats::pnorm(z[, i]))
  }

  return(sim_df)
}


#' Apply derived values and constraints to samples
#'
#' Repeatedly samples and removes invalid until the desired number of valid
#' samples is reached. This also calculates derived values
#'
#' @inheritParams common_internal
#' @inheritParams tidyabc_common
#' @concept abc_prior_s3
#' @param sampler_fn a function that creates random samples
#' @param ... passed on to `sampler_fn`
#'
#' @returns a data frame of samples in MVN (prefixed `abc_mvn_`) and parameter
#' space.
#' @keywords internal
#' @unit
#' p = new_abc_prior(
#'   .dists = list(
#'     mean = as.dist_fns("norm",4,2),
#'     sd = as.dist_fns("gamma",2)
#'   ),
#'   .derived = list(
#'     shape ~ mean^2 / sd^2,
#'     rate ~ mean / sd^2
#'   ),
#'   .constraints = list(
#'     ~ mean > sd
#'   )
#' )
#'
#' s = .sample_constrained(p,1000)
#'
#' testthat::expect_equal(
#'   colnames(s),
#'   c("abc_mvn_mean", "mean", "abc_mvn_sd", "sd", "shape", "rate")
#' )
#'
#' testthat::expect_equal(all(s$mean > s$sd), TRUE)
.sample_constrained = function(
  proposal_list,
  n_sims,
  sampler_fn = .sample_priors,
  ...
) {
  sim_df = dplyr::tibble()
  while (nrow(sim_df) < n_sims) {
    delta = n_sims - nrow(sim_df)
    tmp = sampler_fn(proposal_list, n_sims = delta, ...)

    # 2 sided formulae are derived values:
    for (form in proposal_list@derived) {
      tgt = rlang::f_lhs(form)
      expr = rlang::f_rhs(form)
      tmp = tmp %>% dplyr::mutate(!!tgt := !!expr)
    }

    for (constr in proposal_list@constraints) {
      expr = rlang::f_rhs(constr)
      tmp = tmp %>% dplyr::filter(!!expr)
    }

    sim_df = dplyr::bind_rows(sim_df, tmp)
  }
  return(sim_df)
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
#   nms = priors_list@params
#   param_nms = priors_list@params
#   derived = unname(priors_list[nms == ""])
#
#   # derived values:
#   for (form in priors_list@derived) {
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

# Boilerplate S3 methods: ----

#' @describeIn abc_prior Create a prior from a named list of `dist_fns`
#' @inheritDotParams new_abc_prior
#' @inherit new_abc_prior return
#' @export
#' @concept abc_prior_s3
as.abc_prior = function(x, ...) {
  if (is.abc_prior(x)) {
    return(x)
  }
  if (!is.list(x)) {
    stop("Only lists can be cast to `abc_priors`")
  }
  if (!all(nzchar(names(x)))) {
    stop("All items in the `x` list must be named.")
  }
  if (!all(sapply(x, is.dist_fns))) {
    stop("All items in the `x` list must be `dist_fns` objects")
  }
  return((new_abc_prior(x, ...)))
}

#' @describeIn abc_prior Test is an `abc_prior`
#' @export
#' @concept abc_prior_s3
is.abc_prior = function(x, ...) {
  return(inherits(x, "abc_prior"))
}

#' @describeIn abc_prior Format an `abc_prior`
#' @export
#' @concept abc_prior_s3
format.abc_prior = function(x, ...) {
  paste0(
    c(
      "Parameters: ",
      sapply(names(x), function(nm) {
        sprintf("* %s: %s", nm, x[[nm]]$name)
      }),
      if (length(x@constraints)) {
        c(
          "Constraints:",
          sapply(x@constraints, function(f) {
            sprintf("* %s", format(rlang::f_rhs(f)))
          })
        )
      } else {
        NULL
      },
      if (length(x@derived)) {
        c(
          "Derived values:",
          sapply(x@derived, function(f) {
            sprintf(
              "* %s = %s",
              format(rlang::f_lhs(f)),
              format(rlang::f_rhs(f))
            )
          })
        )
      } else {
        NULL
      }
    ),
    collapse = "\n"
  )
}

#' @describeIn abc_prior Print an `abc_prior`
#' @export
#' @concept abc_prior_s3
print.abc_prior = function(x, ...) {
  cat(format(x, ...))
}

#' @describeIn abc_prior Plot an `abc_prior`
#' @export
#' @concept abc_prior_s3
plot.abc_prior = function(x, ...) {
  for (nm in names(x)) {
    x[[nm]]$name = nm
  }
  plot(as.dist_fns_list(x)) + ggplot2::facet_wrap(~name)
}

#' Extract named attribute from a `abc_prior`
#'
#' @inheritParams abc_prior
#' @param y item to retrieve
#' @returns an attribute value for `x`
#' @export
#' @concept abc_prior_s3
#' @name at.dist_fns
`@.abc_prior` = function(x, y) {
  if (is.character(y)) {
    ylab = y
  } else {
    ylab = deparse(substitute(y))
  }
  return(attr(x, ylab))
}

#' Support for auto suggests on an `abc_prior`s
#' @inheritParams abc_prior
#' @param pattern a regular expression
#' @returns the names of the attributes
#' @exportS3Method utils::.AtNames abc_prior
#' @keywords internal
#' @concept abc_prior_s3
.AtNames.abc_prior = function(x, pattern) {
  return(utils::.DollarNames(attributes(x), pattern))
}
