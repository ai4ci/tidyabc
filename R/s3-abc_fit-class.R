#' `abc_fit` S3 class
#'
#' A class holding the output of a single ABC model fitting, either as a
#' single ABC rejection round or after a set of SMC waves.
#'
#' @param type the type of ABC algorithm
#' @param iterations number of completed iterations
#' @param converged boolean - did the result meet convergence criteria
#' @param wave_df a list of dataframes of wave convergence metrics
#' @param priors_list a named list of priors specified as a `abc_prior` S3 object
#'   (see `priors()`), this can include derived values as unnamed 2-sided
#'   formulae, where the LHS of the formula will be assigned to the value of the
#'   RHS, plus optionally a set of constraints as one sided formulae where the
#'   RHS of the formulae will resolve to a boolean value.
#' @param summ_df the summary of the parameter fits after each wave.
#' @param sim_df the final wave posteriors
#' @param x a `abc_fit` object as output by the `abc_XXX` functions
#' @param object a `abc_fit` object as output by the `abc_XXX` functions
#' @param truth a named numeric vector of known parameter values
#' @param ... passed on to methods
#'
#' @name abc_fit
#' @concept abc_fit_s3
NULL

#' @describeIn abc_fit Create a `abc_fit` object
#'
#' @returns an S3 object of class `abc_fit` this contains the following:
#' - type: the type of ABC algorithm
#' - iterations: number of completed iterations
#' - converged: boolean - did the result meet convergence criteria
#' - waves: a list of dataframes of wave convergence metrics
#' - summary: a dataframe with the summary of the parameter fits after each wave.
#' - priors: the priors for the fit as a `abc_prior` S3 object
#' - posteriors: the final wave posteriors
#'
#' @concept abc_fit_s3
new_abc_fit = function(
  type,
  iterations,
  converged,
  priors_list,
  wave_df,
  summ_df,
  sim_df
) {
  return(
    structure(
      list(
        type = type,
        iterations = iterations,
        converged = converged,
        priors = priors_list,
        waves = wave_df,
        summary = summ_df,
        posteriors = sim_df
      ),
      class = c("abc_fit", "list")
    )
  )
}

#' @describeIn abc_fit S3 format method
#' @export
format.abc_fit = function(x, ...) {
  if (x$iterations > 1) {
    sprintf(
      "ABC %s fit: %d waves - (%s)",
      x$type,
      x$iterations,
      if (x$converged) "converged" else "not yet converged"
    )
  } else {
    sprintf(
      "ABC %s fit: single wave",
      x$type
    )
  }
}

#' @describeIn abc_fit S3 summary method
#' @export
summary.abc_fit = function(object, ..., truth = NULL) {
  if (!is.null(truth)) {
    truth = truth[names(object$priors)]
  }
  cat(paste0(
    c(
      format(object),
      "Parameter estimates:",
      .format_summ(object$summary, truth)
    ),
    collapse = "\n"
  ))
}

.format_summ = function(summ_df, truth = NULL) {
  mean = sd = NULL
  tmp = summ_df %>%
    dplyr::filter(wave == max(wave)) %>%
    dplyr::transmute(
      param,
      mean_sd = sprintf("%1.3f \u00B1 %1.3f", mean, sd),
      median_95_CrI = sprintf(
        "%1.3f [%1.3f \u2014 %1.3f]",
        q.0.5,
        q.0.025,
        q.0.975
      ),
      ESS
    )

  if (!is.null(truth)) {
    truth = tibble::enframe(unlist(truth), name = "param", value = "true")
    tmp = tmp %>% dplyr::left_join(truth, by = "param")
  }

  format(
    tmp
  )
}

#' @describeIn abc_fit S3 summary method
#' @export
tidy.abc_fit = function(x, ...) {
  x$summary %>%
    dplyr::filter(wave == max(wave))
}

#' @describeIn abc_fit S3 print method
#' @export
print.abc_fit = function(x, ...) {
  cat(summary(x, ...))
}

#' @describeIn abc_fit S3 plot method
#' @inheritDotParams plot.dist_fns_list
#' @export
plot.abc_fit = function(x, ..., truth = NULL) {
  fits = x$summary %>%
    dplyr::filter(wave == max(wave)) %>%
    dplyr::pull(density)

  tmp = lapply(names(x$priors), function(nm) {
    tmp = x$priors[[nm]]
    tmp$name = nm
    return(tmp)
  })

  fits = c(as.dist_fns_list(tmp), fits)

  tmp = plot(
    fits,
    plot_quantiles = TRUE,
    mapping = ggplot2::aes(fill = group, colour = group),
    ...
  ) +
    ggplot2::guides(
      fill = ggplot2::guide_none(),
      colour = ggplot2::guide_none()
    ) +
    ggplot2::facet_wrap(~name, scales = "free") +
    ggplot2::scale_color_discrete(
      palette = c("grey30", "red"),
      aesthetics = c("fill", "colour")
    )

  if (!is.null(truth)) {
    truth = truth[unique(x$summary$param)]

    tmp = tmp +
      ggplot2::geom_vline(
        data = tibble::enframe(unlist(truth)),
        mapping = ggplot2::aes(xintercept = value),
        linetype = "dashed"
      )
  }

  tmp
}
