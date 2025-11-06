#' `abc_fit` S3 class
#'
#' A class holding the output of a single ABC model fitting, either as a
#' single ABC rejection round or after a set of SMC waves.
#'
#' @param type the type of ABC algorithm
#' @param iterations number of completed iterations
#' @param converged boolean - did the result meet convergence criteria
#' @param wave_df a list of dataframes of wave convergence metrics
#' @param summ_df the summary of the parameter fits after each wave.
#' @param sim_df the final wave posteriors
#' @param x a `abc_fit` object as output by the `abc_XXX` functions
#' @param truth a named numeric vector of known parameter values
#' @param ... passed on to methods
#'
#' @returns an S3 ojbect of class `abc_fit` this contains the following:
#' - type: the type of ABC algorithm
#' - iterations: number of completed iterations
#' - converged: boolean - did the result meet convergence criteria
#' - waves: a list of dataframes of wave convergence metrics
#' - summary: a dataframe with the summary of the parameter fits after each wave.
#' - posteriors: the final wave posteriors
#'
#' @name abc_fit
#' @concept abc_fit_s3
NULL

#' Create a `abc_fit` object
#'
#' @inherit abc_fit
#' @keywords internal
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
summary.abc_fit = function(x, ...) {
  cat(paste0(
    c(
      format(x),
      "Parameter estimates:",
      .format_summ(x$summary)
    ),
    collapse = "\n"
  ))
}

.format_summ = function(summ_df) {
  format(
    summ_df %>%
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
  )
}

# TODO: format convergence metrics in a usable form
# consider changing comparison stats to first wave.

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
    mapping = ggplot2::aes(fill = group, colour = group)
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
    tmp = tmp +
      ggplot2::geom_vline(
        data = tibble::enframe(unlist(truth)),
        mapping = ggplot2::aes(xintercept = value),
        linetype = "dashed"
      )
  }

  tmp
}
