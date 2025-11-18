#' Plot convergence metrics by wave for SMC and adaptive ABC
#'
#' @inheritParams tidyabc_common
#'
#' @returns a patchwork plot of convergence metrics
#' @export
#' @concept workflow
#'
#' @examples
#'
#' plot_convergence(
#'   example_adaptive_fit()
#' )
#'
plot_convergence = function(fit) {
  tmp = fit$waves %>%
    dplyr::select(-dplyr::any_of(c("summary", "proposal_distribution"))) %>%
    tidyr::unnest(per_param) %>%
    dplyr::filter(dplyr::if_all(.cols = dplyr::everything(), Negate(is.na))) %>%
    dplyr::filter(dplyr::if_all(.cols = dplyr::everything(), ~ .x > 0))

  tmp2 = fit$waves %>%
    dplyr::select(-dplyr::any_of(c("per_param", "proposal_distribution"))) %>%
    tidyr::unnest(summary) %>%
    dplyr::filter(dplyr::if_all(.cols = dplyr::everything(), Negate(is.na)))

  p1 = ggplot2::ggplot(tmp2, ggplot2::aes(x = wave)) +
    ggplot2::geom_line(ggplot2::aes(y = abs_distance)) +
    ggplot2::geom_point(ggplot2::aes(y = abs_distance)) +
    ggplot2::scale_x_continuous(breaks = ~ floor(min(.x)):ceiling(max(.x))) +
    ggplot2::ylab("Summary distance") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ylim(c(0, NA))

  p2 = ggplot2::ggplot(tmp2, ggplot2::aes(x = wave)) +
    ggplot2::geom_line(ggplot2::aes(y = ESS)) +
    ggplot2::geom_point(ggplot2::aes(y = ESS)) +
    ggplot2::geom_hline(
      yintercept = 200,
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::scale_x_continuous(breaks = ~ floor(min(.x)):ceiling(max(.x))) +
    ggplot2::ylab("Effective sample size") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ylim(c(0, NA))

  p3 = ggplot2::ggplot(tmp, ggplot2::aes(x = wave, colour = param)) +
    ggplot2::geom_line(ggplot2::aes(y = IQR_95_redn)) +
    ggplot2::geom_point(ggplot2::aes(y = IQR_95_redn)) +
    ggplot2::geom_hline(
      yintercept = 0.1,
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(breaks = ~ floor(min(.x)):ceiling(max(.x))) +
    ggplot2::ylab("Reduction in 95% CrI") +
    ggplot2::theme(legend.position = "bottom")

  p4 = ggplot2::ggplot(tmp, ggplot2::aes(x = wave, colour = param)) +
    ggplot2::geom_line(ggplot2::aes(y = rel_mean_change * 100)) +
    ggplot2::geom_point(ggplot2::aes(y = rel_mean_change * 100)) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(breaks = ~ floor(min(.x)):ceiling(max(.x))) +
    ggplot2::ylab("% change in median") +
    ggplot2::theme(legend.position = "bottom")

  patchwork::wrap_plots(
    list(p1, p2, p3, p4),
    tag_level = "new",
    guides = "collect",
    axes = "collect",
    ncol = 2
  ) +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = ggplot2::theme(
        legend.position = 'bottom'
      )
    )
}

#' Plot the evolution of the density function by wave for SMC and adaptive ABC
#'
#' @inheritParams tidyabc_common
#' @inheritParams abc_fit
#'
#' @returns a plot of the density functions by wave
#' @export
#' @concept workflow
#' @examples
#'
#' plot_evolution(
#'   example_adaptive_fit(),
#'   example_truth()
#' )
#'
plot_evolution = function(fit, truth = NULL, ...) {
  p1 = plot(
    fit$summary$density,
    plot_quantiles = FALSE,
    mapping = ggplot2::aes(fill = group, colour = group)
  ) +
    ggplot2::facet_wrap(~name, scales = "free") +
    ggplot2::scale_color_viridis_d(
      aesthetics = c("fill", "colour"),
      ...,
      name = "wave"
    )

  if (!is.null(truth)) {
    truth = truth[unique(fit$summary$param)]
    p1 = p1 +
      ggplot2::geom_vline(
        data = tibble::enframe(unlist(truth)),
        mapping = ggplot2::aes(xintercept = value),
        linetype = "dashed",
        colour = "grey40"
      )
  }

  return(p1)
}


#' Spaghetti plot of resampled posterior fits
#'
#' @inheritParams tidyabc_common
#' @param method one of `"auto"`, `"count"`, or `"density"`, or a named vector
#'   one for each component of `obsdata`.
#' @param n the number of simulations to plot
#' @inheritDotParams posterior_resample
#'
#' @returns a patchwork of `ggplot`s
#' @export
#' @concept workflow
#'
#' @examples
#'
#' plot_simulations(
#'   example_obsdata(),
#'   example_adaptive_fit(),
#'   example_sim_fn
#' )
#'
plot_simulations = function(
  obsdata,
  fit,
  sim_fn,
  method = "auto",
  n = 200,
  ...
) {
  simdata = suppressWarnings(fit$posteriors$abc_sim)
  if (is.null(simdata)) {
    simdata = posterior_resample(
      fit$posteriors %>%
        dplyr::slice_min(order_by = abc_summary_distance, n = n),
      # dplyr::slice_max(order_by = abc_weight, n = n),
      sim_fn = sim_fn,
      n_resamples = 1,
      max_samples = n,
      ...
    )
    simdata = simdata$abc_sim
  }

  if (
    !is.data.frame(obsdata) &&
      is.list(obsdata) &&
      length(setdiff(names(obsdata), names(method))) > 0
  ) {
    method = rep(method, length(obsdata))
    names(method) = names(obsdata)
  }

  plts = lapply(names(obsdata), function(nm) {
    #nm = "data1"
    obs = obsdata[[nm]]
    sims = lapply(simdata, function(s) s[[nm]])
    if (is.data.frame(obs)) {
      # TODO: automatic plotting of data frame like inputs
    } else {
      m = method[[nm]]
      r = c(min(floor(obs)), max(ceiling(obs)))
      bw = (r[2] - r[1]) %/% 50 + 1
      if (
        m == "count" || (m == "auto" && is.integer(obs) && length(obs) < 100)
      ) {
        # count like data

        obs = dplyr::tibble(
          x = seq_along(obs),
          count = obs
        )

        sims = dplyr::bind_rows(
          lapply(
            seq_along(sims),
            function(i) {
              dplyr::tibble(
                x = seq_along(sims[[i]]),
                count = sims[[i]],
                i = i
              )
            }
          )
        )

        plt = ggplot2::ggplot() +
          ggplot2::geom_line(
            data = sims,
            mapping = ggplot2::aes(x = x, y = count, group = i),
            colour = ggplot2::alpha("red", 15 / (15 + n))
          ) +
          ggplot2::geom_bar(
            data = obs,
            mapping = ggplot2::aes(x = x, y = count),
            stat = "identity",
            colour = "black",
            fill = NA
          ) +
          ggplot2::coord_cartesian(xlim = r) +
          # ggplot2::facet_wrap(~ rlang::inject(nm))+
          ggplot2::xlab(nm)
      } else {
        # density like data
        obs = dplyr::tibble(x = obs)
        sims = dplyr::bind_rows(
          lapply(
            seq_along(sims),
            function(i) {
              dplyr::tibble(x = sims[[i]], i = i)
            }
          )
        )

        plt = ggplot2::ggplot() +
          ggplot2::geom_density(
            data = sims,
            mapping = ggplot2::aes(x = x, group = i),
            colour = ggplot2::alpha("red", 15 / (15 + n))
          ) +
          ggplot2::geom_histogram(
            data = obs,
            mapping = ggplot2::aes(x = x, y = ggplot2::after_stat(density)),
            colour = "black",
            fill = NA,
            binwidth = bw
          ) +
          ggplot2::coord_cartesian(xlim = r) +
          # ggplot2::facet_wrap(~ rlang::inject(nm))+
          ggplot2::xlab(nm)
      }
    }
    plt
  })

  patchwork::wrap_plots(
    plts,
    tag_level = "new",
    guides = "collect",
    axes = "collect"
  ) +
    patchwork::plot_annotation(
      tag_levels = "A"
    )
}


#' A parameter posterior correlation plot
#'
#' @inheritParams tidyabc_common
#' @inheritParams abc_fit
#'
#' @returns a patchwork of ggplots including density and
#'   2d scatters for each combination of posteriors.
#' @export
#' @concept workflow
#' @examples
#'
#' p = plot_correlations(
#'   example_adaptive_fit(),
#'   example_truth()
#' )
#'
#' p & ggplot2::theme(
#'   axis.title.y = ggplot2::element_text(angle=70,vjust=0.5),
#'   axis.title.x = ggplot2::element_text(angle=20,hjust=0.5)
#' )
#'
plot_correlations = function(posteriors_df, truth = NULL) {
  if (inherits(posteriors_df, "abc_fit")) {
    posteriors_df = posteriors_df$posteriors
  }

  df = posteriors_df %>%
    dplyr::rename(wt = abc_weight) %>%
    dplyr::select(-dplyr::starts_with("abc_"))
  cols = setdiff(colnames(df), c("wt"))
  d = length(cols)
  plts = list()
  n = nrow(posteriors_df)

  max_alpha = 0.8 # .p_from_n(n, 100, 10000) #200 / (200 + n)

  for (i in seq_along(cols)) {
    for (j in seq_along(cols)) {
      # i=1; j=2
      nmx = as.symbol(cols[[j]])
      nmy = as.symbol(cols[[i]])
      plt_df = df %>%
        dplyr::transmute(
          x = !!nmx,
          y = !!nmy,
          wt,
          alph = wt / max(wt) * max_alpha
        )

      rx = stats::quantile(plt_df$x, c(0.01, 0.99)) * c(0.95, 1 / 0.95)
      ry = stats::quantile(plt_df$y, c(0.01, 0.99)) * c(0.95, 1 / 0.95)

      plt = ggplot2::ggplot(plt_df) +
        # ggplot2::facet_grid(
        #   rlang::inject(cols[[j]]) ~ rlang::inject(cols[[i]])
        # ) +
        ggplot2::xlab(if (i == d) cols[[j]] else NULL) +
        ggplot2::ylab(if (j == 1) cols[[i]] else NULL)

      if ((i == d || i == 1)) {
        # keep
        plt = plt + .gg_set_X_angle(90)
      } else {
        plt = plt +
          ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank()
          )
      }

      if ((j == d || j == 1)) {
        # keep
      } else {
        plt = plt +
          ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
          )
      }

      if (i > j) {
        # 2D density plots (unweighted)

        plt = plt +
          ggplot2::geom_density2d(ggplot2::aes(x = x, y = y)) +
          ggplot2::coord_cartesian(xlim = rx, ylim = ry)
      } else if (i < j) {
        # what kind of size for each point? and what transparency?
        # if all the points were evenly spread we would want about 20% of the plot
        # to be dark. size of each point is therefore:
        pointwidth = diff(rx) / sqrt(n)
        pointheight = diff(ry) / sqrt(n)

        plt = plt +
          # ggplot2::geom_point(
          #   ggplot2::aes(x = x, y = y, alpha = wt),
          #   size = 0.5
          # ) +
          ggplot2::geom_tile(
            ggplot2::aes(
              x = x,
              y = y,
              width = pointwidth,
              height = pointheight,
              alpha = alph
            ),
            linewidth = 0
          ) +
          ggplot2::guides(alpha = ggplot2::guide_none()) +
          ggplot2::coord_cartesian(xlim = rx, ylim = ry) +
          ggplot2::scale_y_continuous(position = "right") +
          ggplot2::scale_x_continuous(position = "top") +
          ggplot2::scale_alpha_identity()
      } else {
        plt = plt +
          ggplot2::geom_density(ggplot2::aes(x = x, weight = wt)) +
          ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.text.y.left = ggplot2::element_blank(),
            axis.text.y.right = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
          ) +
          ggplot2::coord_cartesian(xlim = rx)

        if (i == 1) {
          plt = plt +
            ggplot2::scale_x_continuous(position = "top")
        }
      }

      if (!is.null(truth)) {
        xt = truth[[cols[[j]]]]
        yt = truth[[cols[[i]]]]
        if (i == j) {
          plt = plt +
            ggplot2::geom_vline(
              xintercept = yt,
              linetype = "dashed",
              colour = "grey40"
            )
        } else {
          plt = plt +
            ggplot2::geom_point(
              data = dplyr::tibble(
                x = xt,
                y = yt
              ),
              mapping = ggplot2::aes(x = x, y = y),
              colour = "red",
              inherit.aes = FALSE,
              size = 1
            )
        }
      }

      plts = c(plts, plt)
    }
  }
  patchwork::wrap_plots(
    plts,
    tag_level = "new",
    guides = "collect",
    # axes = "collect",
    ncol = length(cols),
    nrow = length(cols)
  )
}
