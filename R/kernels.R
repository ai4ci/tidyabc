#' Generate a non normalised log of kernel for ABC weight calculation
#'
#' \deqn{
#' \alpha_{i,j} = \frac{1}{\beta+1}\\
#' \beta = \alpha +1
#' }
#'
#' @param u the distance (dissimilarilty) between
#' @inheritParams common_internal
#' @param kernel one of "epanechnikov", "uniform", "triangular", "biweight",
#'   "gaussian"
#'
#' @returns the log of the kernel function
#' @keywords internal
.log_kernel = function(
  u,
  epsilon,
  kernel = c(
    c("epanechnikov", "uniform", "triangular", "biweight", "gaussian")
  )
) {
  kernel = match.arg(kernel)

  u = abs(u / epsilon)

  tmp = suppressWarnings(switch(
    kernel,
    "epanechnikov" = ifelse(u > 1, -Inf, log(1 - u^2) + log(3 / 4)),
    "uniform" = ifelse(u > 1, -Inf, 0.5),
    "triangular" = ifelse(u > 1, -Inf, log(1 - u)),
    "biweight" = ifelse(u > 1, -Inf, 3 * log(1 - u^2) + log(15 / 16)),
    "gaussian" = -0.5 * u^2 - log(sqrt(2 * pi))
  ))

  tmp
}
