#' Generate a log of a kernel function for ABC weight calculation
#'
#' \deqn{
#' \text{Uniform: } K(d, \epsilon) = \frac{1}{2\epsilon} \mathbb{I}(d \leq \epsilon)
#' }
#' \deqn{
#' \text{Triangular: } K(d, \epsilon) = \frac{2}{\epsilon}\left(1 - \frac{d}{\epsilon}\right) \mathbb{I}(d \leq \epsilon)
#' }
#' \deqn{
#' \text{Epanechnikov: } K(d, \epsilon) = \frac{3}{4\epsilon}\left(1 - \frac{d^2}{\epsilon^2}\right) \mathbb{I}(d \leq \epsilon)
#' }
#' \deqn{
#' \text{Biweight: } K(d, \epsilon) = \frac{15}{16\epsilon}\left(1 - \frac{d^2}{\epsilon^2}\right)^2 \mathbb{I}(d \leq \epsilon)
#' }
#' \deqn{
#' \text{Gaussian: } K(d, h) = \frac{1}{\sqrt{2\pi h^2}} \exp\left(-\frac{d^2}{2h^2}\right)
#' }
#' where \eqn{d} is the distance between observed and simulated summary
#' statistics, \eqn{\epsilon} is the tolerance level, \eqn{h} is the bandwidth
#' parameter for the Gaussian kernel, and \eqn{\mathbb{I}} is the indicator
#' function taking the value 1 if \eqn{d \leq \epsilon} and 0 otherwise (this
#' applies to all kernels except Gaussian).
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
