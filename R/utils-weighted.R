# # QWEN3 generated
# # need to see what it has done
#
# # Actually, let's use a cleaner PAVA implementation
# weighted_isotonic_regression <- function(y, weights = NULL) {
#   if (is.null(weights)) weights <- rep(1, length(y))
#
#   n <- length(y)
#   if (n <= 1) return(y)
#
#   # PAVA algorithm with weights
#   # Returns isotonic (non-decreasing) regression
#   values <- y
#   weights_vec <- weights
#   n_vec <- rep(1, n)  # number of points in each block (initially 1)
#
#   i <- 1
#   while (i < length(values)) {
#     if (values[i] > values[i + 1]) {
#       # Pool adjacent violators
#       w1 <- weights_vec[i]
#       w2 <- weights_vec[i + 1]
#       pooled_value <- (w1 * values[i] + w2 * values[i + 1]) / (w1 + w2)
#
#       # Update vectors
#       values <- c(values[1:(i-1)], pooled_value, values[(i+2):length(values)])
#       weights_vec <- c(weights_vec[1:(i-1)], w1 + w2, weights_vec[(i+2):length(weights_vec)])
#       n_vec <- c(n_vec[1:(i-1)], n_vec[i] + n_vec[i + 1], n_vec[(i+2):length(n_vec)])
#
#       # Step back to check for new violations
#       if (i > 1) i <- i - 1
#     } else {
#       i <- i + 1
#     }
#   }
#
#   # Expand back to original length
#   result <- rep(values, times = n_vec)
#   return(result)
# }
#
# # Example usage:
# # x <- sort(runif(20, 0, 10))
# # y <- pnorm(x, mean = 5, sd = 2) + rnorm(20, 0, 0.05)  # noisy CDF-like
# # w <- runif(20, 0.5, 2)  # random weights
# # y_smooth <- weighted_monotonic_loess(x, y, weights = w, span = 0.4)
# # plot(x, y, col = "red", pch = 16)
# # lines(sort(x), y_smooth[order(x)], col = "blue", lwd = 2)
#
#
# # More efficient version using splines
# weighted_monotonic_smoother <- function(x, y, weights = NULL, n_knots = 10) {
#   if (is.null(weights)) weights <- rep(1, length(x))
#
#   # Use I-splines (monotonic by construction when coefficients ≥ 0)
#   # For this, we'd need to implement non-negative least squares
#   # Or use existing constrained optimization
#
#   # For now, use the PAVA approach but more efficiently
#   ord <- order(x)
#   x_ord <- x[ord]
#   y_ord <- y[ord]
#   w_ord <- weights[ord]
#
#   # Weighted LOESS at original points (more expensive but more accurate)
#   y_smooth <- numeric(length(x))
#   for (i in seq_along(x_ord)) {
#     y_smooth[i] <- weighted_local_fit(x_ord[i], x_ord, y_ord, w_ord, 0.3, 1)
#   }
#
#   # Apply isotonic regression
#   y_mono <- weighted_isotonic_regression(y_smooth, w_ord)
#
#   # Return in original order
#   result <- numeric(length(x))
#   result[ord] <- y_mono
#   return(result)
# }
#
#
# # Function to fit weighted local polynomial at a point
# weighted_local_fit <- function(x_target, x_data, y_data, w_data, span, degree) {
#   # Find neighbors within span (fraction of data)
#   n <- length(x_data)
#   k <- max(1, floor(span * n))
#
#   # Get k nearest neighbors
#   dists <- abs(x_data - x_target)
#   # Use order to get k closest points
#   ord_dists <- order(dists)
#   close_idx <- ord_dists[seq_len(min(k, n))]
#
#   x_nbrs <- x_data[close_idx]
#   y_nbrs <- y_data[close_idx]
#   w_nbrs <- w_data[close_idx]
#
#   # Compute kernel weights (tricube kernel)
#   max_dist <- max(dists[close_idx])
#   if (max_dist == 0) max_dist <- 1e-10  # Avoid division by zero
#   u <- dists[close_idx] / max_dist
#   kernel_w <- (1 - u^3)^3
#
#   # Combined weights
#   combined_w <- w_nbrs * kernel_w
#
#   # Fit polynomial (degree 0, 1, or 2)
#   if (degree == 0) {
#     # Local weighted average
#     result <- sum(combined_w * y_nbrs) / sum(combined_w)
#   } else {
#     # Local weighted polynomial
#     X <- outer(x_nbrs - x_target, 0:degree, "^")  # [x^0, x^1, ..., x^degree]
#     W <- diag(combined_w)
#     beta <- solve(t(X) %*% W %*% X, t(X) %*% W %*% y_nbrs)
#     result <- beta[1]  # Intercept term (value at x_target)
#   }
#
#   return(result)
# }
