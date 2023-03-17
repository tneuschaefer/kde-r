#' Cross-Validation
#'
#' @description \code{cross_validation} estimates an optimal bandwidth
#'   for kernel density estimation using the cross-validation method
#'
#' @param x numeric vector of the observation sample
#' @param K kernel function used for the kernel density estimation
#' @param n number of bandwidths to be optimized from. \code{cross_validation}
#'   selects a bandwidth contained in  (1/n, 2/n, ..., 1)
#' @param N number of subdivisions used in discretization of integrals
#' @param built_in choose one of the built-in kernels instead of providing one yourself
#' @param na.rm logical; if TRUE, missing values will be removed from x
#'
#' @details The cross-validation method tries to minimize the mean squared
#'   integrated error (MISE) of a kernel density estimator
#'
#'   \code{cross_validation} approximates the estimator-dependent part of the
#'   risk and selects the bandwidth with the minimal risk
#'
#' @return The estimatet optimal bandwidth
#'
#' @seealso \code{\link{kernel_estimator}} for more information about kernel
#'   density estimation
#'
#'   \code{\link{cross_validation}} and for other automatic bandwidth-selection
#'   algorithms
#'
#' @source Comte, F.: Nonparametric Esimation. Spartacus-Idh (2017)
#'
#' @examples
#'
#' x <- stats::rnorm(1000)
#' h <- cross_validation(x, kernel = stats::dnorm)
#' f <- kernel_estimator(x, kernel = stats::dnorm, bandwidth = h)
#'
#' a <- min(x)
#' b <- max(x)
#' ab <- seq(a, b, length.out = 100)
#'
#' plot(ab, stats::dnorm(ab), type = "l")
#' line(ab, f(ab), col = "red")
#' legend("topleft",
#'   legend = c("true", "estimated h"), col = c("black", "red"),
#'   pch = "|"
#' )
#'
#' @include kernel_estimator.R
#'
#' @export
cross_validation <- function(x, kernel = stats::dnorm, n = 40L, N = 100L, 
built_in = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "silverman"), na.rm = FALSE) {
  # remove NA values if na.rm is set to TRUE
  if (na.rm) x <- x[!is.na(x)]

  # ensuring requirements
  stopifnot(
    "x must be numeric" = is.numeric(x),
    "x must not be empty" = length(x) > 0,
    "x contains missing values" = !anyNA(x),
    "kernel must be a function" = is.function(kernel),
    "n must be numeric" = is.numeric(n),
    "n must not be empty" = length(n) > 0,
    "N must be numeric" = is.numeric(N),
    "N must not be empty" = length(N) > 0
  )

  # reformat arguments where possible without throwing an error
  n <- as.integer(abs(n[1]))
  N <- as.integer(abs(N[1]))

  # if built_in was provided use that as the kernel
  if (!missing(built_in)) {
    # match the argument for built_in
    built_in <- match.arg(built_in)
    kernel <- switch(built_in,
      gaussian = stats::dnorm,
      epanechnikov = epanechnikov(),
      rectangular = rectangular(),
      triangular = triangular(),
      biweight = biweight(),
      silverman = silverman()
    )
  }

  a <- min(x)
  b <- max(x)
  ab <- seq(a, b, N)

  m <- length(x)
  estimf <- matrix(0, n, N)
  Mat <- array(0, c(m, m, n))
  crit <- rep(0, n)

  for (k in 1:n) {
    h <- k / n
    f <- kernel_estimator(x, kernel, h)
    estimf[k, ] <- f(ab)
  }

  for (k in 1:n) {
    h <- k / n
    for (i in 1:m) {
      for (j in 1:m) {
        if (i != j) Mat[i, j, k] <- kernel((x[i] - x[j]) / h) / h
      }
    }
  }

  for (k in 1:n) {
    crit[k] <- sum(estimf[k, ]^2) * (b - a) / N - 2 * sum(Mat[, , k]) / n / (n - 1)
  }
  crit <- crit[-1]
  which.min(crit) / n
}
