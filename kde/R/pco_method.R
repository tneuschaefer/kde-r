#' Penalized Comparison to Overfitting
#'
#' @description \code{pco_method} estimates an optimal bandwidth
#'   for kernel density estimator using the PCO method
#'
#' @param x numeric vector of the observation sample
#' @param kernel kernel function used for kernel density estimation
#' @param n number of bandwidths to be optimized from. \code{pco_method} selects
#'   a bandwidth contained in  (1/n, 2/n, ..., 1)
#' @param lambda positive scalar used as tuning parameter
#' @param N number of subdivisions used in discretization of integrals
#' @param built_in choose one of the built-in kernels instead of providing one yourself
#' @param na.rm logical; if TRUE, missing values will be removed from x
#'
#' @details The PCO method tries to minimize an upper bound for the mean
#'   integrated squared error (MISE) of a kernel density estimator
#'
#'   \code{pco_method} calculates the PCO criterion value for every bandwidth,
#'   approximating the risk and selecting the bandwidth with the minimal risk
#'
#'   The bias/variance decomposition is used. Because the bias term depends on
#'   the unknown density, a comparison of the estimator with an associated
#'   bandwidth is used to estimate the bias term. Here the estimator with the
#'   smallest bandwidth is used, this estimator is surely overfitting
#'
#'   Additional a penalty term is computed as the sum of the risk discomposition
#'   variance and the variance of the bias term estimation. In the calculation
#'   the tuning paramter \code{lambda} is used. Recommended is \code{lambda}=1
#'
#'   The PCO criterion is given as the sum of the comparison to the overfitting
#'   and the penalty term. The PCO method tries to find a balance between this
#'   two terms, therefor the name \code{penalized comparison to overfitting}
#'
#' @return The estimated optimal bandwidth
#'
#' @seealso \code{\link{kernel_estimator}} for more information about kernel
#'   density estimation
#'
#'   \code{\link{cross_validation}} and \code{\link{goldenshluger_lepski}} for
#'   other automatic bandwidth-selection algorithms
#'
#' @source Lacour et al, Estimator selection: a new method with applications to
#'   kernel density estimation (2017), https://arxiv.org/abs/1607.05091
#'
#' @examples
#'
#' x <- stats::rnorm(100)
#' h <- pco_method(x, kernel = stats::dnorm)
#' f <- kernel_estimator(x, kernel = stats::dnorm, bandwidth = h)
#'
#' a <- min(x)
#' b <- max(x)
#' ab <- seq(a, b, length.out = 100)
#'
#' plot(ab, stats::dnorm(ab), type = "l")
#' lines(ab, f(ab), col = "red")
#' legend("topleft",
#'   legend = c("true", "estimated h"), col = c("black", "red"),
#'   pch = "|"
#' )
#'
#' @include kernel_estimator.R
#'
#' @export
pco_method <- function(x, kernel = stats::dnorm, n = 40, lambda = 1, N = 100L,
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
    "lambda must be numeric" = is.numeric(lambda),
    "lambda must be greater than 0" = lambda > 0,
    "lambda must not be empty" = length(lambda) > 0,
    "N must be numeric" = is.numeric(N),
    "N must not be empty" = length(N) > 0
  )

  # reformat arguments where possible without throwing an error
  n <- as.integer(abs(n[1]))
  lambda <- lambda[1]
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

  m <- length(x)
  a <- min(x)
  b <- max(x)
  ab <- seq(a, b, N)

  risk <- c()
  h_min <- 1 / n
  f_h_min <- kernel_estimator(x, kernel, h_min)
  K_h_min <- function(x) kernel(x / h) / h

  for (k in 1:n) {
    h <- k / n
    f_h <- kernel_estimator(x, kernel, h)

    bias_estim <- sum((f_h_min(ab) - f_h(ab))^2) * (b - a) / N

    K_h <- function(x) kernel(x / h) / h
    K_h_Norm <- sum(K_h(ab)^2) * (b - a) / N

    var_h <- lambda * K_h_Norm / m

    z <- sum((K_h_min(ab) - K_h(ab))^2) * (b - a) / N
    bias_h <- z / m

    penality <- var_h - bias_h

    l_pco <- bias_estim + penality
    risk <- c(risk, l_pco)
  }
  risk <- risk[-1]
  k <- which.min(risk)
  k / n
}
