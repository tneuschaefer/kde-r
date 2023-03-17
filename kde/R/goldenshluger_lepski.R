#' Goldenshluger Lepski Method
#'
#' @description \code{goldenshluger_lepski} estimates an optimal bandwidth
#'   for kernel density estimation using the Goldenshluger Lepski Method
#'
#' @param x numeric vector of the observation sample
#' @param kernel kernel function used for kernel density estimation
#' @param n number of bandwidths to be optimized from. \code{goldenshluger_lepski}
#'    selects a bandwidth contained in  (1/n, 2/n, ..., 1)
#' @param lambda positive scalar used as tuning parameter
#' @param N number of subdivisions used in discretization of integrals
#' @param built_in choose one of the built-in kernels instead of providing one yourself
#' @param na.rm logical; if TRUE, missing values will be removed from x
#'
#' @details The Goldenshluger and Lepski method tries to minimize an upper bound
#'   for the mean integrated squared error (MISE) of a kernel density estimator
#'
#'   The bias/variance decomposition is used. Because the bias term depends on
#'   the unknown density, a double kernel approach is used to to estimate the
#'   bias.
#'
#'   An upper bound for the variance is estimated, this estemator depends on a
#'   tuning paramter \code{lambda}. Recommended is \code{lambda}=1.2
#'
#' @return The estimatet optimal bandwidth
#'
#' @seealso \code{\link{kernel_estimator}} for more information about kernel
#'   density estimation
#'
#'   \code{\link{cross_validation}} and \code{\link{pco_method}} for other
#'   automatic bandwidth-selection algorithms
#'
#' @source Comte, F.: Nonparametric Esimation. Spartacus-Idh (2017)
#'
#' @examples
#'
#' x <- stats::rnorm(10)
#' h <- goldenshluger_lepski(x, kernel = stats::dnorm)
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
#' @export
goldenshluger_lepski <- function(x, kernel = stats::dnorm, n = 40, lambda = 1.2, N = 100L,
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
    "lambda must be greater than 1" = lambda > 1,
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
      epanechnikov = epanechnikov,
      rectangular = rectangular,
      triangular = triangular,
      biweight = biweight,
      silverman = silverman
    )
  }

  a <- min(x)
  b <- max(x)
  ab <- seq(a, b, length.out = N)

  m <- length(x)

  K_Norm1_sqr <- (sum(abs(kernel(ab))) * (b - a) / N)^2
  K_Norm2_sqr <- (sum(kernel(ab)^2)) * (b - a) / N
  Norm2 <- matrix(0, n, n)
  V <- rep(0, n)
  A <- rep(0, n)

  for (k in 1:n) {
    h <- k / n
    V[k] <- K_Norm1_sqr * K_Norm2_sqr / h / m
  }

  for (i in 1:n) {
    h <- k / n
    f_h <- kernel_estimator(x, kernel = kernel, bandwidth = h)
    for (j in 1:n) {
      h2 <- j / n
      f_h2 <- kernel_estimator(x, kernel = kernel, bandwidth = h2)
      K_h2 <- function(x) kernel(x / h2) / h2
      f_h_h2 <- function(x) {
        res <- rep(0, length(x))
        for (k in 1:length(x)) {
          res[k] <- sum(K_h2(ab) * f_h(x[k] - ab)) * (b - a) / N
        }
        res
      }
      Norm2[i, j] <- sum((f_h_h2(ab) - f_h2(ab))^2) * (b - a) / N
    }
    A[i] <- max(Norm2[i, ] - lambda * V)
  }

  which.min(A + 2 * lambda * V) / N
}
