#' Construct a Kernel Density Estimator
#'
#' @description \code{kernel_estimator} creates a kernel density estimator using
#'   the provided sample, kernel and bandwidth
#'
#' @param X numeric vector of the observation sample
#' @param K kernel function
#' @param h non-negative numeric scalar, the bandwidth of the estimator
#'
#' @return A function which estimates the probability density function of the
#'   sample \code{X}
#'
#' @details The function given by \code{kernel_estimator} is vectorized
#'
#' @seealso \code{\link{cross_validation}} and \code{\link{pco_method}} for
#' automatic bandwidth selection methods
#'
#' @source Comte, F.: Nonparametric Esimation. Spartacus-Idh (2017)
#'
#' @examples
#' X <- rnorm(100)
#'
#' f1 <- kernel_estimator(X, dnorm, 1)
#' f2 <- kernel_estimator(X, dnorm, 0.5)
#' f3 <- kernel_estimator(X, dnorm, 0.1)
#'
#' a <- min(X)
#' b <- max(X)
#' ab <- seq(a, b, length.out = 100)
#'
#' plot(ab, dnorm(ab), type = "l", xlab = "x", ylab = "density")
#' lines(ab, f1(ab), col = "red")
#' lines(ab, f2(ab), col = "blue")
#' lines(ab, f3(ab), col = "green")
#'
#' legend("topleft",
#'   col = c("black", "red", "blue", "green"),
#'   legend = c("true", "h = 1", "h = 0.5", "h = 0.1"), pch = "l"
#' )
#'
#' @export
kernel_estimator <- function(X, K, h = 1) {
  # Sample condition
  stopifnot(
    "X must be a numeric" = is.numeric(X),
    "X must be non empty" = length(X) > 0
  )

  # Kernel condition
  stopifnot("K is not a funtion" = is.function(K))

  # Bandwidth condition
  stopifnot(
    "h must be numeric" = is.numeric(h),
    "h must be non-negativ" = h > 0,
    "h must have length one" = length(h) == 1
  )

  estimator <- function(x) {
    stopifnot(
      "x must be numeric" = is.numeric(x),
      "x miust be non empty" = length(x) > 0
    )

    n <- length(X)
    m <- length(x)
    Mat <- matrix(0, nrow = n, ncol = m)

    for (i in 1:n) {
      for (j in 1:m) {
        U <- (X[i] - x[j]) / h
        Mat[i, j] <- K(U) / h
      }
    }
    apply(Mat, 2, mean)
  }
}
