#' Construct a Kernel Density Estimator
#'
#' @description \code{kernel_estimator} creates a kernel density estimator using
#'   the provided sample, kernel and bandwidth
#'
#' @param x numeric vector of the observation sample
#' @param kernel kernel function, default is Gaussian kernel
#' @param bandwith non-negative numeric scalar, the bandwidth of the estimator
#'  default is 1
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
kernel_estimator <- function(x, kernel = stats::dnorm, bandwith = 1) {
  # Sample condition
  stopifnot(
    "X must be a numeric" = is.numeric(x),
    "X must be non empty" = length(x) > 0
  )

  # Kernel condition
  stopifnot("K is not a funtion" = is.function(kernel))

  # Bandwidth condition
  stopifnot(
    "h must be numeric" = is.numeric(bandwith),
    "h must be non-negativ" = bandwith > 0,
    "h must have length one" = length(bandwith) == 1
  )

  # returning estimator function
  function(t) {
    stopifnot(
      "x must be numeric" = is.numeric(t),
      "x miust be non empty" = length(t) > 0
    )

    n <- length(x)
    m <- length(t)
    Mat <- matrix(0, nrow = n, ncol = m)

    for (i in 1:n) {
      for (j in 1:m) {
        U <- (x[i] - t[j]) / bandwith
        Mat[i, j] <- kernel(U) / bandwith
      }
    }
    apply(Mat, 2, mean)
  }
}
