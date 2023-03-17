#' Construct a Kernel Density Estimator
#'
#' @description \code{kernel_estimator} creates a kernel density estimator using
#'   the provided sample, kernel and bandwidth
#'
#' @param x numeric vector of the observation sample
#' @param kernel kernel function, default is Gaussian kernel
#' @param bandwidth non-negative numeric scalar, the bandwidth of the estimator
#'  default is 1
#' @param built_in choose one of the built-in kernels instead of providing one yourself
#' @param na.rm logical; if TRUE, missing values will be removed from x
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
#' x <- rnorm(100)
#'
#' f1 <- kernel_estimator(x, kernel = dnorm, bandwidth = 1)
#' f2 <- kernel_estimator(x, kernel = dnorm, bandwidth = 0.5)
#' f3 <- kernel_estimator(x, kernel = dnorm, bandwidth = 0.1)
#'
#' a <- min(x)
#' b <- max(x)
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
kernel_estimator <- function(x, kernel = stats::dnorm,
                             bandwidth = 1,
                             built_in = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "silverman"),
                             na.rm = FALSE) {
  # remove NA values if na.rm is set to TRUE
  if (na.rm) x <- x[!is.na(x)]

  # ensuring requirements
  stopifnot(
    "x must be a numeric" = is.numeric(x),
    "x must not be empty" = length(x) > 0,
    "kernel must be a function" = is.function(kernel),
    "bandwidth must be numeric" = is.numeric(bandwidth),
    "bandwidth must be greater than 0" = bandwidth > 0,
    "x contains missing values" = !anyNA(x)
  )

  # in case the provided bandwidth is of length greater than 1
  bandwidth <- bandwidth[1]

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

  # returning estimator function
  function(t) {
    stopifnot(
      "t must be numeric" = is.numeric(t),
      "t must not be empty" = length(t) > 0
    )

    n <- length(x)

    # calculate the value of the estimated probability density for each input value in t
    sapply(t, function(t_i) {
      1 / (n * bandwidth) * sum(kernel((x - t_i) / bandwidth))
    })
  }
}

#' Indicator function
#'
#' @param u real number
#'
#' @return 1 if -1 <= u <= 1, 0 otherwise
ind <- function(u) {
  return(ifelse(abs(u) <= 1, 1, 0))
}


#' Epanechnikov kernel
#'
#' @param u real number
#'
#' @return value of Epanechnikov kernel for u
#' @export
epanechnikov <- function(u) {
  return(0.75 * (1 - u^2) * ind(u))
}

#' Rectangular Kernel
#'
#' @param u real number
#'
#' @return value of Rectangular kernel for u
#' @export
rectangular <- function(u) {
  return(0.5 * ind(u))
}


#' Triangular Kernel
#'
#' @param u real number
#'
#' @return value of Triangular kernel for u
#' @export
triangular <- function(u) {
  return((1 - abs(u)) * ind(u))
}

#' Biweight Kernel
#'
#' @param u real number
#'
#' @return value of Biweight kernel for u
#' @export
biweight <- function(u) {
  return(0.9375 * (1 - u^2)^2 * ind(u))
}

#' Silverman
#'
#' @param u real number
#'
#' @return value of Silverman kernel for u
#' @export
silverman <- function(u) {
  return(0.5 * exp(-abs(u) / sqrt(2)) * sin(abs(u) / sqrt(2) + pi / 4))
}
