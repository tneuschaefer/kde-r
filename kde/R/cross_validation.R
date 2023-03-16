#' Cross-Validation
#'
#' @description \code{cross_validate} estimates an optimal bandwidth
#'   for kernel density estimation using the cross-validation method
#'
#' @param X numeric vector of the observation sample
#' @param K kernel function used for the kernel density estimation
#' @param n number of bandwidths to be optimized from. \code{cross_validate} 
#'   selects a bandwidth contained in  (1/n, 2/n, ..., 1)
#' @param N number of subdivisions used in discretization of integrals
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
#'   X <- rnorm(1000)
#'   h <- cross_validate(X, dnorm)
#'   f <- kernel_estimator(X, dnorm, h)
#'   
#'   a <- min(X)
#'   b <- max(X)
#'   ab <- seq(a, b, length.out = 100)
#'   
#'   plot(ab, dnorm(ab), type = "l")
#'   line(ab, f(ab), col = "red")
#'   legend("topleft", legend = c("true", "estimated h"), col = c("black", "red"),
#'          pch = "|")
#' 
#' @include kernel_estimator.R
#' 
#' @export
cross_validate <- function(X, K, n = 40L, N = 100) {

  #Sample condition
  stopifnot("X must be a numeric" = is.numeric(X),
            "X must be a non empty" = length(X) > 0)

  #Kernel condition
  stopifnot("K is not a function" = is.function(K))

  #Bandwidth condition
  stopifnot("n must be numeric" = is.numeric(n),
            "n must have length one" = length(n) == 1)

  #Subdivision condition
  stopifnot("N must be an integer" = is.integer(n),
            "N must have length one" = length(n) == 1)

  n <- as.integer(n)
  a <- min(X)
  b <- max(X)
  ab <- seq(a, b, N)

  m <- length(X)
  estimf <- matrix(0, n, N)
  Mat <- array(0, c(m, m, n))
  crit <- rep(0, n)

  for(k in 1:n) {
    h <- k/n
    f <- kernel_estimator(X, K, h)
    estimf[k,] <- f(ab)
  }

  for(k in 1:n) {
    h <- k/n
    for(i in 1:m) {
      for(j in 1:m) {
        if(i != j) Mat[i,j,k] <- K((X[i]-X[j])/h)/h
      }
    }
  }

  for(k in 1:n) {
    crit[k] <- sum(estimf[k,]^2)*(b-a)/N - 2*sum(Mat[,,k])/n/(n-1)
  }
  crit <- crit[-1]
  which.min(crit)/n
}
