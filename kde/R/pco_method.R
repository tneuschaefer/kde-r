#' Penalized Comparison to Overfitting
#' 
#' @description \code{pco_method} estimates an optimal bandwidth
#'   for kernel density estimator using the PCO method
#' 
#' @param X numeric vector of the observation sample
#' @param K kernel function used for kernel density estimation
#' @param n number of bandwidths to be optimized from. \code{pco_method} selects
#'   a bandwidth contained in  (1/n, 2/n, ..., 1)  
#' @param lambda positive scalar used as tuning parameter
#' @param N number of subdivisions used in discretization of integrals 
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
#'   \code{\link{cross_validation}} and for other automatic bandwidth-selection
#'   algorithms 
#'   
#' @source Lacour et al, Estimator selection: a new method with applications to
#'   kernel density estimation (2017), https://arxiv.org/abs/1607.05091
#' 
#' 
#' @examples 
#' 
#'   X <- rnorm(1000)
#'   h <- pco_method(X, dnorm)
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
pco_method <- function(X, K, n = 40, lambda = 1, N = 100L) {

  #Sample condition
  stopifnot("X must be a numeric" = is.numeric(X),
            "X must be a non empty" = length(X) > 0)

  #Kernel condition
  stopifnot("K is not a function" = is.function(K))

  #Bandwidth condition
  stopifnot("n must be numeric" = is.numeric(n),
            "n must have length one" = length(n) == 1)

  #Lamda condition
  stopifnot("Lambda must be numeric" = is.numeric(lambda),
            "Lambda must have length one" = length(lambda) == 1,
            "Lambda must be > 0" = lambda > 0)

  #Subdivisions condition
  stopifnot("N must be numeric" = is.numeric(N),
            "N must have length one" = length(N) == 1)

  n <- as.integer(n)
  N <- as.integer(N)
  m <- length(X)
  a <- min(X)
  b <- max(X)
  ab <- seq(a, b, N)

  risk <- c()
  h_min <- 1/n
  f_h_min <- kernel_estimator(X, K , h_min)
  K_h_min <- function(x) K(x/h)/h

  for(k in 1:n) {
    h <- k/n
    f_h <- kernel_estimator(X, K, h)

    bias_estim <- sum( (f_h_min(ab)-f_h(ab))^2 )*(b-a)/N

    K_h <- function(x) K(x/h)/h
    K_h_Norm <- sum(K_h(ab)^2)*(b-a)/N

    var_h <- lambda*K_h_Norm/m

    z <- sum( (K_h_min(ab)-K_h(ab))^2 )*(b-a)/N
    bias_h <- z/m

    penality <- var_h - bias_h

    l_pco <- bias_estim + penality
    risk <- c(risk, l_pco)
  }
  risk <- risk[-1]
  k <- which.min(risk)
  k/n
}