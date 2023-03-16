#' Title
#'
#' @param X
#' @param K
#' @param n
#' @param lambda
#' @param N
#'
#' @return
#' @export
#'
#' @examples
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
  stopifnot("Lamda must be numeric" = is.numeric(kappa),
            "Lamda must have length one" = length(kappa) == 1,
            "Lamda must be > 0" = lambda > 0)

  #Subdivisions condition
  stopifnot("N must be an integer" = is.integer(N),
            "N must have length one" = length(N) == 1)

  n <- as.integer(n)
  m <- length(X)
  a <- min(X)
  b <- max(X)
  ab <- seq(a, b, N)

  risk <- c()
  h_min <- 1/n
  f_h_min <- kernel_estimator(X, K , h_min)
  K_h_min <- function(x) K(x/h)/h

  for(k in 1:n) {
    f_h <- kernel_estimator(X, K, h)

    bias_estim <- sum( (f_h_min(ab)-f_h(ab))^2 )*(b-a)/N

    K_h <- function(x) K(x/h)/h
    K_h_Norm <- sum(K_h(ab)^2)*(b-a)/N

    var_h <- lamda*K_h_Norm/m

    z <- sum( (K_h_min(ab)-K_h(ab))^2 )*(b-a)/N
    bias <- z/m

    penality <- var_h - bias_h

    l_pco <- bias_estim + penality
    risk <- c(risk, l_pco)
  }
  k <- which.min(risk)
  k/n
}
