#' Title
#'
#' @param X
#' @param K
#' @param n
#' @param kappa
#' @param N
#'
#' @return
#' @export
#'
#' @examples
goldenshluger_lepski <- function(x, kernel = stats::dnorm, n = 40, kappa = 1.2, N = 100L) {
  # Sample condition
  stopifnot(
    "X must be a numeric" = is.numeric(x),
    "X must be a non empty" = length(x) > 0
  )

  # Kernel condition
  stopifnot("K is not a function" = is.function(kernel))

  # Bandwidth condition
  stopifnot(
    "n must be numeric" = is.numeric(n),
    "n must have length one" = length(n) == 1
  )

  # Kappa condition
  stopifnot(
    "Kappa must be numeric" = is.numeric(kappa),
    "Kappa must have length one" = length(kappa) == 1,
    "Kappa must be >= 1" = kappa >= 1
  )

  # Subdivisions condition
  stopifnot(
    "N must be an integer" = is.integer(N),
    "N must have length one" = length(N) == 1
  )

  n <- as.integer(n)
  a <- min(x)
  b <- max(x)
  ab <- seq(a, b, length.out = N)

  m <- length(x)
  estimf <- matrix(0, n, N)

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
    f_h <- kernel_estimator(x, kernel, h)
    for (j in 1:n) {
      h2 <- j / n
      f_h2 <- kernel_estimator(x, kernel, h2)
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
    A[i] <- max(Norm2[i, ] - kappa * V)
  }

  which.min(A + 2 * kappa * V) / N
}
