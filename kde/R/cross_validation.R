library(stats)

cross_validate <- function(X, K = dnorm, n = 40L, N = 100) {

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
        if(i != j) Mat[i,j,k] <- K((X[i]+X[j])/h)/h
      }
    }
  }

  for(k in 1:n) {
    crit[k] <- sum(estimf[k,]^2)*(b-a)/N - 2*sum(Mat[,,k])/n/(n-1)
  }
  which.min(crit)/n
}

#Beispiel
set.seed(10)
n <- 1000
X <- rnorm(n)
h <- cross_validate(X)

hnorm_estim <- kernel_estimator(X, dnorm, h)
norm_estim(10)
dnorm(10)
norm_estim(c(0,1,-1))
dnorm(c(0,1,-1))

#Plotten
a <- min(X)
b <- max(X)

x_plot <- seq(a, b, length.out = 100)
x_estim <- norm_estim(x_plot)
x_real <- dnorm(x_plot)

plot(x_plot, x_real, type = "l", col = "red", lwd = 2, xlab = "X", ylab = "Density")
lines(x_plot, x_estim, col = "blue")
