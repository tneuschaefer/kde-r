library(stats)

cross_validate <- function(X, K, n = 40L) {
  stopifnot("X invalid data" = is.numeric(X) && length(X) > 0,
            "K is no funtion" = is.function(K),
            "n is no integer" = is.integer(n))
  
  a <- min(X)
  b <- max(X)
  L <- 100
  ab <- seq(a, b, length.out = L)
  
  n <- as.integer(n)
  m <- length(X)
  Mat <- array(0, c(m, m, n))
  intf <- rep(0, n)
  crit <- rep(0, n)
  
  #integrate using R intern function
  for(k in 1:n) {
    h <- k/n
    f <- kernel_estimator(X, K, h)
    g <- function(x) f(x)^2
    intf[k] <- integrate(g, lower = -Inf, upper = Inf, subdivisions = 2000)$value
  }
  
  l <- rep(1, m)
  R <- diag(l, m, m)
  for(k in 1:n) {
    h <- k/n
    for(i in 1:m) {
      for(j in 1:m) {
        if(i != j) {
          U <- (X[i]-X[j])/h
          Mat[i,j,k] <- K(U)/h
        }
      }
    }
  }
  
  for(k in 1:n) crit[k] <- intf[k] - 2*sum(Mat[,,k])/(n*(n-1))
  k <- which.min(crit)
  k/n
}

#Beispiel
set.seed(10)
n <- 1000
X <- rexp(n)
h <- cross_validate(X, dnorm)

exp_estim <- kernel_estimator(X, dnorm, h)
exp_estim(10)
dexp(10)
exp_estim(c(0,1,-1))
dexp(c(0,1,-1))

#Plotten
a <- min(X)
b <- max(X)

x_plot <- seq(a, b, length.out = 100)
x_estim <- exp_estim(x_plot)
x_real <- dexp(x_plot)

plot(x_plot, x_real, type = "l", col = "red", lwd = 2, xlab = "X", ylab = "Density")
lines(x_plot, x_estim, col = "blue")
