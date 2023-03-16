## cross validation
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

## goldenshluger lepski
#Beispiel
set.seed(10)
n <- 100
X <- rnorm(n)
h <- goldenshluger_lepski(X)

## kernel estimator
#Beispiel
set.seed(10)
n <- 1000
X <- rnorm(n)
norm_estim <- kernel_estimator(X, dnorm, 1)
norm_estim(100)
dnorm(100)
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

## pco method
#Beispiel
set.seed(10)
n <- 1000
X <- rnorm(n)
h <- cross_validate(X)

hnorm_estim <- kernel_estimator(X, dnorm, 0.3)
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

## rejection sampling
#Beispiel
snorm <- rejection_sample(dnorm, rnorm, dnorm)
X <- snorm(100)
mean(X)

sunif <- rejection_sample(dnorm, runif, dunif)
X <- snorm(100)
mean(X)
