

kernel_estimator <- function(X, K, h) {
  stopifnot("X is invalid sample" = is.numeric(X) && length(X) > 0,
            "K is not a function" = is.function(K),
            "h is no number" = is.numeric(h) && length(h) == 1)
  
  estimator <- function(x) {
    stopifnot("x is invalid" = is.numeric(x) && length(x) > 0)
    
    n <- length(X)
    m <- length(x)
    Mat <- matrix(0, nrow = n, ncol = m)
    
    for(i in 1:n) {
      for(j in 1:m) {
        U <- (X[i]-x[j])/h
        Mat[i,j] <- K(U)/h
      }
    }
    apply(Mat, 2, mean)
  }
}

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
