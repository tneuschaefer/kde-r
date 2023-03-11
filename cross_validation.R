# Kerbndichteschätzer Aufgabe 3: 
# Cross-validation bandwidth selection

library(stats)

cross_validation <- function(){
  
  # Simulation of the model
  nobs <- 1000
  alpha <- 3
  beta <- 3
  Z <- (rexp(nobs, rate = 1) - rexp(nobs, rate = 1))/1 # hierfür würde library(stats) benötigt
  
  # Domain of estimation and true function
  a <- min(Z)
  b <- max(Z)
  K <- 50
  ab <- seq(a, b, length.out = K)
  fvraie <- exp(-abs(ab))/2
  
  # Computation of criterion with Epanechnikov kernel
  M <- 40
  estimf <- matrix(0, nrow = M, ncol = K)
  Mat <- array(0, dim = c(nobs, nobs, M))
  Crit <- rep(0, M)
  
  for (k in 1:M) {
    h <- k/M
    T <- abs((matrix(rep(Z, K), nobs, byrow = TRUE) - matrix(rep(ab, nobs), nobs, byrow = FALSE))/h) < 1
    estimf[k, ] <- apply(0.75*(1-((matrix(rep(Z, K), nobs, byrow = TRUE) - matrix(rep(ab, nobs), nobs, byrow = FALSE))/h)^2)*T/h, 2, mean)
  }
  
  Un <- rep(1, nobs)
  R <- diag(Un, nobs, nobs)
  for (k in 1:M) {
    h <- k/M
    for (i in 1:nobs) {
      for (j in 1:nobs) {
        aa <- (abs((Z[i]-Z[j])/h) < 1)
        Mat[i, j, k] <- 0.75*(1-((Z[i]-Z[j])/h)^2)*aa/h-R[i, j]/h
      }
    }
  }
  
  for (k in 1:M) {
    Crit[k] <- sum(estimf[k, ]^2)*(b-a)/K-2*sum(Mat[, , k], na.rm = TRUE)/nobs/(nobs-1)
  }
  
  # Bandwidth selection
  kc <- which.min(Crit)
  estimfinal <- estimf[kc, ]
  
  # Graphical representation
  plot(ab, fvraie, type = 'l', col = 'red', lwd = 1.5, ylim = c(0, 1))
  lines(ab, estimfinal, type = 'l', col = 'blue', lwd = 1.5)
}

cross_validation()
