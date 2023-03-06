# Kerbndichteschätzer Aufgabe 2: 
# Fixed bandwidth estimator

# h = Bandbreite

fix_band_est <- function(h){
  # Simulation of the model
  n <- 1000
  X <- rnorm(n)
  
  # Computation of a kernel estimator with a gaussian and an Epanechnikov kernel
  a <- min(X)
  b <- max(X)
  
  Xplot <- seq(a, b, length.out = 100)
  fvraie <- rep(0, 100)
  fvraie <- exp(-0.5*Xplot^2)/sqrt(2*pi)
  
  fest <- rep(0, 100)
  fest2 <- rep(0, 100)
  Mat <- matrix(0, nrow = 100, ncol = n)
  EPA <- matrix(0, nrow = 100, ncol = n)
  
  for (i in 1:100) {
    for (j in 1:n) {
      U <- (Xplot[i]-X[j])/h
      Mat[i,j] <- exp(-0.5*U^2)/sqrt(2*pi)/h
      T <- (abs((Xplot[i]-X[j])/h) < 1)
      EPA[i,j] <- 0.75*(1-U^2)*T/h
    }
  }
  
  fest <- apply(Mat, 1, mean)
  fest2 <- apply(EPA, 1, mean)
  
  # WIP: change output?
  plot(Xplot, fvraie, type = "l", col = "red", lwd = 2, xlab = "X", ylab = "Density")
  lines(Xplot, fest, col = "blue")
  lines(Xplot, fest2, col = "black", lty = "dashed")
}

