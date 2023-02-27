# Rejection-Sampling-Funktion,
# die gemäß gegebener Dichte zieht

# Main sources:
# https://www.youtube.com/watch?v=kMb4JlvuGlw
# https://rpubs.com/mathetal/rejectsampling

# fnct = Funktion f(x)
# n_samples = Anzahl gewünschter Stichproben
# min = ???
# max = ???

rejection_sampling <- function(fnct, n_samples, min, max){

  # Anzahl benötigter samples
  n <- (n_samples * 5)
  # Konstante c berechnen
  c <- (fnct(max))

  X = runif(n, min, max)
  U = runif(n, min, max)

  accept = c()
  i = 1

  while(i <= n & length(accept) < n_samples){
    test_u = U[i]
    test_x = fnct(X[i])/(c*dunif(X[i],min,max))
    if (test_u <= test_x){
      accept = rbind(accept, X[i])
      i = i + 1
    }
    i = i + 1
  }

  hist(accept)
}

## -----------------------------------------------------------------
# Beispiel:

pi_x <- function(x) {
  new_x = (3/2)*(x^3)+(11/8)*(x^2)+(1/6)*(x)+(1/12)
  return(new_x)
}

rejection_sampling(pi_x, 1000, 0, 1)
