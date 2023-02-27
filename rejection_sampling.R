##! noch zu korrigieren

# Rejection-Sampling-Funktion:

# Main sources:
# https://en.wikipedia.org/wiki/Rejection_sampling
# https://de.wikipedia.org/wiki/Verwerfungsmethode
# https://www.youtube.com/watch?v=kMb4JlvuGlw
# https://rpubs.com/mathetal/rejectsampling

# target_dist = Zielfunktion f(x)
# n_samples = Anzahl gewünschter Stichproben
# prop_dist_min = Minimum der Schätzverteilung
# prop_dist_max = Maximum der Schätzverteilung


rejection_sampling <- function(target_dist, n_samples, prop_dist_min, prop_dist_max){
  
  ##! ERROR, falls Werte fehlen

  # Konstante c berechnen ##!(Maximalwert von Akzeptanzwahrscheinlichkeit)
  c <- (target_dist(prop_dist_max))


  # Generierung von Stichproben
  samples <- numeric(n_samples)
  i = 1

  while (i <= n_samples) {
    # Ziehen einer Stichprobe aus der Schätzverteilung
    x <- runif(1, min = prop_dist_min, max = prop_dist_max)
    u <- runif(1, min = prop_dist_min, max = prop_dist_max)

    if (u <= (target_dist(x)/(c*dunif(x,prop_dist_min,prop_dist_max)))) {
      samples[i] <- x
      i <- i + 1
      # Bei Ablehnung bleibt i gleich und es wird bis zur Akzeptanz neu gezogen.
      # So wird samples nur mit akzeptierten Werten befüllt, bis gewünscht Anzahl an samples erreicht ist.
    }
  }

  hist(samples) ##! möglicherweise andere Ausgabe?
}

## -----------------------------------------------------------------
# Beispiele:

pi_x <- function(x) {
  new_x = (3/2)*(x^3)+(11/8)*(x^2)+(1/6)*(x)+(1/12)
  return(new_x)
}

pi_x <- function(x){

  bla = (3*x^2 + x + 13)
  return(bla)
}

rejection_sampling(pi_x, 1000, 0, 1)

