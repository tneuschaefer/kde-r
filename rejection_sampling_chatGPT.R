# Rejection-Sampling-Function,
# die gemäß gegebener Dichte zieht

##! noch korrigieren


# target_dist_mean: Mittelwert der Normalverteilung der Zielverteilung
# target_dist_sd: Standardabweichung der Normalverteilung der Zielverteilung
# prop_dist_min = Minimum der Schätzverteilung
# prop_dist_max = Maximum der Schätzverteilung
# M = Maximalwert von acceptance_prob
# n_samples = Anzahl der gewünschten Stichproben


rejection_sampling <- function(target_dist_mean, target_dist_sd, prop_dist_min, prop_dist_max, M, n_samples) {

  ##! ERROR, falls Werte fehlen

  # Zielverteilung: Normalverteilung mit Mittelwert und Standardabweichung
  target_dist <- function(x) dnorm(x, mean = target_dist_mean, sd = target_dist_sd)

  # Schätzverteilung: Gleichverteilung auf Intervall [min,max]
  proposal_dist <- function(x) dunif(x, min = prop_dist_min, max = prop_dist_max)

  # Funktion zur Akzeptanz- und Ablehnungswahrscheinlichkeit
  acceptance_prob <- function(x) target_dist(x) / (M * proposal_dist(x))

  # Generierung von Stichproben
  samples <- numeric(n_samples)
  i <- 1
  while (i <= n_samples) {         ##! Video baut hier gewünschte n_samples ein.
    # Ziehen einer Stichprobe aus der Schätzverteilung
    x <- runif(1, min = prop_dist_min, max = prop_dist_max)

    # Akzeptieren oder Ablehnen der Stichprobe mit Hilfe von acceptance_prob
    u <- runif(1)
    if (u < acceptance_prob(x)) {
      samples[i] <- x
      i <- i + 1
      # Bei Ablehnung bleibt i gleich und es wird bis zu Akzeptanz neu gezogen.
      # So wird samples nur mit akzeptierten Werten gefüllt, bis gewünscht Anzahl an samples erreicht ist.
    }
  }
  # Plotten der Zielverteilung und der generierten Stichproben
  hist(samples, breaks = 20, freq = FALSE, main = "Histogram of Generated Samples")
  curve(target_dist, add = TRUE, col = "red", lwd = 2)
}

## -----------------------------------------------------------------
# Beispiel:

rejection_sampling(5, 2, 0, 10, 5, 1000)

