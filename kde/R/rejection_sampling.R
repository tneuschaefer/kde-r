#' Rejection Sampling Factory
#'
#' @description \code{rejection_sampling_factory} creates a function that can draw samples
#'   from a given probabilty density function
#'
#' @param sample_density the probability density function to construct the sample
#' @param proposal_dist the distribution function used to create samples from
#'   \code{density}
#' @param proposal_density the probability density function used to accept or discard samples
#' @param m numeric scalar bigger one satisfying \code{sample_density(x) <= m*density(x)}
#' @param interval numeric interval used to calculate m if no m is given
#'
#' @details Rejection sampling uses \code{distribution} to draw samples and
#'   accepts/rejects these samples according to the densities \code{sample_density}
#'   and \code{density}, such that the resulting samples are \code{sample_density}-distributed.
#'
#'   Many rejected samples result in longer runtimes. To prevent this \code{m} should be
#'   chosen as small as possible, satisfying \code{sample_density(x) <=
#'   m*density(x)} for all \code{x}.
#'
#'   If no m is given \code{rejection_sampling_factory}
#'   searches for a maximum of \code{sample_density/proposal_density} inside \code{interval}
#'
#' @return A function taking a single numeric number \code{n} returning \code{n}
#'   \code{sample_density}-distributed random numbers
#'
#' @source https://en.wikipedia.org/wiki/Rejection_sampling
#'
#' @examples
#' custom_den <- function(x) {
#'   res <- (3 / 2) * (x^3) + (11 / 8) * (x^2) + (1 / 6) * (x) + (1 / 12)
#'   res[x < 0 | x > 1] <- 0
#'   res
#' }
#'
#' custom_sample <- rejection_sampling_factory(custom_den, stats::runif, stats::dunif, 3.2)
#' x <- seq(-0.5, 1.5, by = 0.01)
#' y <- custom_den(x)
#' sample <- custom_sample(100)
#'
#' plot(x, y,
#'   type = "l", main = "cutom_density: (3/2)*(x^3)+(11/8)*(x^2)+(1/6)*(x)+(1/12)",
#'   y_lab = "density"
#' )
#' points(sample, rep(0, 100), col = "red", pch = "o")
#' legend("topleft",
#'   legend = c("custom_density", "sample"), col = c("black", "red"),
#'   pch = c("|", "o")
#' )
#'
#' @export
rejection_sampling_factory <- function(sample_density, proposal_dist, proposal_density, m, interval) {
  # densities and distribution must be functions
  stopifnot(
    "sample_density must be a function" = is.function(sample_density),
    "proposal_dist must be a function" = is.function(proposal_dist),
    "proposal_density must be a function" = is.function(proposal_density)
  )

  if (!missing(m)) {
    # check requirements for m
    stopifnot(
      "m must be numeric" = is.numeric(m),
      "m must not be empty" = length(m) > 0,
      "m must be >= 1" = m >= 1
    )
    # in case the provided m is of length greater than 1
    m <- m[1]
  } else if (!missing(interval)) {
    # use interval to calculate m

    # check requirements for interval
    stopifnot(
      "interval must be numeric" = is.numeric(interval),
      "interval must be of length 2" = length(interval) == 2
    )

    quot <- function(x) sample_density(x) / proposal_density(x)
    m <- unlist(stats::optimize(quot, interval = interval, maximum = TRUE)[2])
  } else {
    stop("either m or interval must be provided")
  }

  # returning this function
  function(n) {
    # checking requirements for n
    stopifnot(
      "n must be numeric" = is.numeric(n),
      "n must be greater than 0" = n > 0
    )
    n <- as.integer(n[1])

    accepted <- c()

    # repeat until n values of the sample distribution have been found
    while (length(accepted) < n) {
      # uniformly pick u
      u <- stats::runif(m * n)

      # sample using the provided proposal distribution
      sample <- proposal_dist(m * n)

      # reject all values that are to large
      accepted <- c(accepted, sample[u * m * proposal_density(sample) < sample_density(sample)])
    }

    # return vector of length n
    accepted[1:n]
  }
}

#' Rejection Sampling
#'
#' @description \code{rejection_sampling} draws samples from a given probabilty density function
#'
#' @param n size of the sample to be drawn
#' @param sample_density the probability density function to construct the sample
#' @param proposal_dist the distribution function used to create samples from
#'   \code{density}
#' @param proposal_density the probability density function used to accept or discard samples
#' @param m numeric scalar bigger one satisfying \code{sample_density(x) <= m*density(x)}
#' @param interval numeric interval used to calculate m if no m is given
#'
#' @details Rejection sampling uses \code{distribution} to draw samples and
#'   accepts/rejects these samples according to the densities \code{sample_density}
#'   and \code{density}, such that the resulting samples are \code{sample_density}-distributed.
#'
#'   Many rejected samples result in longer runtimes. To prevent this \code{m} should be
#'   chosen as small as possible, satisfying \code{sample_density(x) <=
#'   m*density(x)} for all \code{x}.
#'
#'   If no m is given \code{rejection_sampling}
#'   searches for a maximum of \code{sample_density/density} inside \code{interval}
#'
#' @return \code{n} \code{sample_density}-distributed random numbers
#'
#' @examples
#' custom_den <- function(x) {
#'   res <- (3 / 2) * (x^3) + (11 / 8) * (x^2) + (1 / 6) * (x) + (1 / 12)
#'   res[x < 0 | x > 1] <- 0
#'   res
#' }
#'
#' sample <- rejection_sampling(100, custom_den, runif, dunif, 3.2)
#' x <- seq(-0.5, 1.5, by = 0.01)
#' y <- custom_den(x)
#'
#' plot(x, y,
#'   type = "l", main = "cutom_density: (3/2)*(x^3)+(11/8)*(x^2)+(1/6)*(x)+(1/12)",
#'   y_lab = "density"
#' )
#' points(sample, rep(0, 100), col = "red", pch = "o")
#' legend("topleft",
#'   legend = c("custom_density", "sample"), col = c("black", "red"),
#'   pch = c("|", "o")
#' )
#'
#' @export
rejection_sampling <- function(n, sample_density, proposal_dist, proposal_density, m, interval) {
  sampling_function <- rejection_sampling_factory(sample_density, proposal_dist, proposal_density, m, interval)
  sampling_function(n)
}
