#' Rejection sampling
#'
#' @description \code{rejection_sample} creates a function that draws samples
#'   from a given probabilty density function
#'
#' @param sample_density the probability density function to construct the sample
#' @param distribution the distribution function used to create samples from
#'   \code{density}
#' @param density the probability density function used to accept or discard samples
#' @param M numeric scalar bigger one satisfying \code{sample_density(x) <= M*density(x)} 
#' @param interval numeric interval used to calculate M if no M is given
#' 
#' @details Rejection sampling uses \code{distribution} to draw samples and 
#'   accepts/rejects these samples according to the densities \code{sample_density} 
#'   and \code{density}, such that the resulting samples are \code{sample_density}-distributed. 
#'   
#'   Many rejected samples result in longer runtimes. To prevent this \code{M} should be
#'   chosen as small as possible, satisfying \code{sample_density(x) <=
#'   M*density(x)} for all \code{x}. 
#'   
#'   If no M is given \code{rejection_sample}
#'   searches for a maximum of \code{sample_density/density} inside \code{interval}
#'   
#' @return A function taking a single numeric number \code{n} returning \code{n}
#'   \code{sample_density}-distributed random numbers
#'   
#' @source https://en.wikipedia.org/wiki/Rejection_sampling
#' 
#' @examples 
#'   custom_den <- function(x) {
#'     res <- (3/2)*(x^3)+(11/8)*(x^2)+(1/6)*(x)+(1/12)
#'     res[x<0|x>1] <- 0
#'     res
#'   }
#'   
#'   custom_sample <- rejection_sample(custom_den, runif, dunif, 3.2)
#'   x <- seq(-0.5, 1.5, by = 0.01)
#'   y <- custom_den(x)
#'   sample <- custom_sample(100)
#'   
#'   plot(x, y, type = "l", main = "cutom_density: (3/2)*(x^3)+(11/8)*(x^2)+(1/6)*(x)+(1/12)",
#'     y_lab = "density")
#'   points(sample, rep(0, 100), col = "red", pch = "o")
#'   legend("topleft", legend = c("custom_density", "sample"), col = c("black", "red"),
#'     pch = c("|", "o"))
#' 
#' @export
rejection_sample <- function(sample_density, distribution, density, M = NULL,
                             interval = NULL) {

  #Sample_density condition
  stopifnot("Sample_density must be a function" = is.function(sample_density))

  #Distribution condition
  stopifnot("Density must be a function" = is.function(density),
            "Distribution must be a function" = is.function(distribution))
  
  stopifnot("Either M or interval must be != NULL" = !is.null(M) || !is.null(interval))

  #Upper Bound condition
  if(!is.null(M)) stopifnot("M must be numeric" = is.numeric(M),
                            "M must have length one" = length(M) == 1,
                            "Must be > one" = M > 1)
  #Interval condition
  if(is.null(M)) stopifnot("interval must be numeric" = is.numeric(interval),
                           "interval must be of length 2" = length(interval) == 2)


  if(is.null(M)){
    quot <- function(x) sample_density(x)/density(x)
    M <- unlist(stats::optimize(quot, interval = interval, maximum = TRUE)[2])
  }

  if(M <= 1) M <- 1.1

  function(n) {

    #Sampling condition
    stopifnot("n must be numeric" = is.numeric(n),
              "n must havee length one" = length(n) == 1)

    n <- as.integer(n)
    u <- numeric(M*n)
    sample <- numeric(M*n)
    accepted <- c()

    while(length(accepted) < n) {
      u <- stats::runif(M*n)
      sample <- distribution(M*n)
      accepted <- c(accepted, sample[u*M*density(sample) < sample_density(sample)])
    }
    accepted[1:n]
  }
}
