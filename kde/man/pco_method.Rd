% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pco_method.R
\name{pco_method}
\alias{pco_method}
\title{Penalized Comparison to Overfitting}
\source{
Lacour et al, Estimator selection: a new method with applications to
kernel density estimation (2017), https://arxiv.org/abs/1607.05091
}
\usage{
pco_method(
  x,
  kernel = stats::dnorm,
  n = 40,
  lambda = 1,
  N = 100L,
  built_in = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight",
    "silverman"),
  na.rm = FALSE
)
}
\arguments{
\item{x}{numeric vector of the observation sample}

\item{kernel}{kernel function used for kernel density estimation}

\item{n}{number of bandwidths to be optimized from. \code{pco_method} selects
a bandwidth contained in  (1/n, 2/n, ..., 1)}

\item{lambda}{positive scalar used as tuning parameter}

\item{N}{number of subdivisions used in discretization of integrals}

\item{built_in}{choose one of the built-in kernels instead of providing one yourself}

\item{na.rm}{logical; if TRUE, missing values will be removed from x}
}
\value{
The estimated optimal bandwidth
}
\description{
\code{pco_method} estimates an optimal bandwidth
for kernel density estimator using the PCO method
}
\details{
The PCO method tries to minimize an upper bound for the mean
integrated squared error (MISE) of a kernel density estimator

\code{pco_method} calculates the PCO criterion value for every bandwidth,
approximating the risk and selecting the bandwidth with the minimal risk

The bias/variance decomposition is used. Because the bias term depends on
the unknown density, a comparison of the estimator with an associated
bandwidth is used to estimate the bias term. Here the estimator with the
smallest bandwidth is used, this estimator is surely overfitting

Additional a penalty term is computed as the sum of the risk discomposition
variance and the variance of the bias term estimation. In the calculation
the tuning paramter \code{lambda} is used. Recommended is \code{lambda}=1

The PCO criterion is given as the sum of the comparison to the overfitting
and the penalty term. The PCO method tries to find a balance between this
two terms, therefor the name \code{penalized comparison to overfitting}
}
\examples{

x <- stats::rnorm(100)
h <- pco_method(x, kernel = stats::dnorm)
f <- kernel_estimator(x, kernel = stats::dnorm, bandwidth = h)

a <- min(x)
b <- max(x)
ab <- seq(a, b, length.out = 100)

plot(ab, stats::dnorm(ab), type = "l")
lines(ab, f(ab), col = "red")
legend("topleft",
  legend = c("true", "estimated h"), col = c("black", "red"),
  pch = "|"
)

}
\seealso{
\code{\link{kernel_estimator}} for more information about kernel
density estimation

\code{\link{cross_validation}} and \code{\link{goldenshluger_lepski}} for
other automatic bandwidth-selection algorithms
}
