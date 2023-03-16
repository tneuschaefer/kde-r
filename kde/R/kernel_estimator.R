#' Title
#'
#' @param X
#' @param K
#' @param h
#'
#' @return
#' @export
#'
#' @examples
kernel_estimator <- function(X, K = stats::dnorm, h = 1) {

  #Sample condition
  stopifnot("X must be a numeric" = is.numeric(X),
            "X must be non empty" = length(X) > 0)

  #Kernel condition
  stopifnot("K is not a funtion" = is.function(K))

  #Bandwidth condition
  stopifnot("h must be numeric" = is.numeric(h),
            "h must be bigger zero" = h > 0,
            "h must have length one" = length(h) == 1)

  estimator <- function(x) {
    stopifnot("x must be numeric" = is.numeric(x),
              "x miust be non empty" = length(x) > 0)

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
