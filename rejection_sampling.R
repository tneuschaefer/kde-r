

rejection_sample <- function(n, sample_density, distribution, density, M = NULL) {
  stopifnot() #WIP
  
  quot <- function(x) {
    if(density(x) < 0) return(0)
    else return(sample_density(x) / density(x))
  }
  
  if(is.null(M)){
    M <- unlist(optimize(quot, c(-100,100), maximum = TRUE)[1]) #WIP
  }
  
  values <- c()
  i <- 1
  while(i <= n) {
    u <- runif(1)
    y <- distribution(1) #distribution must allow this form of generating values
    
    if(u*M*density(y) < sample_density(y)) { # density(y) < 0
      values <- c(values, y)
      i <- i+1
    }
  }
  return(values)
}

#Example

norm <- function(n) {
  return(rnorm(n, 10, 1))
}

pi_x <- function(x) {
  (3/2)*(x^3)+(11/8)*(x^2)+(1/6)*(x)+(1/12)
}

unif <- function(n) {
  return(runif(n))
}
