

rejection_sample <- function(sample_density, distribution, density,
                             M = NULL, lower = NULL, upper = NULL) {
  
  #Sample_density condition
  stopifnot("Sample_density must be a function" = is.function(sample_density))
  
  #Distribution condition
  stopifnot("Density must be a function" = is.function(density),
            "Distribution must be a function" = is.function(distribution))
  
  #Upper Bound condition
  if(!is.null(M)) stopifnot("M must be numeric" = is.numeric(M),
                            "M must have length one" = length(M) == 1,
                            "Must be > one" = M > 1)
  #Interval condition
  
  
  quot <- function(x) {
    if(density(x) < 0) return(1)
    else return(sample_density(x) / density(x))
  }
  
  if(is.null(M)){
    M <- unlist(optimize(quot, c(min(X),max(X)), maximum = TRUE)[1]) #WIP
  }
  
  if(M < 1) M <- 1.1
  
  function(n) {
    
    #Sampling condition
    stopifnot("n must be numeric" = is.numeric(n),
              "n must havee length one" = length(n) == 1)
    
    n <- as.integer(n)
    u <- numeric(M*n)
    sample <- numeric(M*n)
    accepted <- c()
    
    while(length(accepted) < n) {
      u <- runif(M*n)
      sample <- distribution(M*n)
      accepted <- c(accepted, sample[u*M*density(sample) < sample_density(sample)])
    }
    accepted[1:n]
  }
}

#Beispiel
snorm <- rejection_sample(dnorm, rnorm, dnorm)
X <- snorm(100)
mean(X)

sunif <- rejection_sample(dnorm, runif, dunif)
X <- snorm(100)
mean(X)
