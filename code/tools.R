# Conversion functions for the lognormal ------------------------------------------------------

logMu <- function(mu, sigma){
  
  
  log(mu^2/sqrt(mu^2 + sigma^2))
  
  }

logSigma <- function(mu, sigma){
  
  sqrt(log(1+sigma^2/mu^2))
  
}



# SE from CI ----------------------------------------------------------------------------------

  SEfromCI <- function(mu, lower, upper){
    
    logSE <- (log(upper) - log(lower))/(2*1.96)
    
    sigma <- sqrt((exp(logSE^2)-1)*mu^2)
    
    sigma
      
  }


# CI sampler ----------------------------------------------------------------------------------

  sampleLN <- function(mu, sigma, nSamp){
    
    if(mu == 0){
      outSamp <- rep(0, nSamp)
      } else {
    outSamp <- rlnorm(nSamp, logMu(mu, sigma), logSigma(mu, sigma))
      }
    
    outSamp
    
  }
