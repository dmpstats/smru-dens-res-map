# Conversion functions for the lognormal ------------------------------------------------------

logMu <- function(mu, sigma){
  
  
  log(mu^2/sqrt(mu^2 + sigma^2))
  
  }

logSigma <- function(mu, sigma){
  
  log(1+sigma^2/mu^2)
  
}
