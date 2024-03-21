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



# GAM fit function for bootstrapping --------------------------------------
#' Fits a monotone GAM, no intercept (although zero intercept not guaranteed)

gamFit <- function(inData, inRES){
  
  workingFit <- scam::scam(Density ~ s(RES, bs = "mpi")-1, data = inData)
  
  resGridPred <- scam::predict.scam(workingFit, newdata = inRES)
  
  outData <- inRES %>%
    mutate(Pred = resGridPred)
  
  outData
  
}


# GAM fit function for bootstrapping --------------------------------------
#' Fits a monotone GAM assuming log link tweedie distribution - p needs to be specified
#' Fits can be sensitive to this - p tending towards 1 is a Poisson-like dist, towards 2 is gamma
#' No attempt to constrain the intercept in this construction

gamFitTweedie <- function(inData, inP, inRES){
  
  workingFit <- scam::scam(Density ~ s(RES, bs = "mpi"), data = inData, family = Tweedie(p = inP, link = "log"))
  
  resGridPred <- scam::predict.scam(workingFit, newdata = inRES, type = "response")
  
  outData <- inRES %>%
    mutate(Pred = resGridPred)
  
  outData
  
}
