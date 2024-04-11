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



# overarching function for batch submission -----------------------------------------------------------------------



fitFun <- function(speciesName, locationName, nBoot = 200, inData = workingData, inRES = workingRES, modelForm = "tweedie", inSeed = 345){
  
  dataList <- split(inData, inData$blockID)
  
  set.seed(inSeed)
  
  sampleList <- lapply(dataList, function(q){sampleLN(q$Density[1], q$workingSE[1], nBoot)})
  
  sampleDF <- plyr::ldply(sampleList) %>%
    rename(Survey = .id)
  
  speciesSamples <- inData %>% left_join(sampleDF)
  
  speciesSamples <- speciesSamples %>% 
    select(-Density) %>%
    pivot_longer(names_to = "sampleID", values_to = "Density", V1:last_col()) 
  
  speciesList <- split(speciesSamples, speciesSamples$sampleID)
  
  plan(multicore, workers = 2)
  
  modelFit <- function(inList){
    
    p <- progressor(steps = length(inList))
    
    if(modelForm == "tweedie"){  
      
      future_map(inList, \(x) {p(); gamFitTweedie(x, inP = 1.2, inRES = data.frame(RES = seq(0, 1, by = 0.01)))})
      
    } else {
      
      future_map(inList, \(x) {p(); gamFit(x, inRES = data.frame(RES = seq(0, 1, by = 0.01)))})
      
    }
    
  }
  
  
  with_progress({
    fittedList <- modelFit(speciesList)
  })
  
  
  fittedDF <- fittedList %>%
    bind_rows()
  
  fittedMatrix <- matrix(fittedDF$Pred, ncol = nBoot)
  
  bootEsts <- t(apply(fittedMatrix, 1, function(q){quantile(q, probs = c(0.025, 0.5, 0.975))}))
  
  bootSE <- apply(fittedMatrix, 1, sd)
  
  resFits <- as.data.frame(bootEsts) %>%
    mutate(RES = seq(0, 1, by = 0.01), SE = bootSE) %>%
    rename(lower = `2.5%`, med = `50%`, upper = `97.5%`) %>%
    mutate(med = ifelse(med < 0, min(abs(med)), med),
           CV = SE/med) %>%
    #CV = ifelse(CV > 2, 2, CV)) %>%
    arrange(RES)
  
  plottingDF <- resFits
  
  bootPlot <- ggplot(plottingDF) +
    ggthemes::theme_fivethirtyeight() +
    geom_line(aes(RES, med), size = 2, alpha = 0.7, col = "purple") +
    geom_ribbon(aes(x = RES, ymin = lower, ymax = upper), fill = "purple", alpha = 0.2) +
    ggtitle("Density as function of RES", paste0(speciesName, " whale: bootstrapped monotone spline fits"))
  
  ggsave(plot = bootPlot, here(paste0("docs/images/", speciesName, "_", locationName, "_bootplot_", nBoot, ".png")), units = "cm", width = 30, height = 20)
  saveRDS(plottingDF, here(paste0("data/plotting components/", speciesName, "_", locationName, "_plotElements",".rds")))
  
  
  # Create RES grid predictions -----------------------------------------------------------------
  
  resFits <- resFits %>% select(med, RES, CV) %>%
    rename(PredDensity = med) %>%
    mutate(RES = round(RES, 2))
  
  predictionOutput <- inRES %>% left_join(resFits, by = "RES")
  
  predictionOutput <- predictionOutput %>%
    mutate(PredDensity = ifelse(RES == 0, 0, PredDensity),
           CV = ifelse(RES == 0, NA, CV))
  
  summary(predictionOutput)
  
  write.csv(predictionOutput, file = here(paste0("data/predictions/", speciesName, "_", locationName, "_predictions.csv")), row.names = F)
  
}



# Right whale hack ------------------------------------------------------------------------------------------------



# overarching function for batch submission -----------------------------------------------------------------------



fitFun_RW <- function(speciesName, locationName, nBoot = 200, inData = workingData, inRES = workingRES, modelForm = "tweedie", inSeed = 345){
  
  dataList <- split(inData, inData$blockID)
  
  set.seed(inSeed)
  
  sampleList <- lapply(dataList, function(q){sampleLN(q$Density[1], q$workingSE[1], nBoot)})
  
  sampleDF <- plyr::ldply(sampleList) %>%
    rename(Survey = .id)
  
  speciesSamples <- inData %>% left_join(sampleDF)
  
  speciesSamples <- speciesSamples %>% 
    select(-Density) %>%
    pivot_longer(names_to = "sampleID", values_to = "Density", V1:last_col()) 
  
  speciesList <- split(speciesSamples, speciesSamples$sampleID)
  
  plan(multicore, workers = 2)
  
  modelFit <- function(inList){
    
    p <- progressor(steps = length(inList))
    
    if(modelForm == "tweedie"){  
      
      future_map(inList, \(x) {p(); gamFitTweedie(x, inP = 1.2, inRES = data.frame(RES = seq(0, 1, by = 0.01)))})
      
    } else {
      
      future_map(inList, \(x) {p(); gamFit(x, inRES = data.frame(RES = seq(0, 1, by = 0.01)))})
      
    }
    
  }
  
  
  with_progress({
    fittedList <- modelFit(speciesList)
  })
  
  
  fittedDF <- fittedList %>%
    bind_rows()
  
  fittedMatrix <- matrix(fittedDF$Pred, ncol = nBoot)
  
  bootEsts <- t(apply(fittedMatrix, 1, function(q){quantile(q, probs = c(0.025, 0.5, 0.975))}))
  
  bootSE <- apply(fittedMatrix, 1, sd)
  
  resFits <- as.data.frame(bootEsts) %>%
    mutate(RES = seq(0, 1, by = 0.01), SE = bootSE) %>%
    rename(lower = `2.5%`, med = `50%`, upper = `97.5%`) %>%
    mutate(lower = if_else(lower < 0, 0, lower),
           upper = if_else(upper < 0, 0, upper),
           med = ifelse(med < 0, 0, med),
           CV = SE/med,
           CV = if_else(med == 0, NA, CV)         
    ) %>% 
    arrange(RES)
  
  plottingDF <- resFits
  
  bootPlot <- ggplot(plottingDF) +
    ggthemes::theme_fivethirtyeight() +
    geom_line(aes(RES, med), size = 2, alpha = 0.7, col = "purple") +
    geom_ribbon(aes(x = RES, ymin = lower, ymax = upper), fill = "purple", alpha = 0.2) +
    ggtitle("Density as function of RES", paste0(speciesName, " whale: bootstrapped monotone spline fits"))
  
  ggsave(plot = bootPlot, here(paste0("docs/images/", speciesName, "_", locationName, "_bootplot_", nBoot, ".png")), units = "cm", width = 30, height = 20)
  saveRDS(plottingDF, here(paste0("data/plotting components/", speciesName, "_", locationName, "_plotElements",".rds")))
  
  
  # Create RES grid predictions -----------------------------------------------------------------
  
  resFits <- resFits %>% select(med, RES, CV) %>%
    rename(PredDensity = med) %>%
    mutate(RES = round(RES, 2))
  
  predictionOutput <- inRES %>% left_join(resFits, by = "RES")
  
  predictionOutput <- predictionOutput %>%
    mutate(PredDensity = ifelse(RES == 0, 0, PredDensity),
           CV = ifelse(RES == 0, NA, CV))
  
  summary(predictionOutput)
  
  write.csv(predictionOutput, file = here(paste0("data/predictions/", speciesName, "_", locationName, "_predictions.csv")), row.names = F)
  
}



