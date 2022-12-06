# ---------------------------------------------------------------------------------------------
#' Preliminary analysis for the Bearded Seal. 
#'
#'



# Preamble ------------------------------------------------------------------------------------

  library(tidyverse)
  library(mgcv)
  
  source("code/tools.R")  
  
  rawRes <- readxl::read_xlsx("data/Bearded_seal_Density_RES_Data.xlsx", sheet = "Bearded s(P1)")
  
  # Coding of missing values are "-", several numeric fields interpreted as char
  
  beardSealData <- rawRes %>%
    mutate(across(c(`Uncertainty value`, `95% CI uncertainty low (per km2)`, `95% CI uncertainty high (per km2)`), as.numeric),
           `Uncertainty measure` = ifelse(`Uncertainty measure` == "-", NA, `Uncertainty measure`)) %>%
    rename(RES = `Species relative suitability index`, Density = `Density estimate (per km2)`,
           LowerCI = `95% CI uncertainty low (per km2)`,
           UpperCI = `95% CI uncertainty high (per km2)`) # gam doesn't like these names

  

# Preliminary explorations --------------------------------------------------------------------

  # plot basic data to regress
  regPlot <- ggplot(beardSealData, x = RES, y = Density) + 
    geom_point(aes(RES, Density), alpha = 0.4, position = position_jitter()) +
    geom_smooth(aes(RES, Density)) +
    ggthemes::theme_fivethirtyeight()

  regPlot  

  
  
  

# Naive modelling -----------------------------------------------------------------------------

  
  beardGAM <- gam(Density ~ s(RES, bs = "cs"), data = beardSealData, family = gaussian(link = "log"))
  plot(beardGAM)
  plot(beardSealData$RES, predict(beardGAM, type = "response"))
  
  
  # monotonicity constraints
  monGAM <- scam::scam(Density ~ s(RES, bs = "mpi")-1, data = beardSealData)
  
  plot(monGAM)
  
  predObj <- predict(monGAM, se.fit = T)
  monGAMPred <- data.frame(pred = predObj$fit, SE = predObj$se.fit) %>%
    mutate(lower = pred - 2*SE, upper = pred + 2*SE)
  
  monPredData <- beardSealData %>% 
    bind_cols(monGAMPred)
  
  
  monPlot <- ggplot(monPredData) +
    geom_point(aes(RES, Density), alpha = 0.4) +
    geom_line(aes(RES, pred), size = 2, col = "purple", alpha = 0.6) +
    ggthemes::theme_fivethirtyeight()
    
  monPlot
  
  
  

# Adding survey uncertainty -------------------------------------------------------------------
#' Here devise resampling for the different sorts of uncertainty that are present in the survey data
#' Note, there are single measures and upper/lower 95% CIS
#' 

  table(beardSealData$`Uncertainty measure`, useNA = "always")
  
  # check that where there is not a single measure, we do have a CI
  test <- beardSealData %>% filter(is.na(`Uncertainty measure`)) 
  any(is.na(test$UpperCI))
    
  # note some have CIs and alterative measure. Will use CIs when available
  # %-age CV looks to just be CV * 100
  # convert to CVs, extract as CIs so single sampling method required
  # everything in the data with a CV provided a CI as well
  
  beardSealData <- beardSealData %>%
    mutate(`Uncertainty value` = ifelse(`Uncertainty measure` == "% CV", `Uncertainty value`/100, `Uncertainty value`),
           `Uncertainty measure` = ifelse(`Uncertainty measure` == "% CV", "CV", `Uncertainty measure`),
           `Uncertainty value` = ifelse(`Uncertainty value` == 7000, 0.0402, `Uncertainty value`), # hacks as need abundance to convert
           `Uncertainty value` = ifelse(`Uncertainty value` == 5000, 0.04, `Uncertainty value`), # hacks as need abundance to convert
           `Uncertainty measure` = ifelse(`Uncertainty measure` == "SE value for abundance", "SE", `Uncertainty measure`),
           `Uncertainty value` = ifelse(`Uncertainty measure` == "CV", `Uncertainty value`*Density, `Uncertainty value`),
           `Uncertainty measure` = ifelse(`Uncertainty measure` == "CV", "SE", `Uncertainty measure`),
           `Uncertainty measure` = ifelse(is.na(`Uncertainty value`), "SE", `Uncertainty measure`),
            workingSE = `Uncertainty value`,
            workingSE = ifelse(is.na(workingSE), SEfromCI(Density, LowerCI, UpperCI), workingSE)
           )
  
  # obtain CIs from lognormal using mean, SE
  
  beardSealData <- beardSealData %>%
      mutate(workingLowerCI = qlnorm(0.025, logMu(Density, workingSE), logSigma(Density, workingSE)),
             workingUpperCI = qlnorm(0.975, logMu(Density, workingSE), logSigma(Density, workingSE)),
             lowerError = LowerCI - workingLowerCI, 
             upperError = UpperCI - workingUpperCI,
             blockID = paste(beardSealData$`Survey ID`, beardSealData$`Area/Segment`, sep = ":")
      )

  # sample for each of the area/segments collectively (as not independent)  
  
  dataList <- split(beardSealData, beardSealData$blockID)
  
  sampleList <- lapply(dataList, function(q){sampleLN(q$Density[1], q$workingSE[1], 10)})
  
  sampleDF <- plyr::ldply(sampleList) %>%
    mutate(`Area/Segment` = str_extract(.id, "(?<=:).+"),
           `Survey ID` = str_extract(.id, ".+(?=:)")) %>%
    select(-.id)
    
  
  beardSealSamples <- beardSealData %>% left_join(sampleDF)
  
  
  beardSealSamples %>% 
    select(-Density) %>%
    pivot_longer(names_to = "sampleID", values_to = "Density", V1:V10) %>%
    group_by(sampleID) %>%
    nest() %>%
    mutate(PredVals =  map_df(data, function(q){modFit <- scam::scam(Density ~ s(RES, bs = "mpi")-1, data = q); data.frame(fitted = modFit$fitted.values)})) %>%
    mutate(PredVals = map2(data, PredVals, bind_cols))
  
    
  

# Scratch -------------------------------------------------------------------------------------

  
  test <- seq(0.2, 0.6, length = 100)
  
  testD <- dlnorm(test, logMu(0.5, 0.5), logSigma(0.5, 0.5))
  
  qlnorm(c(0.025, 0.975), logMu(0.372, 0.2*0.372), logSigma(0.372, 0.2*0.372))
  
  plnorm(c(0.252, 0.549), logMu(0.372,  0.07425), logSigma(0.372,  0.07425))
  
  nRep<- 1000
  testMean <- numeric(nRep)
  testSD <- numeric(nRep)
  testLogSD <- numeric(nRep)
  
  for(i in 1:nRep){
    n <- 100
    testSamp <- rlnorm(n, logMu(0.008, 0.011), logSigma(0.008, 0.011))
    testMean[i] <- mean(testSamp)
    testSD[i] <- sd(testSamp)
    testLogSD <- sd(log(testSamp))
    
  }
  
  mean(testMean)
  mean(testSD)

  
  beardSealData
  
  refpts <- c(0.252, 0.372, 0.549)

  rriskDistributions::fit.perc(c(0.025, 0.5, 0.975), refpts)  
  