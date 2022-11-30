# ---------------------------------------------------------------------------------------------
#' Preliminary analysis for the Bearded Seal. 
#'
#'



# Preamble ------------------------------------------------------------------------------------

  library(tidyverse)
  library(mgcv)
  
  rawRes <- readxl::read_xlsx("data/Bearded_seal_Density_RES_Data.xlsx", sheet = "Bearded s(P1)")
  
  # Coding of missing values are "-", several numeric fields interpreted as char
  
  beardSealData <- rawRes %>%
    mutate(across(c(`Uncertainty value`, `95% CI uncertainty low (per km2)`, `95% CI uncertainty high (per km2)`), as.numeric)) %>%
    rename(RES = `Species relative suitability index`, Density = `Density estimate (per km2)`) # gam doesn't like these names
  

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
  