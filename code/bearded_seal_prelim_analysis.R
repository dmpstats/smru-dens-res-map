# ---------------------------------------------------------------------------------------------
#' Preliminary analysis for the Bearded Seal. 
#'
#'



# Preamble ------------------------------------------------------------------------------------

  library(tidyverse)
  
  rawRes <- readxl::read_xlsx("data/Bearded_seal_Density_RES_Data.xlsx", sheet = "Bearded s(P1)")
  
  # Coding of missing values are "-", several numeric fields interpreted as char
  
  beardSealData <- rawRes %>%
    mutate(across(c(`Uncertainty value`, `95% CI uncertainty low (per km2)`, `95% CI uncertainty high (per km2)`), as.numeric))
  

# Preliminary explorations --------------------------------------------------------------------

  # plot basic data to regress
  regPlot <- ggplot(beardSealData, x = `Species relative suitability index`, y = `Density estimate (per km2)`) + 
    geom_point(aes(`Species relative suitability index`, `Density estimate (per km2)`)) +
    geom_smooth(aes(x = `Species relative suitability index`, y = `Density estimate (per km2)`)) +
    ggthemes::theme_fivethirtyeight()

  regPlot  

  
  
  

# Naive modelling -----------------------------------------------------------------------------

  beardGAM <-   
  