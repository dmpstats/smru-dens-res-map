# ---------------------------------------------------------------------------------------------
#' Analysis file for the SAR data
#' 
#' 


# Preamble ------------------------------------------------------------------------------------

  library(tidyverse)

  source("code/tools.R")  
  
  abundTable <- readxl::read_xlsx("data/SAR data.xlsx", sheet = "Abundance estimates", na = "-")
  
  NESRES <- readxl::read_xlsx("data/SAR data.xlsx", sheet = "Elephant seal", na = "-")
  harpRES <- readxl::read_xlsx("data/SAR data.xlsx", sheet = "Harp seal", na = "-")
  hoodRES <- readxl::read_xlsx("data/SAR data.xlsx", sheet = "Hooded seal", na = "-")
  


# Calc SEs ----------------------------------------------------------------

  abundTable <- abundTable %>%
    mutate(calcSE = SEfromCI(`Abundance estimate`, `Lower CI`, `Upper CI`),
           SE = ifelse(is.na(SE), calcSE, SE),
           CV = ifelse(is.na(CV), SE/`Abundance estimate`, CV))

# NES  --------------------------------------------------------------------

  NESSum_Res <- sum(NESRES$`Overall Pr`*NESRES$Area)
  NESAbund <- abundTable %>% filter(Species == "Northern elephant seal") 
  
  NES_predictions <- NESRES %>%
    mutate(weight = (`Overall Pr`*Area)/NESSum_Res,
           predAbund = weight * NESAbund$`Abundance estimate`,
           predDensity = predAbund/Area,
           CV = NESAbund$CV[1]) %>%
    rename(RES = `Overall Pr`)
  
  # check calcs
  sum(NES_predictions$predAbund); NESAbund$`Abundance estimate`
  
  
  NES_predictions <- NES_predictions %>%
    select(Species:RES, predDensity, CV)
    
  write_csv(NES_predictions, "data/predictions/Northern_Elephant_Seal_predictions.csv")
    


# Hooded seal -------------------------------------------------------------

  hoodAbund <- abundTable %>%
    filter(Species == "Hooded seal")
  
  dataList <- split(hoodAbund, hoodAbund$Stock)
  
  set.seed(4835)
  
  abundResamples <- lapply(dataList, function(q){sampleLN(q$`Abundance estimate`, q$SE, 1000)}) %>%
    plyr::ldply() %>%
    rename(Stock = .id) %>%
    pivot_longer(names_to = "rep",  values_to = "abund", -Stock) %>%
    pivot_wider(names_from = Stock, values_from = abund) %>%
    mutate(combinedAbund = `Greenland Sea` + `NW Atlantic`)
    
  hoodCombinedAbund <- tibble(`Abundance estimate` = mean(abundResamples$combinedAbund), 
                                  SE = sd(abundResamples$combinedAbund)) %>%
    mutate(CV = SE/`Abundance estimate`)
    
  hoodSum_Res <- sum(hoodRES$`Overall Pr`*hoodRES$Area)
  
    hood_predictions <- hoodRES %>%
      mutate(weight = (`Overall Pr`*Area)/hoodSum_Res,
             predAbund = weight * hoodCombinedAbund$`Abundance estimate`,
             predDensity = predAbund/Area,
             CV = hoodCombinedAbund$CV) %>%
      rename(RES = `Overall Pr`)
    
    # check calcs
    sum(hood_predictions$predAbund); hoodCombinedAbund$`Abundance estimate`
    
    
    hood_predictions <- hood_predictions %>%
      select(Species:RES, predDensity, CV)
    
    write_csv(hood_predictions, "data/predictions/Hooded_Seal_predictions.csv")
    

# Harp seals --------------------------------------------------------------

 
     
    harpAbund <- abundTable %>%
      filter(Species == "Harp seal")
    
    dataList <- split(harpAbund, harpAbund$Stock)
    
    set.seed(4835)
    
    abundResamples <- lapply(dataList, function(q){sampleLN(q$`Abundance estimate`, q$SE, 1000)}) %>%
      plyr::ldply() %>%
      rename(Stock = .id) %>%
      pivot_longer(names_to = "rep",  values_to = "abund", -Stock) %>%
      pivot_wider(names_from = Stock, values_from = abund) %>%
      mutate(combinedAbund = `Greenland Sea` + `NW Atlantic` + `Barents Sea/White Sea`)
    
    harpCombinedAbund <- tibble(`Abundance estimate` = mean(abundResamples$combinedAbund), 
                                SE = sd(abundResamples$combinedAbund)) %>%
      mutate(CV = SE/`Abundance estimate`)
    
    harpSum_Res <- sum(harpRES$`Overall Pr`*harpRES$Area)
    
    harp_predictions <- harpRES %>%
      mutate(weight = (`Overall Pr`*Area)/harpSum_Res,
             predAbund = weight * harpCombinedAbund$`Abundance estimate`,
             predDensity = predAbund/Area,
             CV = harpCombinedAbund$CV) %>%
      rename(RES = `Overall Pr`)
    
    # check calcs
    sum(harp_predictions$predAbund); harpCombinedAbund$`Abundance estimate`
    
    
    harp_predictions <- harp_predictions %>%
      select(Species:RES, predDensity, CV)
    
    write_csv(harp_predictions, "data/predictions/Harp_Seal_predictions.csv")
    