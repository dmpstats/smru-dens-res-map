# ---------------------------------------------------------------------------------------------
#' Analysis file for the narwhal
#' 
#' 


# Preamble ------------------------------------------------------------------------------------

  library(tidyverse)
  library(mgcv)
  
  source("code/tools.R")  
  
  # note missing values in spreadsheet indicated by "-" in some cases (cols J and K)

  rawRes <- readxl::read_xlsx("data/MM_Density_RES_Data.xlsx", sheet = "Narwhal (P2)", na = "-")
  
  
  resGrid <- read_csv("data/Narwhal.csv") %>%
    rename(RES = `Overall Probability`)
  
  # Coding of missing values are "-", several numeric fields interpreted as char
  
  narwhalData <- rawRes %>%
    mutate(across(c(`Uncertainty value`, `95% CI uncertainty low (per km2)`, `95% CI uncertainty high (per km2)`), as.numeric),
           `Uncertainty measure` = ifelse(`Uncertainty measure` == "-", NA, `Uncertainty measure`)) %>%
    rename(RES = `Species relative suitability index`, Density = `Density estimate (per km2)`,
           LowerCI = `95% CI uncertainty low (per km2)`,
           UpperCI = `95% CI uncertainty high (per km2)`) 

# Adding survey uncertainty -------------------------------------------------------------------
#' Here devise resampling for the different sorts of uncertainty that are present in the survey data
#' Note, there are single measures and upper/lower 95% CIS
#' 

  table(narwhalData$`Uncertainty measure`, useNA = "always")
  
  # check that where there is not a single measure, we do have a CI
  test <- narwhalData %>% filter(is.na(`Uncertainty measure`)) 
  any(is.na(test$UpperCI))
  
  narwhalData <- narwhalData %>% 
    mutate(`Uncertainty measure` = ifelse(`Uncertainty measure` == "NA", NA, `Uncertainty measure`)) %>% 
    filter(!(is.na(`Uncertainty measure`) & is.na(UpperCI))) 
  
  
  table(narwhalData$`Uncertainty measure`)
    
  # note some have CIs and alterative measure. Will use CIs when available
  # %-age CV looks to just be CV * 100
  # convert to CVs, extract as CIs so single sampling method required
  # everything in the data with a CV provided a CI as well
  
  narwhalData <- narwhalData %>%
    mutate(`Uncertainty value` = ifelse(str_detect(`Uncertainty measure`, "% CV for abundance"), `Uncertainty value`/100, `Uncertainty value`),
           `Uncertainty measure` = ifelse(str_detect(`Uncertainty measure`, "CV for abundance"), "CV", `Uncertainty measure`), 
           `Uncertainty value` = ifelse(str_detect(`Uncertainty measure`, "CV"), `Uncertainty value`*Density, `Uncertainty value`),
           `Uncertainty measure` = ifelse(str_detect(`Uncertainty measure`, "CV"), "SE", `Uncertainty measure`),
            workingSE = `Uncertainty value`,
            workingSE = ifelse(is.na(workingSE), SEfromCI(Density, LowerCI, UpperCI), workingSE)
           )
  
  # obtain CIs from lognormal using mean, SE
  
  narwhalData <- narwhalData %>%
      mutate(workingLowerCI = qlnorm(0.025, logMu(Density, workingSE), logSigma(Density, workingSE)),
             workingUpperCI = qlnorm(0.975, logMu(Density, workingSE), logSigma(Density, workingSE)),
             lowerError = LowerCI - workingLowerCI, 
             upperError = UpperCI - workingUpperCI,
             blockID = paste(narwhalData$`Survey ID`, narwhalData$`Area/Segment`, sep = ":")
      )

  # sample for each of the area/segments collectively (as not independent)  
  
  dataList <- split(narwhalData, narwhalData$blockID)
  
  set.seed(4835)
  
  sampleList <- lapply(dataList, function(q){sampleLN(q$Density[1], q$workingSE[1], 500)})
  
  sampleDF <- plyr::ldply(sampleList) %>%
    mutate(`Area/Segment` = str_extract(.id, "(?<=:).+"),
           `Survey ID` = str_extract(.id, ".+(?=:)")) %>%
    select(-.id)
    
  
  narwhalSamples <- narwhalData %>% left_join(sampleDF)
  
  
  narwhalSamples <- narwhalSamples %>% 
    select(-Density) %>%
    pivot_longer(names_to = "sampleID", values_to = "Density", V1:V500) 
  
  
  narwhalList <- split(narwhalSamples, narwhalSamples$sampleID)

  fittedList <- lapply(narwhalList, gamFit, inRES = data.frame(RES = seq(0, 1, by = 0.01))) 
  
  fittedDF <- fittedList %>% 
    bind_rows() 
  
  fittedMatrix <- matrix(fittedDF$Pred, ncol = 500)
  
  test <- t(apply(fittedMatrix, 1, function(q){quantile(q, probs = c(0.025, 0.5, 0.975))}))
  
  testSE <- apply(fittedMatrix, 1, sd)
  
  resFits <- as.data.frame(test) %>%
    mutate(RES = seq(0, 1, by = 0.01), SE = testSE) %>%
    rename(lower = `2.5%`, med = `50%`, upper = `97.5%`) %>%
    mutate(med = ifelse(med < 0, min(abs(med)), med),
            CV = SE/med) %>%
    arrange(RES)
  
  plottingDF <- resFits 
  
  ggplot(plottingDF) + 
    ggthemes::theme_fivethirtyeight() +
    geom_line(aes(RES, med), size = 2, alpha = 0.7, col = "purple") +
    geom_ribbon(aes(x = RES, ymin = lower, ymax = upper), fill = "purple", alpha = 0.2) + 
    ggtitle("Fitted function", "Narwhal: observed densities & monotone fit") 
  

# Create RES grid predictions -----------------------------------------------------------------

  resFits <- resFits %>% select(med, RES, CV) %>%
    rename(PredDensity = med) %>%
    mutate(RES = round(RES, 2))
  
  predictionOutput <- resGrid %>% left_join(resFits, by = "RES")
  
  predictionOutput <- predictionOutput %>%
    mutate(PredDensity = ifelse(RES == 0, 0, PredDensity),
           CV = ifelse(RES == 0, NA, CV))
   
  summary(predictionOutput)
  
  
  
  write.csv(predictionOutput, file = "data/predictions/narwhal_predictions.csv", row.names = F)

