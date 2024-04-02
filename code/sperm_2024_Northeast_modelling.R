# ---------------------------------------------------------------------------------------------
#' Analysis file for the sperm whales
#' 
#' 


# Preamble ------------------------------------------------------------------------------------

library(tidyverse)
library(mgcv)
library(furrr)
library(progressr)

source("code/tools.R")  
speciesName <- "Sperm"  
locationName <- "northeast"
nBoot <- 500

  # note missing values in spreadsheet indicated by "-" in some cases (cols J and K)

  rawRes <- read_csv("NE Atlantic Data/Sperm_whale.csv") %>% 
    filter(str_detect(Location, "Northeast"))
  
  
  resGrid <- read_csv("NE Atlantic Data/RES spreadsheets/Sperm whale (Physeter macrocephalus) - Native range.csv")
  
  # Coding of missing values are "-", several numeric fields interpreted as char
  
  spermData <- rawRes %>% 
    rename(UpperCI = CI_95_high,
           LowerCI = CI_95_low) 
    


# Adding survey uncertainty -------------------------------------------------------------------
#' Here devise resampling for the different sorts of uncertainty that are present in the survey data
#' Note, there are single measures and upper/lower 95% CIS
#' 

  table(spermData$Uncertainty_measure, useNA = "always")
  
  # check that where there is not a single measure, we do have a CI
  test <- spermData %>% filter(is.na(Uncertainty_measure)) 
  any(is.na(test$UpperCI))
  all(is.na(test$UpperCI))
  
  test %>% filter(!is.na(LowerCI))
  
  spermData <- spermData %>% 
    mutate(Uncertainty_measure = ifelse(is.na(Uncertainty_measure) & !is.na(LowerCI), "Confint", Uncertainty_measure)) %>%
    filter(!(is.na(Uncertainty_measure) & is.na(UpperCI))) 
  
  table(spermData$Uncertainty_measure, useNA = "always")
  
    
  # note some have CIs and alterative measure. Will use CIs when available
  # %-age CV looks to just be CV * 100
  # convert to CVs, extract as CIs so single sampling method required
  # everything in the data with a CV provided a CI as well
  
  spermData <- spermData %>%
    mutate(Uncertainty = ifelse(str_detect(Uncertainty_measure, "% CV abundance"), Uncertainty/100, Uncertainty), 
           Uncertainty_measure = ifelse(str_detect(Uncertainty_measure, "% CV for abundance"), "CV", Uncertainty_measure), 
           Uncertainty = ifelse(str_detect(Uncertainty_measure, "CV"), Uncertainty*Density, Uncertainty),
           Uncertainty_measure = ifelse(str_detect(Uncertainty_measure, "CV"), "SE", Uncertainty_measure),
           Uncertainty = ifelse(str_detect(Uncertainty_measure, "Variance of density"), sqrt(Uncertainty), Uncertainty),
           Uncertainty_measure = ifelse(str_detect(Uncertainty_measure, "Variance of density"), "SE", Uncertainty_measure),
           workingSE = Uncertainty,
           workingSE = ifelse(Uncertainty_measure == "Confint", SEfromCI(Density, LowerCI, UpperCI), workingSE)
           ) 
  
  # obtain CIs from lognormal using mean, SE
  
  spermData <- spermData %>%
      mutate(workingLowerCI = qlnorm(0.025, logMu(Density, workingSE), logSigma(Density, workingSE)),
             workingUpperCI = qlnorm(0.975, logMu(Density, workingSE), logSigma(Density, workingSE)),
             lowerError = LowerCI - workingLowerCI, 
             upperError = UpperCI - workingUpperCI,
             blockID = Survey
      ) 

  # sample for each of the area/segments collectively (as not independent)  
  
  dataList <- split(spermData, spermData$blockID)
  
  set.seed(4835)
  
  sampleList <- lapply(dataList, function(q){sampleLN(q$Density[1], q$workingSE[1], nBoot)})
  
  sampleDF <- plyr::ldply(sampleList) %>%
    rename(Survey = .id)
  
  spermSamples <- spermData %>% left_join(sampleDF)
  
  
  spermSamples <- spermSamples %>% 
    select(-Density) %>%
    pivot_longer(names_to = "sampleID", values_to = "Density", V1:last_col()) 
  
  
  spermList <- split(spermSamples, spermSamples$sampleID)
  
  
  # fittedList <- lapply(spermList, gamFit, inRES = data.frame(RES = seq(0, 1, by = 0.01))) 
  
  plan(multicore, workers = 10)
  
  testFN <- function(inList){
  
  p <- progressor(steps = length(inList))
  
  future_map(inList, \(x) {p(); gamFitTweedie(x, inP = 1.2, inRES = data.frame(RES = seq(0, 1, by = 0.01)))})
  
  }
  
  
  with_progress({
    fittedList <- testFN(spermList)
  })
  
  
  fittedDF <- fittedList %>%
    bind_rows()
  
  fittedMatrix <- matrix(fittedDF$Pred, ncol = nBoot)
  
  test <- t(apply(fittedMatrix, 1, function(q){quantile(q, probs = c(0.025, 0.5, 0.975))}))
  
  testSE <- apply(fittedMatrix, 1, sd)
  
  resFits <- as.data.frame(test) %>%
    mutate(RES = seq(0, 1, by = 0.01), SE = testSE) %>%
    rename(lower = `2.5%`, med = `50%`, upper = `97.5%`) %>%
    mutate(med = ifelse(med < 0, min(abs(med)), med),
           CV = SE/med) %>%
    #CV = ifelse(CV > 2, 2, CV)) %>%
    arrange(RES)
  
  plottingDF <- resFits
  
  bootPlot <- ggplot(plottingDF) +
    ggthemes::theme_fivethirtyeight() +
    #geom_point(data = ribbonsealData, aes(RES, Density), size = 2, alpha = 0.2) +
    geom_line(aes(RES, med), size = 2, alpha = 0.7, col = "purple") +
    geom_ribbon(aes(x = RES, ymin = lower, ymax = upper), fill = "purple", alpha = 0.2) +
    ggtitle("Density as function of RES", paste0(speciesName, " whale: bootstrapped monotone spline fits"))
  
  bootPlot
  
  ggsave(paste0("docs/images/", speciesName, "_", locationName, "_bootplot_", nBoot, ".png"), units = "cm", width = 30, height = 20)
  saveRDS(plottingDF, paste0("data/plotting components/", speciesName, "_", locationName, "_plotElements",".rds"))
  
  
  
  # Create RES grid predictions -----------------------------------------------------------------
  
  resFits <- resFits %>% select(med, RES, CV) %>%
    rename(PredDensity = med) %>%
    mutate(RES = round(RES, 2))
  
  predictionOutput <- resGrid %>% left_join(resFits, by = "RES")
  
  predictionOutput <- predictionOutput %>%
    mutate(PredDensity = ifelse(RES == 0, 0, PredDensity),
           CV = ifelse(RES == 0, NA, CV))
  
  summary(predictionOutput)
  
  
  
  write.csv(predictionOutput, file = paste0("data/predictions/", speciesName, "_", locationName, "_predictions.csv"), row.names = F)
  