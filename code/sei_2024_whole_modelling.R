# ---------------------------------------------------------------------------------------------
#' Analysis file for the sei whales data
#' 2024 - NE atlantic analyses
#'  - Location: whole
#' 


# Preamble ------------------------------------------------------------------------------------
library(tidyverse)
library(mgcv)
library(furrr)

source("code/tools.R")  
speciesName <- "Sei"  
locationName <- "whole"
nBoot <- 100

# note missing values in spreadsheet indicated by "-" in some cases (cols J and K)

rawRes <- read_csv("NE Atlantic Data/Sei_whale.csv") 


resGrid <- read_csv("NE Atlantic Data/RES spreadsheets/Sei whale (Balaenoptera borealis) - Native range.csv")

# Coding of missing values are "-", several numeric fields interpreted as char

modellingData <- rawRes %>%
  rename(UpperCI = CI_95_high,
         LowerCI = CI_95_low)




# Adding survey uncertainty -------------------------------------------------------------------
#' Here devise resampling for the different sorts of uncertainty that are present in the survey data
#' Note, there are single measures and upper/lower 95% CIS
#'

  table(modellingData$Uncertainty_measure, useNA = "always")

  # check that where there is not a single measure, we do have a CI
  test <- modellingData %>% filter(is.na(Uncertainty_measure))
  any(is.na(test$UpperCI))
  all(is.na(test$UpperCI))

  test %>% filter(!is.na(LowerCI))
  test %>% filter(is.na(LowerCI)) %>% summarise(sum(Density))

  #temp <- test %>% filter(is.na(LowerCI)); plot(temp$RES)

  modellingData <- modellingData %>%
    mutate(Uncertainty_measure = ifelse(is.na(Uncertainty_measure) & !is.na(LowerCI), "Confint", Uncertainty_measure)) %>%
    filter(!(is.na(Uncertainty_measure) & is.na(UpperCI)))

  table(modellingData$Uncertainty_measure, useNA = "always")


  # note some have CIs and alterative measure. Will use CIs when available
  # %-age CV looks to just be CV * 100
  # convert to CVs, extract as CIs so single sampling method required
  # everything in the data with a CV provided a CI as well

  modellingData <- modellingData %>%
    mutate(Uncertainty = ifelse(str_detect(Uncertainty_measure, "% CV abundance"), Uncertainty/100, Uncertainty),
           Uncertainty_measure = ifelse(str_detect(Uncertainty_measure, "% CV for abundance"), "CV", Uncertainty_measure),
           Uncertainty = ifelse(str_detect(Uncertainty_measure, "CV"), Uncertainty*Density, Uncertainty),
           Uncertainty_measure = ifelse(str_detect(Uncertainty_measure, "CV"), "SE", Uncertainty_measure),
           Uncertainty = ifelse(str_detect(Uncertainty_measure, "Variance of density"), sqrt(Uncertainty), Uncertainty),
           Uncertainty_measure = ifelse(str_detect(Uncertainty_measure, "Variance of density"), "SE", Uncertainty_measure),
           workingSE = Uncertainty,
           workingSE = ifelse(Uncertainty_measure == "Confint", SEfromCI(Density, LowerCI, UpperCI), workingSE)
           ) %>%
    filter(Uncertainty_measure != "SE abundance",
           !is.na(Uncertainty))

  # obtain CIs from lognormal using mean, SE

  modellingData <- modellingData %>%
      mutate(workingLowerCI = qlnorm(0.025, logMu(Density, workingSE), logSigma(Density, workingSE)),
             workingUpperCI = qlnorm(0.975, logMu(Density, workingSE), logSigma(Density, workingSE)),
             lowerError = LowerCI - workingLowerCI,
             upperError = UpperCI - workingUpperCI,
             blockID = Survey
      )

  # sample for each of the area/segments collectively (as not independent)

  dataList <- split(modellingData, modellingData$blockID)

  set.seed(4835)

  sampleList <- lapply(dataList, function(q){sampleLN(q$Density[1], q$workingSE[1], nBoot)})

  sampleDF <- plyr::ldply(sampleList) %>%
    rename(Survey = .id)

  speciesSamples <- modellingData %>% left_join(sampleDF)


  speciesSamples <- speciesSamples %>%
    select(-Density) %>%
    pivot_longer(names_to = "sampleID", values_to = "Density", V1:V100)


  sampleList <- split(speciesSamples, speciesSamples$sampleID)


  plan(multicore, workers = 10)


  fittedList <- future_map(sampleList, \(x) gamFitTweedie(x, inP = 1.2, inRES = data.frame(RES = seq(0, 1, by = 0.01))))


  #fittedList <- lapply(sampleList, gamFit, inRES = data.frame(RES = seq(0, 1, by = 0.01)))

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