# ---------------------------------------------------------------------------------------------
#' Analysis file for the Ribbon Seal
#' 
#' 


# Preamble ------------------------------------------------------------------------------------

  library(tidyverse)
  library(mgcv)
  
  source("code/tools.R")  
  
  # note missing values in spreadsheet indicated by "-" in some cases (cols J and K)

  rawRes <- readxl::read_xlsx("data/MM_Density_RES_Data.xlsx", sheet = "Ribbon seal (P2)", na = "-")
  
  
  resGrid <- read_csv("data/Ribbon_seal.csv") %>%
    rename(RES = `Overall Probability`)
  
  # Coding of missing values are "-", several numeric fields interpreted as char
  
  ribbonsealData <- rawRes %>%
    mutate(across(c(`Uncertainty value`, `95% CI uncertainty low (per km2)`, `95% CI uncertainty high (per km2)`), as.numeric),
           `Uncertainty measure` = ifelse(`Uncertainty measure` == "-", NA, `Uncertainty measure`)) %>%
    rename(RES = `Species relative suitability index`, Density = `Density estimate (per km2)`,
           LowerCI = `95% CI uncertainty low (per km2)`,
           UpperCI = `95% CI uncertainty high (per km2)`) 


# Adding survey uncertainty -------------------------------------------------------------------
#' Here devise resampling for the different sorts of uncertainty that are present in the survey data
#' Note, there are single measures and upper/lower 95% CIS
#' 

  table(ribbonsealData$`Uncertainty measure`, useNA = "always")
  
  # check that where there is not a single measure, we do have a CI
  test <- ribbonsealData %>% filter(is.na(`Uncertainty measure`)) 
  any(is.na(test$UpperCI))
  
  ribbonsealData <- ribbonsealData %>% 
    mutate(`Uncertainty measure` = ifelse(`Uncertainty measure` == "NA", NA, `Uncertainty measure`)) %>%
    filter(!(is.na(`Uncertainty measure`) & is.na(UpperCI))) 
    
  # note some have CIs and alterative measure. Will use CIs when available
  # %-age CV looks to just be CV * 100
  # convert to CVs, extract as CIs so single sampling method required
  # everything in the data with a CV provided a CI as well
  # 117,000 abund, SE 25,000 gives 25/117 of density of 0.152519704 = 0.03258968
  # 38,000 abund, SE 4,000 gives 4/38 of desnity of 0.049536314 = 0.005214349
  # 
  
  ribbonsealData <- ribbonsealData %>%
    mutate(#`Uncertainty value` = ifelse(str_detect(`Uncertainty measure`, "% CV for abundance"), `Uncertainty value`/100, `Uncertainty value`), 
            `Uncertainty value` = ifelse(`Uncertainty value` == 25000, 0.03258968, `Uncertainty value`), # hacks as need abundance to convert
            `Uncertainty value` = ifelse(`Uncertainty value` == 4000, 0.005214349, `Uncertainty value`), # hacks as need abundance to convert
            `Uncertainty measure` = ifelse(str_detect(`Uncertainty measure`, "for abundance"), "SE", `Uncertainty measure`),
            workingSE = `Uncertainty value`,
            workingSE = ifelse(is.na(workingSE), SEfromCI(Density, LowerCI, UpperCI), workingSE)
           )
  
  # obtain CIs from lognormal using mean, SE
  
  ribbonsealData <- ribbonsealData %>%
      mutate(workingLowerCI = qlnorm(0.025, logMu(Density, workingSE), logSigma(Density, workingSE)),
             workingUpperCI = qlnorm(0.975, logMu(Density, workingSE), logSigma(Density, workingSE)),
             lowerError = LowerCI - workingLowerCI, 
             upperError = UpperCI - workingUpperCI,
             blockID = paste(ribbonsealData$`Survey ID`, ribbonsealData$`Area/Segment`, sep = ":")
      )

  # sample for each of the area/segments collectively (as not independent)  
  
  dataList <- split(ribbonsealData, ribbonsealData$blockID)
  
  set.seed(4835)
  
  sampleList <- lapply(dataList, function(q){sampleLN(q$Density[1], q$workingSE[1], 500)})
  
  sampleDF <- plyr::ldply(sampleList) %>%
    mutate(`Area/Segment` = str_extract(.id, "(?<=:).+"),
           `Survey ID` = str_extract(.id, ".+(?=:)")) %>%
    select(-.id)
    
  
  ribbonsealSamples <- ribbonsealData %>% left_join(sampleDF)
  
  
  ribbonsealSamples <- ribbonsealSamples %>% 
    select(-Density) %>%
    pivot_longer(names_to = "sampleID", values_to = "Density", V1:V500) 
  
  
  ribbonsealList <- split(ribbonsealSamples, ribbonsealSamples$sampleID)

  fittedList <- lapply(ribbonsealList, gamFit, inRES = data.frame(RES = seq(0, 1, by = 0.01))) 
  
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
  
  bootPlot <- ggplot(plottingDF) + 
    ggthemes::theme_fivethirtyeight() +
     geom_line(aes(RES, med), size = 2, alpha = 0.7, col = "purple") +
    geom_ribbon(aes(x = RES, ymin = lower, ymax = upper), fill = "purple", alpha = 0.2) + 
    ggtitle("Density as function of RES", "Ribbon seal: bootstrapped monotone spline fits") 
  
  bootPlot
  
  ggsave("docs/images/ribbonSeal_bootplot.png", units = "cm", width = 30, height = 20)
  
# Create RES grid predictions -----------------------------------------------------------------

  resFits <- resFits %>% select(med, RES, CV) %>%
    rename(PredDensity = med) %>%
    mutate(RES = round(RES, 2))
  
  predictionOutput <- resGrid %>% left_join(resFits, by = "RES")
  
  predictionOutput <- predictionOutput %>%
    mutate(PredDensity = ifelse(RES == 0, 0, PredDensity),
           CV = ifelse(RES == 0, NA, CV))
   
  summary(predictionOutput)
  
  
  
  write.csv(predictionOutput, file = "data/predictions/ribbonseal_predictions.csv", row.names = F)

