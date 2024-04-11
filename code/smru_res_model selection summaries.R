library(tidyverse)

# fin -------------------------------------------------------------------------------------------------------------

speciesData <- read_csv("NE Atlantic Data/fin_whale.csv") %>% 
  select(Commonname, RES, Density, Study_area)

whole <- read_csv("data/predictions/fin_whole_predictions.csv") %>% 
  select(RES, PredDensity, CV) %>% 
  distinct()


atlantic <- read_csv("data/predictions/fin_Atlantic_predictions.csv") %>% 
  select(RES, PredDensity, CV) %>% 
  distinct()


north <- read_csv("data/predictions/fin_North_predictions.csv") %>% 
  select(RES, PredDensity, CV) %>% 
  distinct()


northEast <- read_csv("data/predictions/fin_northeast_predictions.csv") %>% 
  select(RES, PredDensity, CV) %>% 
  distinct()

compData <- speciesData %>% 
  left_join(whole, by = "RES") %>% 
  rename(whole_pred = PredDensity, whole_CV = CV) %>% 
  left_join(atlantic, by = "RES") %>% 
  rename(atlantic_pred = PredDensity, atlantic_CV = CV) %>% 
  left_join(northEast, by = "RES") %>% 
  rename(NE_pred = PredDensity, NE_CV = CV) %>% 
  left_join(north, by = "RES") %>% 
  filter(Study_area == "Y") %>% 
  rename(north_pred = PredDensity, north_CV = CV) %>% 
  mutate(wholeDiff = (whole_pred - Density)^2,
         atlanticDiff = (atlantic_pred - Density)^2,
         northDiff = (north_pred - Density)^2,
         NEDiff = (NE_pred - Density)^2) %>% 
  select(-contains("pred")) %>% 
  pivot_longer(names_to = "source", values_to = "value", whole_CV:NEDiff) 

  finSummaryData <- compData %>% 
    filter(RES != 0) %>% 
    group_by(source) %>% 
    summarise(Commonname = Commonname[1], meanMetric = mean(value)) %>% 
    mutate(region = str_remove(source, "Diff|_CV"),
           type = if_else(str_detect(source, "Diff"), "SSE", "CV")) %>% 
    select(-source) %>% 
    pivot_wider(names_from = "type", values_from = meanMetric)


# blue -------------------------------------------------------------------------------------------------------------
  
  speciesData <- read_csv("NE Atlantic Data/blue_whale.csv") %>% 
    select(Commonname, RES, Density, Study_area)
  
  whole <- read_csv("data/predictions/blue_whole_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  atlantic <- read_csv("data/predictions/blue_Atlantic_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  north <- read_csv("data/predictions/blue_North_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  northEast <- read_csv("data/predictions/blue_northeast_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  compData <- speciesData %>% 
    left_join(whole, by = "RES") %>% 
    rename(whole_pred = PredDensity, whole_CV = CV) %>% 
    left_join(atlantic, by = "RES") %>% 
    rename(atlantic_pred = PredDensity, atlantic_CV = CV) %>% 
    left_join(northEast, by = "RES") %>% 
    rename(NE_pred = PredDensity, NE_CV = CV) %>% 
    left_join(north, by = "RES") %>% 
    filter(Study_area == "Y") %>% 
    rename(north_pred = PredDensity, north_CV = CV) %>% 
    mutate(wholeDiff = (whole_pred - Density)^2,
           atlanticDiff = (atlantic_pred - Density)^2,
           northDiff = (north_pred - Density)^2,
           NEDiff = (NE_pred - Density)^2) %>% 
    select(-contains("pred")) %>% 
    pivot_longer(names_to = "source", values_to = "value", whole_CV:NEDiff) 
  
  blueSummaryData <- compData %>% 
    filter(RES != 0) %>% 
    group_by(source) %>% 
    summarise(Commonname = Commonname[1], meanMetric = mean(value)) %>% 
    mutate(region = str_remove(source, "Diff|_CV"),
           type = if_else(str_detect(source, "Diff"), "SSE", "CV")) %>% 
    select(-source) %>% 
    pivot_wider(names_from = "type", values_from = meanMetric)
  
  


# bowhead -------------------------------------------------------------------------------------------------------------
#' Note there are no data in the study region
 
  speciesData <- read_csv("NE Atlantic Data/bowhead_whale.csv") %>% 
    select(Commonname, RES, Density, Study_area)
  
  whole <- read_csv("data/predictions/bow_whole_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  atlantic <- read_csv("data/predictions/Bow_atlantic_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  north <- read_csv("data/predictions/Bow_North_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  northEast <- read_csv("data/predictions/Bow_Northeast_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  compData <- speciesData %>% 
    left_join(whole, by = "RES") %>% 
    rename(whole_pred = PredDensity, whole_CV = CV) %>% 
    left_join(atlantic, by = "RES") %>% 
    rename(atlantic_pred = PredDensity, atlantic_CV = CV) %>% 
    left_join(northEast, by = "RES") %>% 
    rename(NE_pred = PredDensity, NE_CV = CV) %>% 
    left_join(north, by = "RES") %>% 
    #filter(Study_area == "Y") %>% 
    rename(north_pred = PredDensity, north_CV = CV) %>% 
    mutate(wholeDiff = (whole_pred - Density)^2,
           atlanticDiff = (atlantic_pred - Density)^2,
           northDiff = (north_pred - Density)^2,
           NEDiff = (NE_pred - Density)^2) %>% 
    select(-contains("pred")) %>% 
    pivot_longer(names_to = "source", values_to = "value", whole_CV:NEDiff) 
  
  bowheadSummaryData <- compData %>% 
    filter(RES != 0) %>% 
    group_by(source) %>% 
    summarise(Commonname = Commonname[1], meanMetric = mean(value)) %>% 
    mutate(region = str_remove(source, "Diff|_CV"),
           type = if_else(str_detect(source, "Diff"), "SSE", "CV")) %>% 
    select(-source) %>% 
    pivot_wider(names_from = "type", values_from = meanMetric)
  
  
# brydes -------------------------------------------------------------------------------------------------------------
#' note no north east or data in region

  speciesData <- read_csv("NE Atlantic Data/brydes_whale.csv") %>% 
    select(Commonname, RES, Density, Study_area)
  
  whole <- read_csv("data/predictions/Bryde_Whole_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
 
  
  compData <- speciesData %>% 
    left_join(whole, by = "RES") %>% 
    rename(whole_pred = PredDensity, whole_CV = CV) %>% 
    mutate(wholeDiff = (whole_pred - Density)^2) %>% 
    select(-contains("pred")) %>% 
    pivot_longer(names_to = "source", values_to = "value", whole_CV:wholeDiff) 
  
  brydesSummaryData <- compData %>% 
    filter(RES != 0) %>% 
    group_by(source) %>% 
    summarise(Commonname = Commonname[1], meanMetric = mean(value)) %>% 
    mutate(region = str_remove(source, "Diff|_CV"),
           type = if_else(str_detect(source, "Diff"), "SSE", "CV")) %>% 
    select(-source) %>% 
    pivot_wider(names_from = "type", values_from = meanMetric)
  
  
# humpback -------------------------------------------------------------------------------------------------------------
  
  speciesData <- read_csv("NE Atlantic Data/humpback_whale.csv") %>% 
    select(Commonname, RES, Density, Study_area)
  
  whole <- read_csv("data/predictions/Hump_Whole_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  atlantic <- read_csv("data/predictions/Hump_Atlantic_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  north <- read_csv("data/predictions/Hump_North_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  northEast <- read_csv("data/predictions/Hump_Northeast_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  compData <- speciesData %>% 
    left_join(whole, by = "RES") %>% 
    rename(whole_pred = PredDensity, whole_CV = CV) %>% 
    left_join(atlantic, by = "RES") %>% 
    rename(atlantic_pred = PredDensity, atlantic_CV = CV) %>% 
    left_join(northEast, by = "RES") %>% 
    rename(NE_pred = PredDensity, NE_CV = CV) %>% 
    left_join(north, by = "RES") %>% 
    filter(Study_area == "Y") %>% 
    rename(north_pred = PredDensity, north_CV = CV) %>% 
    mutate(wholeDiff = (whole_pred - Density)^2,
           atlanticDiff = (atlantic_pred - Density)^2,
           northDiff = (north_pred - Density)^2,
           NEDiff = (NE_pred - Density)^2) %>% 
    select(-contains("pred")) %>% 
    pivot_longer(names_to = "source", values_to = "value", whole_CV:NEDiff) 
  
  humpbackSummaryData <- compData %>% 
    filter(RES != 0) %>% 
    group_by(source) %>% 
    summarise(Commonname = Commonname[1], meanMetric = mean(value)) %>% 
    mutate(region = str_remove(source, "Diff|_CV"),
           type = if_else(str_detect(source, "Diff"), "SSE", "CV")) %>% 
    select(-source) %>% 
    pivot_wider(names_from = "type", values_from = meanMetric)
  
  
# minke -------------------------------------------------------------------------------------------------------------
  
  speciesData <- read_csv("NE Atlantic Data/minke_whale.csv") %>% 
    select(Commonname, RES, Density, Study_area)
  
  whole <- read_csv("data/predictions/minke_whole_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  atlantic <- read_csv("data/predictions/minke_Atlantic_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  north <- read_csv("data/predictions/minke_North_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  northEast <- read_csv("data/predictions/minke_northeast_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  compData <- speciesData %>% 
    left_join(whole, by = "RES") %>% 
    rename(whole_pred = PredDensity, whole_CV = CV) %>% 
    left_join(atlantic, by = "RES") %>% 
    rename(atlantic_pred = PredDensity, atlantic_CV = CV) %>% 
    left_join(northEast, by = "RES") %>% 
    rename(NE_pred = PredDensity, NE_CV = CV) %>% 
    left_join(north, by = "RES") %>% 
    filter(Study_area == "Y") %>% 
    rename(north_pred = PredDensity, north_CV = CV) %>% 
    mutate(wholeDiff = (whole_pred - Density)^2,
           atlanticDiff = (atlantic_pred - Density)^2,
           northDiff = (north_pred - Density)^2,
           NEDiff = (NE_pred - Density)^2) %>% 
    select(-contains("pred")) %>% 
    pivot_longer(names_to = "source", values_to = "value", whole_CV:NEDiff) 
  
  minkeSummaryData <- compData %>% 
    filter(RES != 0) %>% 
    group_by(source) %>% 
    summarise(Commonname = Commonname[1], meanMetric = mean(value)) %>% 
    mutate(region = str_remove(source, "Diff|_CV"),
           type = if_else(str_detect(source, "Diff"), "SSE", "CV")) %>% 
    select(-source) %>% 
    pivot_wider(names_from = "type", values_from = meanMetric)
  
  
# North atlantic right whale -------------------------------------------------------------------------------------------------------------
#' Note only north, and none in region

  speciesData <- read_csv("NE Atlantic Data/North_Atlantic_right_whale.csv") %>% 
    select(Commonname, RES, Density, Study_area)
  
 
  north <- read_csv("data/predictions/Right_north_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  compData <- speciesData %>% 
    left_join(north, by = "RES") %>% 
    #filter(Study_area == "Y") %>% 
    rename(north_pred = PredDensity, north_CV = CV) %>% 
    mutate(northDiff = (north_pred - Density)^2) %>% 
    select(-contains("pred")) %>% 
    pivot_longer(names_to = "source", values_to = "value", north_CV:northDiff) 
  
  northRightSummaryData <- compData %>% 
    filter(RES != 0) %>% 
    group_by(source) %>% 
    summarise(Commonname = Commonname[1], meanMetric = mean(value)) %>% 
    mutate(region = str_remove(source, "Diff|_CV"),
           type = if_else(str_detect(source, "Diff"), "SSE", "CV")) %>% 
    select(-source) %>% 
    pivot_wider(names_from = "type", values_from = meanMetric)
  
  
# sei -------------------------------------------------------------------------------------------------------------
  
  speciesData <- read_csv("NE Atlantic Data/sei_whale.csv") %>% 
    select(Commonname, RES, Density, Study_area)
  
  whole <- read_csv("data/predictions/sei_whole_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  atlantic <- read_csv("data/predictions/sei_Atlantic_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  north <- read_csv("data/predictions/sei_North_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  northEast <- read_csv("data/predictions/sei_northeast_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  compData <- speciesData %>% 
    left_join(whole, by = "RES") %>% 
    rename(whole_pred = PredDensity, whole_CV = CV) %>% 
    left_join(atlantic, by = "RES") %>% 
    rename(atlantic_pred = PredDensity, atlantic_CV = CV) %>% 
    left_join(northEast, by = "RES") %>% 
    rename(NE_pred = PredDensity, NE_CV = CV) %>% 
    left_join(north, by = "RES") %>% 
    filter(Study_area == "Y") %>% 
    rename(north_pred = PredDensity, north_CV = CV) %>% 
    mutate(wholeDiff = (whole_pred - Density)^2,
           atlanticDiff = (atlantic_pred - Density)^2,
           northDiff = (north_pred - Density)^2,
           NEDiff = (NE_pred - Density)^2) %>% 
    select(-contains("pred")) %>% 
    pivot_longer(names_to = "source", values_to = "value", whole_CV:NEDiff) 
  
  seiSummaryData <- compData %>% 
    filter(RES != 0) %>% 
    group_by(source) %>% 
    summarise(Commonname = Commonname[1], meanMetric = mean(value)) %>% 
    mutate(region = str_remove(source, "Diff|_CV"),
           type = if_else(str_detect(source, "Diff"), "SSE", "CV")) %>% 
    select(-source) %>% 
    pivot_wider(names_from = "type", values_from = meanMetric)
  
  
# sperm -------------------------------------------------------------------------------------------------------------
  
  speciesData <- read_csv("NE Atlantic Data/sperm_whale.csv") %>% 
    select(Commonname, RES, Density, Study_area)
  
  whole <- read_csv("data/predictions/sperm_whole_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  atlantic <- read_csv("data/predictions/sperm_Atlantic_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  north <- read_csv("data/predictions/sperm_North_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  
  northEast <- read_csv("data/predictions/sperm_northeast_predictions.csv") %>% 
    select(RES, PredDensity, CV) %>% 
    distinct()
  
  compData <- speciesData %>% 
    left_join(whole, by = "RES") %>% 
    rename(whole_pred = PredDensity, whole_CV = CV) %>% 
    left_join(atlantic, by = "RES") %>% 
    rename(atlantic_pred = PredDensity, atlantic_CV = CV) %>% 
    left_join(northEast, by = "RES") %>% 
    rename(NE_pred = PredDensity, NE_CV = CV) %>% 
    left_join(north, by = "RES") %>% 
    filter(Study_area == "Y") %>% 
    rename(north_pred = PredDensity, north_CV = CV) %>% 
    mutate(wholeDiff = (whole_pred - Density)^2,
           atlanticDiff = (atlantic_pred - Density)^2,
           northDiff = (north_pred - Density)^2,
           NEDiff = (NE_pred - Density)^2) %>% 
    select(-contains("pred")) %>% 
    pivot_longer(names_to = "source", values_to = "value", whole_CV:NEDiff) 
  
  spermSummaryData <- compData %>% 
    filter(RES != 0) %>% 
    group_by(source) %>% 
    summarise(Commonname = Commonname[1], meanMetric = mean(value)) %>% 
    mutate(region = str_remove(source, "Diff|_CV"),
           type = if_else(str_detect(source, "Diff"), "SSE", "CV")) %>% 
    select(-source) %>% 
    pivot_wider(names_from = "type", values_from = meanMetric)
  
  
  

# Combine all -----------------------------------------------------------------------------------------------------

outputTable <- finSummaryData %>% 
    bind_rows(blueSummaryData, bowheadSummaryData, brydesSummaryData, 
              humpbackSummaryData, minkeSummaryData, northRightSummaryData, 
              seiSummaryData, spermSummaryData) %>% 
    rename(dataBasis = region)

writexl::write_xlsx(outputTable, "docs/NEModelling_modelselection.xlsx")
