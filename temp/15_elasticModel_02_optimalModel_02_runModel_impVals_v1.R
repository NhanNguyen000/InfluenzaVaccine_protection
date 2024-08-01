rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}

# check NAs percentage
get.rmProteins <- function(dat, NA_cutoff) {
  # Description: identify variable have more missing value (NA) than the NA threshold (NA cutoff)
  
  # Arguments: 
  # dat - a numeric data table with sample in row (sample name is rowname), variable per column
  # selected_sample - name of the selected samples
  
  # Returns: variable name which have more NA than the NA_cutoff
  rmProteins <- which(colSums(is.na(dat))/nrow(dat) > NA_cutoff)
  return(rmProteins)
}

# load data =======================================================================
load("/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~convert_protectees(.x)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "protectee"))) %>%
  mutate(sex = ifelse(sex == "m", 0, 1))

# impute protein data --------------------------------
cutoff <- 0.2
rmProteins <- list()
for (season in names(protein_Dat)) {
  rmProteins[[season]] <- get.rmProteins(protein_Dat[[season]], NA_cutoff = cutoff)
} 
rmProteins # season 2019 and 2020, 14 proteins have a lot of NA --> can use 292/306 protein
rmProtein_names <- unique(rmProteins %>% lapply(function(x) names(x)) %>% unlist())

proteinDat_impute <- list()
for (season in names(protein_Dat)) {
  proteinDat_impute[[season]] <- protein_Dat[[season]] %>%
    select(-rmProtein_names) %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))
}

# prepare input data ------------------------------------------------------
inputDat <- proteinDat_impute %>% purrr::reduce(rbind) %>% 
  cbind(mebo_Dat %>% purrr::reduce(rbind)) %>%
  as.data.frame %>% rownames_to_column("name")

dat_temp <- metadata_healthy %>% left_join(inputDat)

dat_temp %>% count(season, H1N1_reclassify) # LL group size >= 3 in all 4 seasons -> can use all 4 seasons
dat_temp %>% count(season, H3N2_reclassify) # LL group size >= 3 in season 2015-> can use season 2015
dat_temp %>% count(season, B_reclassify) # LL group size >= 3 in all 2 seasons 2014 and 2015 -> can use all 2 seasons
dat_temp %>% count(season, Bvictoria_reclassify) # LL group size < 3 in all 2 seasons 2019 and 2020 --> can not use 
dat_temp %>% count(season, Byamagata_reclassify) # LL group size >= 3 in seasons 2020 --> can use season 2020


# Train elastic model --------------------------------
set.seed(123) #set the seed to make your partition reproducible

## prepare train-test --------------------------------
trainSet <- dat_temp %>% filter(season == "2015") %>% 
  rename("reclassify" = "H1N1_reclassify", "ab_T1" = "H1N1_T1")
trainSet %>% count(reclassify)

## input variable --------------------------------
load("impVars_elasticModels.RData")

inputVariables <- c("age", "sex", impVars)

# prepare input and validation sets
train_inputSet <- trainSet[, inputVariables]
train_predOut <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()

# train the model -------------------------------------
netFit_models_impVars <- list()
for (i in c(1:10)) {
  netFit_models_impVars[[i]] <- train(x = train_inputSet,
                                          y = train_predOut,
                                          method = "glmnet", 
                                          #metric = 'Accuracy', # for classification
                                          metric = "Kappa", # similar to classification accuracy but it is useful to normalize the imbalance in classes
                                          #metric = "ROC",
                                          tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                                                               .lambda = seq(0,1,by=0.01)),
                                          trControl = trainControl(method="repeatedcv",
                                                                   number=5,
                                                                   repeats=20,
                                                                   summaryFunction=twoClassSummary, 
                                                                   classProbs=T,
                                                                   savePredictions = T),
  )
}

save(netFit_models_impVars, file = "netFit_models_impVars.RData")


