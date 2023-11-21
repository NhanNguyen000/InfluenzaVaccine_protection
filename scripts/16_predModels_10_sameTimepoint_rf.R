rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)

get.ablevel <- function(input) {
  # Aim: convert the ab titer level into low (< 40) and high (>= 40) groups 
  outcome = factor(ifelse(input < 40, "L", ifelse(input >= 40, "H", input)), 
                   levels = c("L", "H"))
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
  left_join(cohorts$donorSample_all %>% filter(time %in% c("d0", "d28"))) %>%
  filter(condition == "Healthy") %>%
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

## prepare train-test, season 2015, d0--------------------------------
dat_2015_d0 <- metadata_healthy %>% left_join(inputDat) %>% 
  filter(season == "2015", time == "d0")

trainSet <- list()
trainSet$H1N1 <- dat_2015_d0 %>% mutate(ab = get.ablevel(H1N1_d0))
trainSet$H3N2 <- dat_2015_d0 %>% mutate(ab = get.ablevel(H3N2_d0))
trainSet$B <- dat_2015_d0 %>% mutate(ab = get.ablevel(B_d0))

trainSet %>% lapply(function(x) x%>% count(ab))

## prepare validation-test 1, season 2015, d28 --------------------------------
dat_2015_d28 <- metadata_healthy %>% left_join(inputDat) %>% 
  filter(season == "2015", time == "d28")

valiSet_2015 <- list()
valiSet_2015$H1N1 <- dat_2015_d28 %>% mutate(ab = get.ablevel(H1N1_d28))
valiSet_2015$H3N2 <- dat_2015_d28 %>% mutate(ab = get.ablevel(H3N2_d28))
valiSet_2015$B <- dat_2015_d28 %>% mutate(ab = get.ablevel(B_d28))

valiSet_2015 %>% lapply(function(x) x%>% count(ab))

## prepare validation-test, all 4 season --------------------------------
valiSets <- list()
for (timepoint in c("d0", "d28")) {
  dat_temp <- metadata_healthy %>% left_join(inputDat) %>% 
    filter(time == timepoint)
  
  # for H1N1
  for (year in c("2014", "2015", "2019", "2020")) {
    valiSets[[paste0("H1N1_", year, "_", timepoint)]] <- dat_temp %>% 
      filter(season == year) %>% 
      rename("ab_strain" =  paste0("H1N1_", timepoint)) %>% mutate(ab = get.ablevel(ab_strain))
  }
  
  # for H3N2
  for (year in c("2014", "2015", "2019", "2020")) {
    valiSets[[paste0("H3N2_", year, "_", timepoint)]] <- dat_temp %>% 
      filter(season == year) %>% 
      rename("ab_strain" =  paste0("H3N2_", timepoint)) %>% mutate(ab = get.ablevel(ab_strain))
  }
  
  # for B
  for (year in c("2014", "2015")) {
    valiSets[[paste0("B_", year, "_", timepoint)]] <- dat_temp %>% 
      filter(season == year) %>% 
      rename("ab_strain" =  paste0("B_", timepoint)) %>% mutate(ab = get.ablevel(ab_strain))
  }
  
  # for Bvictoria
  for (year in c("2019", "2020")) {
    valiSets[[paste0("Bvictoria_", year, "_", timepoint)]] <- dat_temp %>% 
      filter(season == year) %>% 
      rename("ab_strain" =  paste0("Bvictoria_", timepoint)) %>% mutate(ab = get.ablevel(ab_strain))
  }
  
  # for Byamagata
  for (year in c("2019", "2020")) {
    valiSets[[paste0("Byamagata_", year, "_", timepoint)]] <- dat_temp %>% 
      filter(season == year) %>% 
      rename("ab_strain" =  paste0("Byamagata_", timepoint)) %>% mutate(ab = get.ablevel(ab_strain))
  }
}

valiSets %>% lapply(function(x) x%>% count(ab))

## input variable --------------------------------
proNames <- colnames(proteinDat_impute$iMED_2014)
meboNames <- colnames(mebo_Dat$ZirFlu_2019)

inputVars <- c("age", "sex", proNames, meboNames)

# prepare input and validation sets
train_inputSet <- trainSet %>% lapply(function(x) x[, inputVars])
train_predOut <- trainSet %>% lapply(function(x) x$ab)

validation_inputSets <- valiSets %>% lapply(function(x) x[, inputVars])
validation_predOuts <- valiSets %>% lapply(function(x) x$ab)

# Run models (run in Slum/server) ==========================================
set.seed(123) #set the seed to make your partition reproducible

#selected_metric = 'Accuracy' # for classification
selected_metric =  "Kappa" # similar to classification accuracy but it is useful to normalize the imbalance in classes
selected_crossVal <- trainControl(method="repeatedcv", number = 5, 
                                  repeats = 10)

method_list <- c("rf", "svmRadial", "nb", "knn", "nnet")

nreps <- 10
for (selected_method in method_list) {
  netFits <- list()
  netFit_pred <- list()
  for (strain in names(train_inputSet)) {
    print(paste0("Strain ", strain, " 2015"))
    for (i in 1:nreps) {
      netFits[[strain]][[i]] <- train(x = train_inputSet[[strain]], y = train_predOut[[strain]],
                            method = selected_method, metric = selected_metric,
                            tuneLength =  30,
                            trControl = selected_crossVal)
      
      for (valiSet in names(validation_inputSets)) {
        val.pred_temp <- predict(netFits[[strain]][[i]], newdata = validation_inputSets[[valiSet]], type = 'raw')
        netFit_pred[[strain]][[valiSet]][[i]] <- confusionMatrix(val.pred_temp, validation_predOuts[[valiSet]]) # Training accuracy
      }
      print(paste0("Run ", i, " time."))
    }
  }
  save(netFits, netFit_pred, file = paste0("res.predModel_sameTime", method, ".RData"))
}