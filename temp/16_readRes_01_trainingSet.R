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
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~convert_protectees(.x)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "protectee"))) %>%
  mutate(sex = ifelse(sex == "m", 0, 1))

# Run models (run in Slum/server) ==========================================
set.seed(123) #set the seed to make your partition reproducible

netFit_pred_trainSet <- list()
nreps <- 10

## Run prediction models: "knn", "nnet", "rf_tuneLength", "rf", "smv" --------------------------------
### impute protein data --------------------------------
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

## prepare input data ------------------------------------------------------
inputDat <- proteinDat_impute %>% purrr::reduce(rbind) %>% 
  cbind(mebo_Dat %>% purrr::reduce(rbind)) %>%
  as.data.frame %>% rownames_to_column("name")

dat_temp <- metadata_healthy %>% left_join(inputDat)

trainSet <- dat_temp %>% filter(season == "2015") %>% 
  rename("reclassify" = "H1N1_reclassify", "ab_d0" = "H1N1_d0")
trainSet %>% count(reclassify)

### input variable --------------------------------
proNames <- colnames(proteinDat_impute$iMED_2014)
meboNames <- colnames(mebo_Dat$ZirFlu_2019)

inputVars <- c("age", "sex", "ab_d0", proNames, meboNames)

# prepare input and validation sets
train_inputSet <- trainSet[, inputVars]
train_predOut <- trainSet[, c("reclassify")]

### run the models --------------------------------
models <- c("knn", "nnet", "glmnet", "rf_tuneLength", "rf", "smv")
for (nameModel in models) {
  load(paste0("res.predModel_", nameModel, ".RData"))
  
  for (i in 1:nreps) {
    val.pred_temp <- predict(netFits[[i]], 
                             newdata = train_inputSet, 
                             type = 'raw')
    
    netFit_pred_trainSet[[nameModel]][[i]] <- confusionMatrix(val.pred_temp, train_predOut) # Training accuracy
    
    print(paste0("Run ", i, " time."))
  }
} 

## Run prediction models: "xgboost" --------------------------------

## prepare input data ------------------------------------------------------
inputDat <- protein_Dat %>% purrr::reduce(rbind) %>% 
  cbind(mebo_Dat %>% purrr::reduce(rbind)) %>%
  as.data.frame %>% rownames_to_column("name")

dat_temp <- metadata_healthy %>% left_join(inputDat)

trainSet <- dat_temp %>% filter(season == "2015") %>% 
  rename("reclassify" = "H1N1_reclassify", "ab_d0_log2" = "H1N1_d0_log2")

### input variable --------------------------------
proNames <- colnames(protein_Dat$iMED_2014)
meboNames <- colnames(mebo_Dat$ZirFlu_2019)

inputVars <- c("age", "sex", "ab_d0_log2", proNames, meboNames)

# prepare input and validation sets
train_inputSet <- trainSet[, inputVars] %>% as.matrix()
train_predOut <- trainSet[, c("reclassify")]

### run the models --------------------------------
load("res.predModel_xgboost.RData")
for (i in 1:nreps) {
  val.pred_temp <- predict(netFits[[i]], newdata = train_inputSet) %>%  
    round() %>% as.data.frame() %>% rename("reclassify_numb" = ".") %>%
    mutate(reclassify = ifelse(reclassify_numb == 0, "LL", "protectee"),
           reclassify = factor(reclassify, levels = c("LL", "protectee")))
  
  netFit_pred_trainSet$xgboost[[i]] <- confusionMatrix(val.pred_temp$reclassify, train_predOut) # Training accuracy
  
  print(paste0("Run ", i, " time."))
}

# model predictions in validation set ====================================
netFit_parameters_total <- netFit_pred_trainSet %>%
  lapply(function(x) x %>% 
           lapply(function(y) y$byClass %>% 
                    as.data.frame %>% rownames_to_column("parameter")) %>% 
           bind_rows(.id = "run") %>% rename("value" = ".")) %>%
  bind_rows(.id = "model")

netFit_pred_params <- netFit_parameters_total %>% 
  filter(parameter %in% c("Sensitivity", "Specificity")) %>%
  pivot_wider(names_from = parameter, values_from = value)

# plot
netFit_pred_params %>%
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) + 
  geom_jitter(aes(col = model), size = 3, alpha = 0.9) + 
  geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
  ggtitle("Prediction accuracy in train set") +
  theme_classic()