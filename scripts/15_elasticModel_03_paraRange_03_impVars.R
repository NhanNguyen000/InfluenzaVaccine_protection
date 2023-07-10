rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(glmnet)

# load data =======================================================================
load("input.elasticModel.RData")

## input variable ----------------------------------------------------------------
inputVars <- c("age", "sex", "ab_T1", proNames, meboNames)

## prepare input and validation sets --------------------------------
train_inputSet <- trainSet[, inputVars]
train_predOut <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()

validation_inputSets <- valiSets %>% lapply(function(x) x[, inputVars])
validation_predOuts <- valiSets %>% 
  lapply(function(x) x[, c("reclassify")] %>% as.vector %>% unlist %>% as.factor)

# Run models ==========================================
set.seed(123) #set the seed to make your partition reproducible

#selected_metric = 'Accuracy' # for classification
selected_metric =  "Kappa" # similar to classification accuracy but it is useful to normalize the imbalance in classes

## run models to identify the important input variables ---------------------------------
selected_lamda_range <- seq(0, 0.5, by = 0.001)

defaut_alpha_range <- seq(0.1, 0.9, by = 0.1)
defaut_crossVal <- trainControl(method="repeatedcv", number = 5, repeats = 5)

nreps <- 100
netFits_temp <- list()
for (i in 1:nreps) {
  netFits_temp[[1]] <- train(x = train_inputSet, y = train_predOut,
                       method = "glmnet", metric = selected_metric,
                       tuneGrid=expand.grid(.alpha = defaut_alpha_range,
                                            .lambda = selected_lamda_range),
                       trControl = defaut_crossVal)
  
  print(paste0("Run ", i, " time."))
}

save(netFits_temp, file = "netFits_temp_varImp.RData")

# importantce variables
weightedVals_models <- netFits_temp %>% 
  lapply(function(x) varImp(x)) %>%
  lapply(function(x) x$importance) %>%
  imap(~.x %>% rename_with(function(x) paste(x, .y, sep = "_"))) %>%
  lapply(function(x) x %>% rownames_to_column("valName")) %>%
  purrr::reduce(full_join)

weightedVals_avg <- weightedVals_models %>% 
  column_to_rownames("valName") %>% 
  rowMeans() %>% as.data.frame()

varimp<- varImp(netFit_models_allDatTrain[[1]])
plot(varimp, main="Variable Importance")

varImp <- varImp(netFit_models_allDatTrain[[1]])
varImp$importance <- weightedVals_avg %>% filter(. > 0)
plot(varImp, main="Variable Importance")

impVars <- rownames(varImp$importance)

