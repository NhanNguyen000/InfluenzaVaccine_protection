rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(glmnet)

# load data =======================================================================
load("input.elasticModel.RData")
load("impVars_top.RData")

## input variable ----------------------------------------------------------------
#inputVars <- c("age", "sex", "ab_T1", proNames, meboNames)
inputVars <- c("age", "sex", "ab_T1", impVars_top)

## prepare input and validation sets---------------------------------------------------------------
train_inputSet <- trainSet[, inputVars]
train_predOut <- trainSet[, c("reclassify")] %>% 
  as.vector() %>% unlist()

validation_inputSets <- valiSets %>% lapply(function(x) x[, inputVars])

validation_predOuts <- valiSets %>% 
  lapply(function(x) x[, c("reclassify")] %>% as.vector() %>% unlist() %>% as.factor())
## prepare input and validation sets --------------------------------
train_inputSet <- trainSet[, inputVars]
train_predOut <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()

validation_inputSets <- valiSets %>% lapply(function(x) x[, inputVars])
validation_predOuts <- valiSets %>% 
  lapply(function(x) x[, c("reclassify")] %>% as.vector %>% unlist %>% as.factor)

# Run models (run in Slum/server) ==========================================
set.seed(123) #set the seed to make your partition reproducible

#selected_metric = 'Accuracy' # for classification
selected_metric =  "Kappa" # similar to classification accuracy but it is useful to normalize the imbalance in classes

## run models to identify the important input variables ---------------------------------
selected_lamda_range <- seq(0, 0.5, by = 0.001)
selected_alpha_range <- seq(0, 1, by = 0.001)
selected_crossVal <- trainControl(method="repeatedcv", number = 5, 
                                  repeats = 30)
#                                  repeats = 100)

nreps <- 1
netFits <- list()
netFit_pred <- list()
for (i in 1:nreps) {
  netFits[[i]] <- train(x = train_inputSet, y = train_predOut,
                        method = "glmnet", metric = selected_metric,
                        tuneGrid=expand.grid(.alpha = selected_alpha_range,
                                             .lambda = selected_lamda_range),
                        trControl = selected_crossVal)
  
  for (valiSet in names(validation_inputSets)) {
    val.pred_temp <- predict(netFits[[i]], newdata = validation_inputSets[[valiSet]], type = 'raw')
    netFit_pred[[valiSet]][[i]] <- confusionMatrix(val.pred_temp, validation_predOuts[[valiSet]]) # Training accuracy
  }
  
  print(paste0("Run ", i, " time."))
}

# save the result -----------------------
#save(netFits, netFit_pred, file = "res.elasticModel_optimal.RData")
save(netFits, netFit_pred, file = "res.elasticModel_optimal_rep30_10.RData")