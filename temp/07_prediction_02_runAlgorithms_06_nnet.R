rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(pROC)

# Run models (run in Slum/server) ==========================================
set.seed(123) #set the seed to make your partition reproducible

## input variable ----------------------------------------------------------------
load("processedDat/prediction_input.RData")
inputVars <- c("age", "sex", "ab_d0", proNames, meboNames)

# prepare input and validation sets
train_inputSet <- trainSet[, inputVars]
train_predOut <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()

validation_inputSets <- valiSets %>% lapply(function(x) x[, inputVars])

validation_predOuts <- valiSets %>% 
  lapply(function(x) x[, c("reclassify")] %>% 
           as.vector() %>% unlist() %>% as.factor())


## parameters for the algorithms--------------------------------------------------
#selected_metric = 'Accuracy' # for classification
selected_metric =  "Kappa" # similar to classification accuracy but it is useful to normalize the imbalance in classes
selected_crossVal <- trainControl(method="repeatedcv", number = 5, 
                                  repeats = 10)

selected_method <- "nnet"
nreps <- 10
netFits <- list()
netFit_pred <- list()
for (i in 1:nreps) {
  # train model
  netFits[[i]] <- train(x = train_inputSet, 
                        y = train_predOut,
                        method = selected_method, 
                        metric = selected_metric,
                        tuneLength =  30,
                        trControl = selected_crossVal)
  
  # prediction for train dataset
  val.pred_temp <- predict(netFits[[i]], 
                           newdata = train_inputSet, 
                           type = 'raw')
  
  netFit_pred$trainSet[[i]] <- confusionMatrix(
    val.pred_temp, train_predOut %>% as.factor()) # prediction accuracy
  
  # prediction for validation dataset
  for (valiSet in names(validation_inputSets)) {
    val.pred_temp <- predict(netFits[[i]], 
                             newdata = validation_inputSets[[valiSet]], 
                             type = 'raw')
    
    netFit_pred[[valiSet]][[i]] <- confusionMatrix(
      val.pred_temp, validation_predOuts[[valiSet]])  # prediction accuracy
  }
  
  print(paste0("Run ", i, " time."))
}

# save the result -----------------------
save(netFits, netFit_pred, file = "processedDat/prediction_res_nnet.RData")
