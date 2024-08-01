rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}

get_best_result = function(caret_fit) {
  # Aim: get the bet predict models after training models with the data
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

# load data =======================================================================
load("processedDat/cohorts_dat.RData")

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
selected_crossVal <- trainControl(method="repeatedcv", 
                                  number = 5, 
                                  repeats = 30,
                                  summaryFunction=twoClassSummary, 
                                  classProbs=T,
                                  savePredictions = T)

## run prediction models -------------------------------------------------
method_list <- c("rf", "svmRadial", "nb", "knn", "nnet")
nreps <- 1 #10

print("Start run models.")

for (selected_method in method_list) {
  netFits <- list()
  netFit_pred <- list()
  netFit_crosVar <- list()
  coefs_Tab <- list()
  
  for (i in 1:nreps) {
    netFits[[selected_method]][[i]] <- train(x = train_inputSet, 
                                             y = train_predOut,
                                             method = selected_method, 
                                             metric = selected_metric,
                                             tuneLength =  30,
                                             trControl = selected_crossVal)
    
    netFit_model[[i]] <- netFit # save model
    netFit_crosVar[[i]] <- get_best_result(netFit) # Test accuracy
    
    for (valiSet in names(validation_inputSets)) {
      val.pred_temp <- predict(netFits[[selected_method]][[i]], 
                               newdata = validation_inputSets[[valiSet]], 
                               type = 'raw')
      
      netFit_pred[[selected_method]][[valiSet]][[i]] <- confusionMatrix(
        val.pred_temp, validation_predOuts[[valiSet]]) # Training accuracy
    }
    
    # Coefficients
    coefs <- coef(object = netFit$finalModel, s = netFit$bestTune$lambda)
    coefs_Tab[[i]] <- coefs[, 1] %>% as.data.frame()
    
    # running time
    print(paste0("Run ", i, " time."))
  }
  res <- list("netFit_model" = netFit_model, 
              "netFit_crosVar" = netFit_crosVar, 
              "netFit_pred" = netFit_pred,
              "coefs_Tab" = coefs_Tab)
  
  save(res, file = paste0("processedDat/prediction_res_", method, ".RData"))
}

# run algorithms which could not put in loop ---------------
# xgboost
library(xgboost)

netFits <- list()
netFit_pred <- list()

for (i in 1:nreps) {
  netFits[[i]] <- xgboost(data = train_inputSet, label = train_predOut$reclassify_numb, 
                          max.depth = 2, eta = 1, nthread = 2, nrounds = 2, 
                          objective = "binary:logistic")
  
  for (valiSet in names(validation_inputSets)) {
    val.pred_temp <- predict(netFits[[i]], newdata = validation_inputSets[[valiSet]]) %>%  
      round() %>% as.data.frame() %>% rename("reclassify_numb" = ".") %>%
      mutate(reclassify = ifelse(reclassify_numb == 0, "LL", "protectee"),
             reclassify = factor(reclassify, levels = c("LL", "protectee")))
    
    netFit_pred[[valiSet]][[i]] <- confusionMatrix(val.pred_temp$reclassify, 
                                                   validation_predOuts[[valiSet]]$reclassify) # Training accuracy
  }
  
  print(paste0("Run ", i, " time."))
}
save(res, file = "processedDat/prediction_res_xgboost.RData")
