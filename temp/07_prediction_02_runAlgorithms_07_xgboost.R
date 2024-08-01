rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(xgboost)

# Run models (run in Slum/server) ==========================================
set.seed(123) #set the seed to make your partition reproducible

## input variable ----------------------------------------------------------------
load("processedDat/predictInput_H1N1trainSet.RData")
inputVars <- c("age", "sex", "ab_d0", proNames, meboNames)

# prepare input and validation sets
train_inputSet <- trainSet[, inputVars] %>% as.matrix()

train_predOut <- trainSet %>% 
  mutate(reclassify_numb = ifelse(reclassify == "LL", 0, 1)) %>%
  dplyr::select(reclassify, reclassify_numb)


validation_inputSets <- valiSets %>% 
  lapply(function(x) x[, inputVars]) %>% as.matrix()

validation_predOuts <- valiSets %>% 
  lapply(function(x)  x %>% 
           mutate(reclassify_numb = ifelse(reclassify == "LL", 0, 1)) %>%
           dplyr::select(reclassify, reclassify_numb))


## parameters for the algorithms--------------------------------------------------
nreps <- 10

netFits <- list()
netFit_pred <- list()
for (i in 1:nreps) {
  # train model
  netFits[[i]] <- xgboost(data = train_inputSet, 
                          label = train_predOut$reclassify_numb, 
                          max.depth = 2, eta = 1, nthread = 2, nrounds = 2, 
                          objective = "binary:logistic")
  
  # prediction for train dataset
  val.pred_temp <- predict(netFits[[i]], 
                           newdata = train_inputSet) %>%  
    round() %>% as.data.frame() %>% 
    rename("reclassify_numb" = ".") %>%
    mutate(reclassify = ifelse(reclassify_numb == 0, "LL", "protectee"),
           reclassify = factor(reclassify, levels = c("LL", "protectee")))
  
  netFit_pred$trainSet[[i]] <- confusionMatrix(
    val.pred_temp$reclassify, train_predOut$reclassify) # prediction accuracy
  
  # prediction for validation dataset
  for (valiSet in names(validation_inputSets)) {
    val.pred_temp <- predict(netFits[[i]], 
                             newdata = validation_inputSets[[valiSet]]) %>%  
      round() %>% as.data.frame() %>% 
      rename("reclassify_numb" = ".") %>%
      mutate(reclassify = ifelse(reclassify_numb == 0, "LL", "protectee"),
             reclassify = factor(reclassify, levels = c("LL", "protectee")))
    
    netFit_pred[[valiSet]][[i]] <- confusionMatrix(
      val.pred_temp$reclassify, validation_predOuts[[valiSet]]$reclassify) # prediction accuracy
  }
  
  print(paste0("Run ", i, " time."))
}

# save the result -----------------------
save(netFits, netFit_pred, file = "processedDat/prediction_res_xgboost.RData")
