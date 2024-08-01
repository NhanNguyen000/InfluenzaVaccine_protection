rm(list = ls())

library(tidyverse)
library(magrittr)
library(xgboost)
library(caret)
library(pROC)

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
  lapply(function(x) x[, inputVars] %>% as.matrix())

validation_predOuts <- valiSets %>% 
  lapply(function(x)  x %>% 
           mutate(reclassify_numb = ifelse(reclassify == "LL", 0, 1)) %>%
           dplyr::select(reclassify, reclassify_numb))

## parameters for the algorithms--------------------------------------------------
predictOutcome <- list()

# train model
print("Train model")
predictModel <- xgboost(data = train_inputSet, 
                        label = train_predOut$reclassify_numb, 
                        max.depth = 2, eta = 1, nthread = 2, nrounds = 2, 
                        objective = "binary:logistic")

# prediction for train dataset
print("Validation")
val.pred_temp <- predict(predictModel, 
                         newdata = train_inputSet) %>%  
  round() %>% as.data.frame() %>% 
  rename("reclassify_numb" = ".") %>%
  mutate(reclassify = ifelse(reclassify_numb == 0, "LL", "protectee"),
         reclassify = factor(reclassify, levels = c("LL", "protectee")))

predictOutcome$trainSet <- confusionMatrix(
  val.pred_temp$reclassify, train_predOut$reclassify) # prediction accuracy

# prediction for validation dataset
for (valiSet in names(validation_inputSets)) {
  val.pred_temp <- predict(predictModel, 
                           newdata = validation_inputSets[[valiSet]]) %>%  
    round() %>% as.data.frame() %>% 
    rename("reclassify_numb" = ".") %>%
    mutate(reclassify = ifelse(reclassify_numb == 0, "LL", "protectee"),
           reclassify = factor(reclassify, levels = c("LL", "protectee")))
  
  predictOutcome[[valiSet]] <- confusionMatrix(
    val.pred_temp$reclassify, validation_predOuts[[valiSet]]$reclassify) # prediction accuracy
}

# ROC outcome
predict_ROC <- list()

val.pred_temp_v2 <- predict(predictModel, 
                            newdata = train_inputSet, 
                            type = 'prob')
predict_ROC$trainSet <- roc(train_predOut$reclassify, val.pred_temp_v2)

for (valiSet in names(validation_inputSets)) {
  val.pred_temp_v2 <- predict(predictModel, 
                              newdata = validation_inputSets[[valiSet]], 
                              type = 'prob')
  
  predict_ROC[[valiSet]] <- roc(
    validation_predOuts[[valiSet]]$reclassify, val.pred_temp_v2) # prediction accuracy
}


# save the result -----------------------
save(predictModel, predictOutcome, file = "processedDat/predictOutput_H1N1trainSet_xgboost.RData")
