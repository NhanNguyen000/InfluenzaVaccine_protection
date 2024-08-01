rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(pROC)

# Run models (run in Slum/server) ==========================================
set.seed(123) #set the seed to make your partition reproducible

## input variable ----------------------------------------------------------------
load("processedDat/predictInput_H1N1trainSet.RData")
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
selected_crossVal <- trainControl(method="repeatedcv", number = 10, 
                                  repeats = 10)

selected_method <- "nnet"
predictOutcome <- list()

# train model
print("Train model")
predictModel <- train(x = train_inputSet, 
                      y = train_predOut,
                      method = selected_method, 
                      metric = selected_metric,
                      tuneLength =  30,
                      trControl = selected_crossVal)

# prediction for train dataset
print("Validation")
val.pred_temp <- predict(predictModel, 
                         newdata = train_inputSet, 
                         type = 'raw')

predictOutcome$trainSet <- confusionMatrix(
  val.pred_temp, train_predOut %>% as.factor()) # prediction accuracy

# prediction for validation dataset
for (valiSet in names(validation_inputSets)) {
  val.pred_temp <- predict(predictModel, 
                           newdata = validation_inputSets[[valiSet]], 
                           type = 'raw')
  
  predictOutcome[[valiSet]] <- confusionMatrix(
    val.pred_temp, validation_predOuts[[valiSet]]) # prediction accuracy
}

# ROC outcome
predict_ROC <- list()

val.pred_temp_v2 <- predict(predictModel, 
                            newdata = train_inputSet, 
                            type = 'prob')
predict_ROC$trainSet <- roc(train_predOut %>% as.factor(), val.pred_temp_v2$LL)

for (valiSet in names(validation_inputSets)) {
  val.pred_temp_v2 <- predict(predictModel, 
                              newdata = validation_inputSets[[valiSet]], 
                              type = 'prob')
  
  predict_ROC[[valiSet]] <- roc(
    validation_predOuts[[valiSet]], val.pred_temp_v2$LL) # prediction accuracy
}


# save the result -----------------------
save(predictModel, predictOutcome, file = "processedDat/predictOutput_H1N1trainSet_nnet.RData")
