# still could not run teh code, issue with the train()

rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
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
selected_crossVal <- trainControl(method="repeatedcv", number = 5, 
                                  repeats = 10)

selected_method <- "nb"
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
    val.pred_temp, train_predOut) # prediction accuracy
  
  # prediction for validation dataset
  for (valiSet in names(validation_inputSets)) {
    val.pred_temp <- predict(netFits[[i]], 
                             newdata = validation_inputSets[[valiSet]], 
                             type = 'raw')
    
    netFit_pred[[valiSet]][[i]] <- confusionMatrix(
      val.pred_temp, validation_predOuts[[valiSet]]) # prediction accuracy
  }
  
  print(paste0("Run ", i, " time."))
}

# save the result -----------------------
save(netFits, netFit_pred, file = "processedDat/prediction_res_nb.RData")


library(caret)
library(e1071)
library(klaR)

train_inputSet <- trainSet[, c("reclassify", inputVars)]
train_inputSet$reclassify <- as.factor(train_inputSet$reclassify)
a <- train(reclassify ~ ., data = train_inputSet ,
           method = "nb",
           trControl = trainControl(method = "cv", number = 10))

# test -------------
# Load the dataset (using iris as an example)
data(iris)

# Set up training control
train_control <- trainControl(method = "cv", number = 10)  # 10-fold cross-validation

# Train the model
set.seed(123)  # For reproducibility
nb_model <- train(Species ~ ., data = iris, method = "nb", trControl = train_control)

# Print the model
print(nb_model)

# Example validation set (replace this with your actual validation set)
validation_set <- iris[1:10, ]  # Using the first 10 rows as an example

# Predict using the trained model on the validation set
predictions <- predict(nb_model, newdata = validation_set)


# Print predictions
print(predictions)

# Compare predictions with actual values (if known)
actual_values <- validation_set$Species
comparison <- data.frame(Actual = actual_values, Predicted = predictions)
print(comparison)
