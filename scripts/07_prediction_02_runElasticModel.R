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

get.elasticModel <- function(train_inputSet, train_predOut, 
                             validation_inputSets, validation_predOuts) {
  # Aim: run the model 10 times first and 50 time latter
  n <- 50
  netFit_model <- list()
  res_crosVal <- list()
  res_preVal <- list()
  coefs_Tab <- list()
  for (i in 1:n) {
    netFit <- train(x = train_inputSet,
                    y = train_predOut,
                    method = "glmnet", 
                    #metric = 'Accuracy', # for classification
                    metric = "Kappa", # similar to classification accuracy but it is useful to normalize the imbalance in classes
                    tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                                         .lambda = seq(0,1,by=0.01)),
                    trControl = trainControl(method="repeatedcv",
                                             number=5,
                                             repeats=20,
                                             summaryFunction=twoClassSummary, 
                                             classProbs=T,
                                             savePredictions = T))
    netFit_model[[i]] <- netFit # save model
    res_crosVal[[i]] <- get_best_result(netFit) # Test accuracy
    
    for (valiSet in names(validation_inputSets)) {
      val.pred_temp <- predict(netFit, newdata = validation_inputSets[[valiSet]], type = 'raw')
      res_preVal[[valiSet]][[i]] <- confusionMatrix(val.pred_temp, validation_predOuts[[valiSet]]) # Training accuracy
    }
    
    
    # Coefficients
    coefs <- coef(object = netFit$finalModel, s = netFit$bestTune$lambda)
    coefs_Tab[[i]] <- coefs[, 1] %>% as.data.frame()
    
    # running time
    print(paste0("Run ", i, " time."))
    
  }
  return(list("netFit_model" = netFit_model, 
              "res_crosVal" = res_crosVal, 
              "res_preVal" = res_preVal,
              "coefs_Tab" = coefs_Tab))
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

# Run elastic model ---------------------------------------------------------------
set.seed(123) #set the seed to make your partition reproducible


## input variable ----------------------------------------------------------------
load("processedDat/prediction_input.RData")

inputVariables <- list()
inputVariables$age_sex_proteins <- c("age", "sex", proNames)
inputVariables$age_sex_metabolites <- c("age", "sex", meboNames)
inputVariables$age_sex_abD0_proteins <- c("age", "sex", "ab_d0", proNames)
inputVariables$age_sex_abD0_metabolites <- c("age", "sex", "ab_d0", meboNames)
inputVariables$age_sex_abD0_proteins_metabolites <- c("age", "sex", "ab_d0", proNames, meboNames)
inputVariables$proteins_metabolites <- c(proNames, meboNames)

res.elasticModel <- list()
for (inputVal in names(inputVariables)) {
  # prepare input and validation sets
  train_inputSet <- trainSet[, inputVariables[[inputVal]]]
  train_predOut <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()
  
  validation_inputSets <- valiSets %>% 
    lapply(function(x) x[, inputVariables[[inputVal]]])

  validation_predOuts <- valiSets %>% 
    lapply(function(x) x[, c("reclassify")] %>% 
             as.vector() %>% unlist() %>% as.factor())
  
  # run the model
  res.elasticModel[[inputVal]] <- get.elasticModel(train_inputSet, train_predOut, 
                                                   validation_inputSets, validation_predOuts)
}

# save the result -----------------------
save(res.elasticModel, 
     file = "processedDat/prediction_resElasticModel.RData")
