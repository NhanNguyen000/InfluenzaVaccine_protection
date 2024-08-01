# NOTE: After run this code files, we found the optimal lamda value is in the range [0, 0.5]
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

## identify the range to search for the optimal lamda value (lamda can have range [0, infinity])-----------
defaut_alpha_range <- seq(0.1, 0.9, by = 0.1)
defaut_crossVal <- trainControl(method="repeatedcv", number = 5, repeats = 10)

# lamda range [0, 10]
netFit_temp <- train(x = train_inputSet, y = train_predOut,
                     method = "glmnet", metric = selected_metric,
                     tuneGrid=expand.grid(.alpha = defaut_alpha_range,
                       .lambda = seq(0, 10, by = 0.1)),
                     trControl = defaut_crossVal)

cowplot::plot_grid(plot(netFit_temp, metric = "Kappa"), 
                   plot(netFit_temp, metric = "Accuracy"), ncol = 1)


# lamda range [0, 5]
netFit_temp <- train(x = train_inputSet, y = train_predOut,
                     method = "glmnet", metric = selected_metric,
                     tuneGrid=expand.grid(.alpha = defaut_alpha_range,
                                          .lambda = seq(0, 5, by = 0.1)),
                     trControl = defaut_crossVal)

cowplot::plot_grid(plot(netFit_temp, metric = "Kappa"), 
                   plot(netFit_temp, metric = "Accuracy"), ncol = 1)

# lamda range [0, 1]
netFit_temp <- train(x = train_inputSet, y = train_predOut,
                     method = "glmnet", metric = selected_metric,
                     tuneGrid=expand.grid(.alpha = defaut_alpha_range,
                                          .lambda = seq(0, 1, by = 0.01)),
                     trControl = defaut_crossVal)

cowplot::plot_grid(plot(netFit_temp, metric = "Kappa"), 
                   plot(netFit_temp, metric = "Accuracy"), ncol = 1)

# lamda range [0, 0.5]
netFit_temp <- train(x = train_inputSet, y = train_predOut,
                     method = "glmnet", metric = selected_metric,
                     tuneGrid=expand.grid(.alpha = defaut_alpha_range,
                                          .lambda = seq(0, 0.5, by = 0.005)),
                     trControl = defaut_crossVal)

cowplot::plot_grid(plot(netFit_temp, metric = "Kappa"), 
                   plot(netFit_temp, metric = "Accuracy"), ncol = 1)

