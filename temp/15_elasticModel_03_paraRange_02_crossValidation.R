# NOTE: After run this code files, we found the optimal k-cros validatioin value is >= 10, but 5 is already quite good
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

## identify the value of alpha ---------------------------------
selected_lamda_range <- seq(0, 0.5, by = 0.001)

defaut_alpha_range <- seq(0.1, 0.9, by = 0.1)
defaut_crossVal <- trainControl(method="repeatedcv", number = 5, repeats = 10)

netFit_temp <- train(x = train_inputSet, y = train_predOut,
                     method = "glmnet", metric = selected_metric,
                     tuneGrid=expand.grid(.alpha = defaut_alpha_range,
                                          .lambda = selected_lamda_range),
                     trControl = defaut_crossVal)

# cowplot::plot_grid(plot(netFit_temp, metric = "Kappa"), 
#                    plot(netFit_temp, metric = "Accuracy"), ncol = 1)

netFit_temp$bestTune # alpha = 0.1

## use the value of alpha to identify the optimal value for the k-cross validation --------------------------------
reps <- 50
nfolds <- list()
nfolds$rep30 <- nfolds$rep20 <- nfolds$rep10 <- nfolds$rep5 <- nfolds$rep3 <- rep(NA, nreps)
for (i in 1:nreps) {
  nfolds$rep30[i] <- cv.glmnet(train_inputSet %>% as.matrix, train_predOut, 
                            alpha = 0.1, family = "binomial", type.measure = "auc",
                            nfolds = 10)$lambda.min %>% log()
  
  nfolds$rep20[i] <- cv.glmnet(train_inputSet %>% as.matrix, train_predOut, 
                            alpha = 0.1, family = "binomial", type.measure = "auc",
                            nfolds = 10)$lambda.min %>% log()
  
  nfolds$rep10[i] <- cv.glmnet(train_inputSet %>% as.matrix, train_predOut, 
                       alpha = 0.1, family = "binomial", type.measure = "auc",
                       nfolds = 10)$lambda.min %>% log()
  
  nfolds$rep5[i] <- cv.glmnet(train_inputSet %>% as.matrix, train_predOut, 
                       alpha = 0.1, family = "binomial", type.measure = "auc", 
                       nfolds = 5)$lambda.min %>% log()
  
  nfolds$rep3[i] <- cv.glmnet(train_inputSet %>% as.matrix, train_predOut, 
                       alpha = 0.1, family = "binomial", type.measure = "auc", 
                       nfolds = 3)$lambda.min %>% log()
}

boxplot(nfolds$rep3, nfolds$rep5, nfolds$rep10, nfolds$rep20, nfolds$rep30, 
        names = c("nfolds=3", "nfolds=5", "nfolds=10", "nfolds=20", "nfolds=30"), 
        ylab = "log(lambda.min)")
abline(h = log(0.5), col = "red", lty = 5, lwd = 2)
text(x = 1.5, y = -0.5, "log(lamda = 0.5)", col = "red")

