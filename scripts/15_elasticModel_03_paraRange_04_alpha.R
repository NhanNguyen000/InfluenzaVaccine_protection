rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(glmnet)
library(latticeExtra) 
library(RColorBrewer)

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
selected_crossVal <- trainControl(method="repeatedcv", number = 5, repeats = 30)
defaut_alpha_range <- seq(0, 1, by = 0.01)

nreps <- 10
netFits <- list()
for (i in 1:nreps) {
  netFits[[i]] <- train(x = train_inputSet, y = train_predOut,
                        method = "glmnet", metric = selected_metric,
                        tuneGrid=expand.grid(.alpha = defaut_alpha_range,
                                             .lambda = selected_lamda_range),
                        trControl = selected_crossVal)
  
  print(paste0("Run ", i, " time."))
}

# save the result -----------------------
save(netFits, file = "netFits_temp_alpha.RData")

# optimal alpha range ==========================================
#load data
load("netFits_temp_alpha.RData")
netFits_1 <- netFits

load("netFits_temp_alpha_5.RData")
netFits_5 <- netFits

load("netFits_temp_alpha_10.RData")
netFits_10 <- netFits

netFits <- list()
netFits$run_1 <- netFits_1
netFits$run_5 <- netFits_5
netFits$run_10 <- netFits_10

netFits_temp <- netFits %>% purrr::flatten()

# check the Kappa value within lambda and alpha ranges
alpha_range <- netFits_temp %>% 
  lapply(function(x) x$results) %>%
  bind_rows(.id = "run")

levelplot(Kappa ~ alpha * lambda, alpha_range, 
          col.regions = colorRampPalette(brewer.pal(8, "Purples"))(100),
          main = list("Kappa value", side = 1, line = 0.5)) 

alpha_range %>% 
  ggplot(aes(x = alpha, y = Kappa, group = alpha)) + 
  geom_boxplot() + 
  theme_classic()

alpha_range %>% filter(alpha < 0.3) %>% 
  ggplot(aes(x = alpha, y = Kappa, group = alpha)) + 
  geom_boxplot() + 
  theme_classic()

# check the best tune (optimal lambda and alpha value) in 30 repeated models
bestTune_range <- netFits_temp %>% 
  lapply(function(x) x$bestTune) %>%
  bind_rows(.id = "run")

bestTune_range  %>% 
  ggplot(aes(x = alpha, y = lambda, group = alpha)) + 
  geom_boxplot() + geom_jitter(width = 0.001) +
  ggtitle("Best turn parameters") +
  theme_classic()
