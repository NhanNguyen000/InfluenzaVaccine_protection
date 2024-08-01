rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)

# load data =======================================================================
load("/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/netFit_models_allDatTrain.RData")

# importantce variables  =======================================================================

weightedVals_models <- netFit_models_allDatTrain %>% 
  lapply(function(x) varImp(x)) %>%
  lapply(function(x) x$importance) %>%
  imap(~.x %>% rename_with(function(x) paste(x, .y, sep = "_"))) %>%
  lapply(function(x) x %>% rownames_to_column("valName")) %>%
  purrr::reduce(full_join)

weightedVals_avg <- weightedVals_models %>% 
  column_to_rownames("valName") %>% 
  rowMeans() %>% as.data.frame()

varimp<- varImp(netFit_models_allDatTrain[[1]])
plot(varimp, main="Variable Importance")

varImp <- varImp(netFit_models_allDatTrain[[1]])
varImp$importance <- weightedVals_avg %>% filter(. > 0)
plot(varImp, main="Variable Importance")

impVars <- rownames(varImp$importance)
save(impVars, file = "impVars_elasticModels.RData")
