rm(list = ls())

library(tidyverse)
library(caret)
# load the data -------------------------
load("processedDat/prediction_resElasticModel.RData")

# load results from elastic model with combine cohorts
res_crosVal_Tab <- res.elasticModel %>%
  lapply(function(x) x$res_crosVal %>% purrr::reduce(rbind))
res_crosVal_Tab %>% lapply(function(x) range(x$Accuracy))

res_preVal_Tab <- res.elasticModel %>%
  lapply(function(x) x$res_preVal %>% 
           lapply(function(y) y %>% lapply(function(z) z$overall) %>% 
                    purrr::reduce(rbind) %>% as.data.frame()))
res_preVal_Tab %>% 
  lapply(function(x) x %>% lapply(function(y) range(y$Accuracy)))

res_preVal_Tab$age_sex_proteins %>% lapply(function(x) range(x$Accuracy))
res_preVal_Tab$age_sex_abD0_proteins_metabolites %>% lapply(function(x) range(x$Accuracy))
res_preVal_Tab$age_sex_proteins_metabolites %>% lapply(function(x) range(x$Accuracy))

res_preVal_Tab_v2 <- res.elasticModel %>%
  lapply(function(x) x$res_preVal %>% 
           lapply(function(y) y %>% lapply(function(z) z$byClass) %>% 
                    purrr::reduce(rbind) %>% as.data.frame()))
res_preVal_Tab_v2$age_sex_abD0_proteins_metabolites %>% lapply(function(x) range(x$Recall))
res_preVal_Tab_v2$age_sex_proteins_metabolites %>% lapply(function(x) range(x$Recall))


