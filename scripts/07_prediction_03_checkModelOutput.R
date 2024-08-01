rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(glmnet)

# model predictions in validation set ====================================
models <- c("knn", "nnet", "glmnet", "rf_tuneLength", "rf", "smv", "xgboost")
models <- c("knn", "glmnet", "smv", "xgboost")

netFit_pred_params <- list()
for (nameModel in models) {
  load(paste0("processedDat/prediction_res_", nameModel, ".RData"))
  
  netFit_parameters_total <- netFit_pred %>%
    lapply(function(x) x %>% 
             lapply(function(y) y$byClass %>% 
                      as.data.frame %>% rownames_to_column("parameter")) %>% 
             bind_rows(.id = "vali_set") %>% rename("value" = ".")) %>%
    bind_rows(.id = "repeated_run")
  
  netFit_pred_params[[nameModel]] <- netFit_parameters_total %>% 
    filter(parameter %in% c("Sensitivity", "Specificity")) %>%
    pivot_wider(names_from = parameter, values_from = value)
  
  rm(netFits, netFit_pred, netFit_parameters_total)
}

# all models
netFit_pred_params_total <- netFit_pred_params %>% bind_rows(.id = "model")

# prediction on train set
netFit_pred_params_total %>% filter(repeated_run == "trainSet") %>%
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) + 
  geom_jitter(aes(col = model), size = 3, alpha = 0.9) + 
  geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
  ggtitle("Prediction accuracy in train set") +
  theme_classic()

# prediction on validaton set
netFit_pred_params_total %>%
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) + 
  geom_jitter(aes(col = repeated_run, shape = model), size = 3, alpha = 0.9) + 
  geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
  ggtitle("Prediction accuracy in validation cohorts") +
  theme_classic()


# selected models: 
netFit_pred_params_total %>% 
  filter(model %in% c("knn", "nnet", "glmnet", "rf_tuneLength", "rf", "smv")) %>%
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) + 
  geom_jitter(aes(col = repeated_run, shape = model), size = 3, alpha = 0.9) + 
  geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
  ggtitle("Prediction accuracy in validation cohorts") +
  theme_classic()

# per model
nameModel <- "glmnet"
nameModel <- "nnet"
nameModel <- "rf"
nameModel <- "rf_tuneLength"
nameModel <- "smv"
nameModel <- "xgboost"

netFit_pred_params_total %>% 
  filter(model %in% nameModel) %>%
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) + 
  geom_jitter(aes(col = repeated_run), size = 2, alpha = 0.9) + 
  geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
  ggtitle(paste0("Prediction accuracy in validation cohorts - ", nameModel)) +
  theme_classic()

# other plot styles
netFit_pred_params$knn %>%
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) + 
  geom_point(color = "black", size = 3, alpha = 0.9) + 
  geom_point(aes(col = vali_set), size = 2, alpha = 0.9) + 
  geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
  ggtitle("Prediction accuracy in validation cohorts") +
  #ggtitle("Prediction in validation cohorts, normal metabolite data") +
  #ggtitle("Prediction in validation cohorts, zscore metabolite data") +
  theme_classic()
