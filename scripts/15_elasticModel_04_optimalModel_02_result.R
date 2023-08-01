rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(glmnet)

# Prediction outcome - repeat 1 time, data - metric ====================================
#load("res.elasticModel_optimal_rep30_0.RData")
#load("res.elasticModel_optimal_rep30_zscore_0.RData")

## model performance  -------------------------------------------------
netFit_parameters <- netFit_pred %>%
  lapply(function(x) x[[1]]$byClass %>% 
           as.data.frame %>% rownames_to_column("parameter")) %>% 
  bind_rows(.id = "vali_set") %>% rename("value" = ".")

para_plotDat <- netFit_parameters %>% 
  filter(parameter %in% c("Sensitivity", "Specificity")) %>%
  pivot_wider(names_from = parameter, values_from = value)

para_plotDat %>%
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) + 
  geom_point(color = "black", size = 6, alpha = 0.9) + 
  geom_point(aes(col = vali_set), size = 5, alpha = 0.9) + 
  geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
  ggtitle("Prediction accuracy in validation cohorts") +
  #ggtitle("Prediction in validation cohorts, normal metabolite data") +
  #ggtitle("Prediction in validation cohorts, zscore metabolite data") +
  theme_classic()

## prediction ------------------------------------------------------------------
pred_table <- netFit_pred %>%
  lapply(function(x) x[[1]]$table %>% as_tibble) %>% 
  bind_rows(.id = "vali_set") %>%
  mutate(pred_perform = ifelse(Prediction == Reference, TRUE, FALSE))

pred_table %>% 
  ggplot(aes(x = vali_set, y = n, fill = pred_perform)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_classic()

pred_table %>% 
  ggplot(aes(x = Reference, y = n, fill = pred_perform)) + 
  geom_bar(stat = "identity", position = "stack") + 
  facet_grid(~vali_set) +
  theme_bw() + theme(legend.position = "top")


#  Prediction outcome - repeat multiple times ======================================
netFits_total <- list()
netFit_pred_total <- list()

## load data ---------------------------------------------------------------
# load data - metric: Kappa, 19 times 
for (i in c(seq(0, 9), seq(11, 19))) {
  load(paste0("res.elasticModel_optimal_rep30_", i, ".RData"))
  netFits_total[[as.character(i)]] <- netFits
  netFit_pred_total[[as.character(i)]] <- netFit_pred
}

# load data - metric: Accuracy, 10 times
for (i in seq(0, 9)) {
  load(paste0("res.elasticModel_optimal_rep30_", i, "_accuracy.RData"))
  netFits_total[[as.character(i)]] <- netFits
  netFit_pred_total[[as.character(i)]] <- netFit_pred
}

# load data - metric: Kappa, 10 times with metabolite z-score input 
for (i in seq(0, 9)) {
  load(paste0("res.elasticModel_optimal_rep30_zscore_", i, ".RData"))
  netFits_total[[as.character(i)]] <- netFits
  netFit_pred_total[[as.character(i)]] <- netFit_pred
}

# models use z-score metabolite data
for (i in seq(0, 9)) {
  load(paste0("res.elasticModel_allinputs_Zscore_rep10_", i, ".RData")) # 10 time * rep10 = 100 times
  netFits_total[[as.character(i)]] <- netFits
  for (datSet in names(netFit_pred)) {
    netFit_pred_total[[datSet]][[as.character(i)]] <- netFit_pred[[datSet]]
  }
}
netFits_total <- netFits_total %>% flatten()
netFit_pred_total <- netFit_pred_total %>% lapply(function(x) x %>% flatten())

## plot data ---------------------------------------------------------------
netFit_parameters_total <- netFit_pred_total %>%
  lapply(function(x) x %>% 
           lapply(function(y) y[[1]]$byClass %>% 
                    as.data.frame %>% rownames_to_column("parameter")) %>% 
           bind_rows(.id = "vali_set") %>% rename("value" = ".")) %>%
  bind_rows(.id = "repeated_run")

para_plotDat_total <- netFit_parameters_total %>% 
  filter(parameter %in% c("Sensitivity", "Specificity")) %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(repeated_run = paste0(repeated_run, "_", vali_set))

netFit_parameters
para_plotDat_total %>%
  ggplot(aes(x = 1-Specificity, y = Sensitivity, )) + 
  #geom_point(color = "black", size = 6, alpha = 0.9) + 
  geom_jitter(aes(col = vali_set), size = 2, alpha = 0.9) + 
  geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
  ggtitle("Prediction accuracy in validation cohorts") +
  theme_classic()

para_plotDat_total <- netFit_parameters %>% 
  filter(parameter %in% c("Sensitivity", "Specificity")) %>%
  pivot_wider(names_from = parameter, values_from = value)

para_plotDat_total %>%
    ggplot(aes(x = 1-Specificity, y = Sensitivity)) + 
    #geom_point(color = "black", size = 6, alpha = 0.9) + 
    geom_jitter(aes(col = vali_set), size = 2, alpha = 0.9) + 
    geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
    ggtitle("Prediction accuracy in validation cohorts") +
    theme_classic()
  
## plot data for zscore metabolite ---------------------------------------------------------------
netFit_parameters_total <- netFit_pred_total %>%
  lapply(function(x) x %>% 
           lapply(function(y) y$byClass %>% 
                    as.data.frame %>% rownames_to_column("parameter")) %>% 
           bind_rows(.id = "repeated_run") %>% rename("value" = ".")) %>%
  bind_rows(.id = "vali_set")

para_plotDat_total <- netFit_parameters_total %>% 
  filter(parameter %in% c("Sensitivity", "Specificity")) %>%
  pivot_wider(names_from = parameter, values_from = value)

para_plotDat_total %>% 
  #filter(vali_set == "H1N1_2015_train") %>% 
  filter(repeated_run %in% H1N1_2015_goodPred$repeated_run) %>% 
  ggplot(aes(x = 1-Specificity, y = Sensitivity, )) + 
  #geom_point(color = "black", size = 2, alpha = 0.9) + 
  geom_jitter(aes(col = vali_set), size = 2, alpha = 0.9) + 
  #geom_jitter(aes(col = vali_set), size = 2, alpha = 0.9) + 
  geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
  #ggtitle("Prediction accuracy in train set - H1N1 2015") +
  ggtitle("Prediction accuracy in cohorts") +
  theme_classic()

H1N1_2015_goodPred <- para_plotDat_total %>% 
  filter(vali_set == "H1N1_2015_train") %>% filter(Sensitivity>0.75)

netFit_pred_total$H1N1_2015_train[[7]]
netFit_pred_total$H1N1_2014[[7]]
netFit_pred_total$H1N1_2019[[7]]
netFit_pred_total$H3N2_2015[[7]]
netFit_pred_total$B_2015[[7]]
netFit_pred_total$B_2014[[7]]
