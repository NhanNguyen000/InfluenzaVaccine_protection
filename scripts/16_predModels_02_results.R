rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(glmnet)

# load data ====================================
load("res.predModel_rf.RData")


## model performance  -------------------------------------------------
netFit_parameters <- netFit_pred %>%
  lapply(function(x) x %>% 
           lapply(function(y) y$byClass %>% 
                    as.data.frame %>% rownames_to_column("parameter")) %>%
         bind_rows(.id = "run")) %>% 
  bind_rows(.id = "vali_set") %>% rename("value" = ".")

para_plotDat <- netFit_parameters %>% 
  filter(parameter %in% c("Sensitivity", "Specificity")) %>%
  pivot_wider(names_from = parameter, values_from = value)

para_plotDat %>%
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) + 
  geom_point(color = "black", size = 3, alpha = 0.9) + 
  geom_point(aes(col = vali_set), size = 2, alpha = 0.9) + 
  geom_abline(slope = 1, linetype = "dashed")+
  xlim(0, 1) + ylim(0, 1) + 
  ggtitle("Prediction accuracy in validation cohorts") +
  #ggtitle("Prediction in validation cohorts, normal metabolite data") +
  #ggtitle("Prediction in validation cohorts, zscore metabolite data") +
  theme_classic()

