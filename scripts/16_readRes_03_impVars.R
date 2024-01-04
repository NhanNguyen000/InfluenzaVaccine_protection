rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)

# load data : the weights of variables per model ===================================================
models <- c("glmnet", "nnet", "rf_tuneLength", "rf", "smv") # no: "knn", "xgboost"

weightedVars <- list()
for (nameModel in models) {
  load(paste0("res.predModel_", nameModel, ".RData"))
  
  weightedVars[[nameModel]] <- netFits %>% 
    lapply(function(x) varImp(x)) %>%
    lapply(function(x) x$importance) %>%
    imap(~.x %>% rename_with(function(x) paste(x, .y, sep = "_"))) %>%
    lapply(function(x) x %>% rownames_to_column("valName")) %>%
    purrr::reduce(full_join)
  
  rm(netFits, netFit_pred)
}

weightedVals_avg <- weightedVars %>% 
  lapply(function(x) x %>% 
           column_to_rownames("valName") %>% 
           rowMeans() %>% as.data.frame())

# Importance variables  ========================================================
# select important and consistent variable
impVars_names <- weightedVals_avg %>% 
  lapply(function(x) x %>% arrange(desc(.)) %>% top_n(., n = 100) %>% rownames()) # get the important variable

library(venn)
venn(impVars_names)

consistVars <- intersect(impVars_names$glmnet,
                         intersect(intersect(impVars_names$nnet, impVars_names$rf_tuneLength),
                                   intersect(impVars_names$rf, impVars_names$smv)))

# weight of variables
weightedVals_longDat <- weightedVars %>% 
  lapply(function(x) x %>% 
           pivot_longer(!valName, names_to = "run", values_to = "weight")) %>%
  bind_rows(.id = "model")

consistVars_longDat <- weightedVals_longDat %>% filter(valName %in% consistVars) 

# boxplot for variableds that are important and consistent across models
consistVars_longDat %>% 
  mutate(valName = factor(valName, levels = consistVars)) %>% 
  ggplot(aes(x = weight, y = valName)) +
  geom_boxplot() + geom_jitter(aes(color = model), size = 1.5) +
  theme_classic() + ggtitle("Important and consistent variable")
