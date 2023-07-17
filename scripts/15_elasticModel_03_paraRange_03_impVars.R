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

# Run models (run in Slum/server) ==========================================
set.seed(123) #set the seed to make your partition reproducible

#selected_metric = 'Accuracy' # for classification
selected_metric =  "Kappa" # similar to classification accuracy but it is useful to normalize the imbalance in classes

## run models to identify the important input variables ---------------------------------
selected_lamda_range <- seq(0, 0.5, by = 0.001)

defaut_alpha_range <- seq(0.1, 0.9, by = 0.1)
defaut_crossVal <- trainControl(method="repeatedcv", number = 5, repeats = 5)

nreps <- 200
netFits_temp <- list()
for (i in 1:nreps) {
  netFits_temp[[i]] <- train(x = train_inputSet, y = train_predOut,
                       method = "glmnet", metric = selected_metric,
                       tuneGrid=expand.grid(.alpha = defaut_alpha_range,
                                            .lambda = selected_lamda_range),
                       trControl = defaut_crossVal)
  
  print(paste0("Run ", i, " time."))
}

save(netFits_temp, file = "netFits_temp_varImp.RData")

# Importance variables ==========================================
load("netFits_temp_varImp.RData")

impVars <- netFits_temp %>% 
  lapply(function(x) varImp(x)) %>%
  lapply(function(x) x$importance) %>%
  imap(~.x %>% rename_with(function(x) paste(x, .y, sep = "_"))) %>%
  lapply(function(x) x %>% rownames_to_column("valName")) %>%
  purrr::reduce(full_join)

impVars_avg <- impVars %>% 
  column_to_rownames("valName") %>% 
  rowMeans() %>% as.data.frame()

impVars_names <- impVars_avg %>% filter(. > 0) %>% arrange(desc(.)) %>% rownames() # get the important variable
impVars_top <- impVars_avg %>% arrange(desc(.)) %>% top_n(., n = 50) %>% rownames()
#save(impVars_top, file = "impVars_top.RData")

impVars_longDat <- impVars %>% 
  pivot_longer(!valName, names_to = "model", values_to = "weight")

# variable important in 1 model
varImp <- varImp(netFits_temp[[1]])
varImp$importance <- impVars_avg %>% filter(. > 0)
plot(varImp, main="Variable Importance")

# average of variable important across models
plot(x = impVars_avg$., y = row.names(impVars_avg), 
     main="Variable Importance") # still need to fix the code

# boxplot for variable important across models
impVars_longDat %>% 
  filter(valName %in% impVars_top) %>%
  mutate(valName = factor(valName, levels = impVars_top)) %>% 
  ggplot(aes(x = weight, y = valName)) +
  geom_boxplot() + geom_jitter(size = 0.3) +
  theme_classic() + ggtitle("Top important variable")

# calculate the contribution per groups: protien, metabolites -----------
impVars_avg_v2 <- impVars_avg %>%
  rownames_to_column("valName") %>%
  mutate(type = ifelse(valName %in% proNames, "protein",
                       ifelse(valName %in% meboNames, "metabolite", valName))) %>%
  rename("avg_coef" = ".")

impVars_avg_v2 %>% count(type)

impVars_avg_v2 %>% 
  ggplot(aes(x = type, y = avg_coef)) +
  geom_boxplot() + #geom_jitter(size = 0.3) +
  theme_classic() + ggtitle("Important variable by types")

impVars_avg_v2 %>% 
  ggplot(aes(avg_coef)) +
  geom_histogram() + #geom_jitter(size = 0.3) +
  theme_classic() + ggtitle("Important variable by types")
