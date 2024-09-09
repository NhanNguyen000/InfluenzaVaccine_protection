rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(xgboost)

# load data : the weights of variables per model ==============================================
models <- c("glmnet", "knn", "nnet", "rf_tuneGrid", "rf_tuneLength", "smv", "xgboost") # 
strains <- c("H1N1", "H3N2", "B")

varImp_list <- list()
for (strain in strains) {
  for (nameModel in models) {
    load(paste0("processedDat/predictOutput_", strain, "trainSet_", nameModel, ".RData"))
    
    if (nameModel == "xgboost") {
      varImp_list[[strain]][[nameModel]] <-  xgb.importance(
        feature_names = predictModel$feature_names, model = predictModel)
    } else {
      varImp_list[[strain]][[nameModel]] <- varImp(predictModel)
    }
  }
}

# important variables in 4 selected models ===================================================
varImp_selectedList <- list()

varImp_selectedList$H1N1_rf_tuneLength <- varImp_list$H1N1$rf_tuneLength
varImp_selectedList$H1N1_knn <- varImp_list$H1N1$knn
#varImp_selectedList$H1N1_nnet <- varImp_list$H1N1$nnet

varImp_selectedList$B_knn <- varImp_list$B$knn
varImp_selectedList$B_smv <- varImp_list$B$smv

# make table inside the list
varImp_temp <- varImp_selectedList %>% 
  lapply(function(x) x$importance %>% 
           dplyr::rename_at(1, ~"weighted") %>%
           tibble::rownames_to_column("varName") %>% 
           dplyr::select(1, 2))

# venn diagram
library(venn)
varImp_topName <- varImp_temp %>% 
  lapply(function(x) x %>% 
           #top_n(50) %>% 
           dplyr::top_n(50) %>% 
           dplyr::select(1) %>% unlist())

varImp_topName_order <- list(
  "B_knn" = varImp_topName$B_knn, 
  "B_smv" = varImp_topName$B_smv,
  "H1N1_rf_tuneLength" = varImp_topName$H1N1_rf_tuneLength,   
  "H1N1_knn" = varImp_topName$H1N1_knn)

# save plot
png("output/prediction_varImp_top50.png")
venn(varImp_topName_order, plotsize = 1.5, 
     ilcs = 1.5, sncs = 1.5, box = FALSE)
dev.off()

# selecte variables
consistVars <- intersect(intersect(varImp_topName$H1N1_rf_tuneLength, varImp_topName$H1N1_knn),
                         intersect(varImp_topName$B_knn, varImp_topName$B_smv))
consistVars

# boxplot for variableds that are important and consistent across models
varImp_topName_table <- varImp_temp %>%
  bind_rows(.id = "model") %>% 
  group_by(model) %>%  
 # top_n(30) %>%
  top_n(50)

vars_count <- varImp_topName_table %>% group_by(varName) %>% count(varName)
selected_vars <- (vars_count %>% filter(n >= 3))$varName

plotDat <- varImp_topName_table  %>% 
  filter(varName %in% selected_vars) 

plotDat %>% 
  ggplot(aes(x = weighted, y = varName)) +
  geom_boxplot() + geom_jitter(aes(color = model), size = 1.5) +
  theme_classic() + ggtitle("Important and consistent variable")

# with order
var_medians <- plotDat %>% 
  group_by(varName) %>% 
  summarize(median_weighted = median(weighted)) %>% 
  arrange(desc(median_weighted))

varImp_plot <- plotDat %>% 
  mutate(varName = factor(varName, levels = var_medians$varName)) %>%
  ggplot(aes(x = weighted, y = varName)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = model), size = 4, alpha = 0.7) +
  theme_classic() +
  theme(text = element_text(size = 18)) +
  ggtitle("Important and consistent variable")


# save plot
varImp_plot 

png("output/prediction_varImp_top50_weight.png", width = 672)
varImp_plot 
dev.off()

# previous code ----------------------------------------------------------------------------------------------------------------------------------------------------------------


models <- c("glmnet", "nnet", "rf_tuneLength", "rf", "smv") # no: "knn", "xgboost"
models <- c("glmnet", "smv")

weightedVars <- list()
for (nameModel in models) {
  load(paste0("processedDat/prediction_res_", nameModel, ".RData"))
  
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

consistVars <- intersect(impVars_names$glmnet, impVars_names$smv)

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
