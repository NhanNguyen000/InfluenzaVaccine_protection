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

# save important variable information
save(varImp_list, file = "processedDat/predictOutcome_varImp.RData")

# important variables in 2 selected models ===================================================
varImp_selectedList <- list()
varImp_selectedList$H1N1_knn <- varImp_list$H1N1$knn
varImp_selectedList$B_knn <- varImp_list$B$knn

# make table inside the list
varImp_temp <- varImp_selectedList %>% 
  lapply(function(x) x$importance %>% 
           dplyr::rename_at(1, ~"weighted") %>%
           tibble::rownames_to_column("varName") %>% 
           dplyr::select(1, 2))

## venn diagram ------------------------------------------------------
library(venn)

varImp_topName <- varImp_temp %>% 
  lapply(function(x) x %>% 
           #top_n(30) %>% 
           dplyr::top_n(50) %>% 
           dplyr::select(1) %>% unlist())

varImp_topName_order <- list(
  "H1N1_knn" = varImp_topName$H1N1_knn,
  "B_knn" = varImp_topName$B_knn)

# save plot
png("output/prediction_varImp_top50_2selecteModels.png")
venn(varImp_topName_order, plotsize = 1.5, 
     ilcs = 3, sncs = 3, box = FALSE, 
     col = c("#9B2727", "#034E91"), zcolor = c("#9B2727", "#034E91"))
dev.off()

##  Variable weights in the models  ------------------------------------------------------
plotDat_varWeights <- varImp_temp%>% 
  lapply(function(x) x %>% tidyr::drop_na(weighted) %>% 
           dplyr::arrange(desc(weighted) ) %>% 
           tibble::rownames_to_column("order"))  %>%
  dplyr::bind_rows(.id = "model") %>%
  dplyr::mutate(order = as.numeric(order),
                model = factor(model, levels = c("H1N1_knn", "B_knn")))

varWeights_plot <- plotDat_varWeights %>% 
  ggplot(aes(x=order, y = weighted, group = model, color = model)) + 
  geom_line(size = 2, alpha = 0.8) + 
  geom_vline(xintercept = 50, linetype="dotted", alpha = 0.5, size=1.5) +
  xlab("List of variables") +  ylab("Variables' weight") +
  theme_classic() + 
  scale_x_continuous(breaks = c(0, 50, 200, 400, 600, 800)) +
  scale_color_manual(values=c( c("#9B2727", "#034E91"))) +
  theme(text = element_text(size = 28),
        legend.position = c(0.3, 0.8))

# save plot
varWeights_plot

png("output/prediction_varWeights_2selecteModels.png", width = 672)
varWeights_plot
dev.off()

# weights of the important variables in 2 selected models ===================================================
# selected variables
consistVars <- intersect(varImp_topName$H1N1_knn,varImp_topName$B_knn)
consistVars

selected_vars <- c(varImp_topName$H1N1_knn,varImp_topName$B_knn) %>% unique()

varImp_topName_table <- varImp_temp %>%
  dplyr::bind_rows(.id = "model") %>% 
  dplyr::group_by(model)

## plot important variables  ------------------------------------------------------
var_order <- varImp_temp$H1N1_knn %>% dplyr::arrange(desc(weighted))

plotDat <- varImp_topName_table  %>% 
  dplyr::filter(varName %in% selected_vars) %>% 
  dplyr::mutate(varName = factor(varName, levels = var_order$varName),
                color_code = ifelse(varName %in% consistVars,  "#034E91", "#898366"))

text_color <- (plotDat %>% filter(model == "H1N1_knn") %>% arrange(desc(weighted)))$color_code

varImp_plot_horizontal <- plotDat %>% ggplot(aes(x = varName))+
  geom_hline(yintercept= 0) +
  geom_hline(yintercept= c(-100, -50, 0, 50, 100), linetype="dashed", alpha = 0.5) +
  geom_col(aes( y = ifelse(model == "H1N1_knn", weighted, -weighted)), 
           fill = plotDat$color_code, width = 0.7, alpha = 0.7) +
  ylab("Variables' weight") + xlab("") +
  scale_y_continuous(labels=abs) +
  annotate("text", x=80, y=c(80, -80), label= c("H1N1_knn", "B_knn"), size = 20/.pt) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                   color = text_color), 
        text = element_text(size = 14),
        axis.title.y = element_text(vjust = -3)
  ) + coord_flip()


varImp_plot_vertical <- plotDat %>% ggplot(aes(y = varName))+
  geom_vline(xintercept= 0) +
  geom_vline(xintercept= c(-100, -50, 0, 50, 100), linetype="dashed", alpha = 0.5) +
  geom_col(aes( x = ifelse(model == "H1N1_knn", weighted, -weighted)), 
           fill = plotDat$color_code, width = 0.7, alpha = 0.7) +
  xlab("Variables' weight") + ylab("") +
  scale_x_continuous(labels=abs) +
  annotate("text", x=c(80, -80), y=80, label= c("H1N1_knn", "B_knn"), size = 20/.pt) +
  theme_classic() + 
  theme(axis.text.y = element_text(color = text_color), 
        text = element_text(size = 20)) 

# save plot
varImp_plot_horizontal

# png("output/prediction_varImp_top50_weight_2selecteModels.png", width = 1200)
# varImp_plot_horizontal 
# dev.off()


varImp_plot_vertical

png("output/prediction_varImp_top50_weight_2selecteModels.png", height = 1400, width = 576)
varImp_plot_vertical
dev.off()
