rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(glmnet)
library(pROC)

# model predictions in validation set ====================================
models <- c("glmnet", "knn", "nnet", "rf_tuneGrid", "rf_tuneLength", "smv", "xgboost")
strains <- c("H1N1", "H3N2", "B")

outcome_prediction <- list()
outcome_ROC <- list()
for (strain in strains) {
  for (nameModel in models) {
    load(paste0("processedDat/predictOutput_", strain, "trainSet_", nameModel, ".RData"))
    
    outcome_prediction[[strain]][[nameModel]] <- predictOutcome %>% 
      lapply(function(x) x$byClass%>% 
               as.data.frame %>% tibble::rownames_to_column("parameter")) %>% 
      dplyr::bind_rows(.id = "vali_set") %>% dplyr::rename("value" = ".") 
    
    outcome_ROC[[strain]][[nameModel]] <- predict_ROC
  }
}

outcomeROC <- outcome_ROC %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% 
                    lapply(function(z) round(auc(z), 2)))) # AUC = Area under the curve

# plot ROC on traning set ------------------------------
roc_list <-  outcome_ROC %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y$trainSet)) 

roc_plotList <- list()

for (strain in strains) {
  roc_plotList[[strain]] <- roc_list[[strain]] %>% 
    ggroc(legacy.axes = TRUE, size = 1, alpha = 0.7) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          text = element_text(size = 16)) +  
    ggtitle(paste0(strain, " - ROC curves of different algorithms"))
  
}

# save plot
roc_plotList$H1N1

png("output/roc_H1N1in2015_trainingSet.png", width = 576)
roc_plotList$H1N1
dev.off()

roc_plotList$H3N2

png("output/roc_H3N2in2015_trainingSet.png", width = 576)
roc_plotList$H3N2
dev.off()

roc_plotList$B

png("output/roc_Bin2015_trainingSet.png", width = 576)
roc_plotList$B
dev.off()


# ROC boxplot for validation set-----------------------------------------------
outcomeROC_valiSet <- outcomeROC %>% 
  lapply(function(x) x%>%  dplyr::bind_rows(.id = "model")) %>%
  dplyr::bind_rows(.id = "trainSet_strain") %>% 
  dplyr::select(-trainSet) %>%
  tidyr::pivot_longer(cols = -c("model", "trainSet_strain"), 
               names_to = "vali_set", values_to = "AUC") %>% 
  tidyr::drop_na(AUC)

## box plot --------------------------------------------------------------------
plotList <- list()
plotList_v2 <- list()

for (strain in strains) {
  plotList[[strain]] <- outcomeROC_valiSet %>% 
    dplyr::filter(trainSet_strain == strain) %>%
    ggplot(aes(x = model, y = AUC)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 3, width = 0.2) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size = 16)) +  
    ggtitle(paste0(strain, " - AUC in validation sets"))
  
  plotList_v2[[strain]] <- outcomeROC_valiSet %>% 
    dplyr::filter(trainSet_strain == strain) %>%
    ggplot(aes(x = model, y = AUC)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(col = vali_set), size = 3, width = 0.2) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size = 16)) +  
    ggtitle(paste0(strain, " - AUC in validation sets"))
}

## box plot, set up the colors ------------------------------------
# blue scheme for same strain different year, gray for differnt strains
plotList_v3 <- list()

for (strain in strains) {
  # prepare data
  dat_temp <- outcomeROC_valiSet %>% 
    dplyr::filter(trainSet_strain == strain)  
  
  # sort the model based on their performance
  order_models <- dat_temp %>% 
    dplyr::group_by(model) %>% 
    dplyr::summarise(mean_AUC = mean(AUC)) %>% dplyr::arrange(desc(mean_AUC))
  
  # set up the colors: blue scheme for same strain different year, gray for differnt strains
  color_sameStrain <- dat_temp$vali_set[grep(strain, dat_temp$vali_set)] %>% unique()
  blue_scheme <- colorRampPalette(c("dodgerblue", "navyblue"))(length(color_sameStrain))
  
  color_diffStrain <- dat_temp$vali_set[-grep(strain, dat_temp$vali_set)] %>% unique()
  gray_scheme <- colorRampPalette(c("gray50", "gray75"))(length(color_diffStrain))
  
  plotList_v3[[strain]] <- dat_temp %>%
    dplyr::mutate(model = factor(model, levels = order_models$model),
           vali_set = factor(vali_set, 
                             levels = c(color_sameStrain, " ", color_diffStrain))) %>%
    ggplot(aes(x = model, y = AUC)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(col = vali_set), size = 5, width = 0.2, alpha = 0.7) +  
    ggtitle(paste0(strain, " - AUC in validation sets")) +
    scale_colour_manual(values = c(blue_scheme, "white", gray_scheme), drop = FALSE) +
    guides(color=guide_legend(ncol=2)) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 30, hjust=1),
          text = element_text(size = 20),
          legend.position = c(1, 1.1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6), 
          legend.title = element_blank(),
          legend.text=element_text(size=16))
}


# save plot
png("output/auc_H1N1in2015_valiSet.png", width = 672)
plotList_v3$H1N1
dev.off()

png("output/auc_H3N2in2015_valiSet.png", width = 720)
cowplot::plot_grid(plotList$H3N2,
                   plotList_v2$H3N2)
dev.off()

png("output/auc_Bin2015_valiSet.png", width = 672)
plotList_v3$B
dev.off()

# scatter plot: sensitivity vs. specificitys -------------------------------------------
outcome_Tab <- outcome_prediction %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% 
                    filter(parameter %in% c("Sensitivity", "Specificity")) %>%
                    pivot_wider(names_from = parameter, values_from = value)) %>%
           bind_rows(.id = "model")) %>% 
  bind_rows(.id = "strain_forTraining")

# plot function
plot_model <- function(dat, name) {
  dat %>% 
    ggplot(aes(x = 1-Specificity, y = Sensitivity)) + 
    geom_jitter(aes(shape = model, col = vali_set), 
                size = 4, alpha = 0.8,
                position = position_jitter(width = 0.01, height = 0.01, seed = 123)) + 
    geom_abline(slope = 1, linetype = "dashed")+
    xlim(0, 1) + ylim(0, 1) + 
    theme_classic() + 
    theme(text = element_text(size = 18)) +  
    ggtitle(paste0(name, " - Prediction accuracy"))
}

# plot
plotList <- list()

plotList$H1N1 <- outcome_Tab  %>% 
  filter(strain_forTraining == "H1N1") %>%
  filter(model %in% c("rf_tuneLength", "knn", "nnet")) %>%
  plot_model(., name = "H1N1")

plotList$H3N2 <- outcome_Tab  %>% 
  filter(strain_forTraining == "H3N2") %>%
  filter(model %in% c("xgboost", "smv")) %>%
  plot_model(., name = "H3N2")

plotList$B <- outcome_Tab  %>% 
  filter(strain_forTraining == "B") %>%
  filter(model %in% c("knn", "smv")) %>%
  plot_model(., name = "B")

# save plot
plotList$H1N1 

png("output/sensitivity_vsSpecificity_H1N1in2015_trainingSet.png", width = 576)
plotList$H1N1
dev.off()

plotList$H3N2

png("output/sensitivity_vsSpecificity_H3N2in2015_trainingSet.png", width = 576)
plotList$H3N2
dev.off()

plotList$B
png("output/sensitivity_vsSpecificity_HBin2015_trainingSet.png", width = 576)
plotList$B
dev.off()
