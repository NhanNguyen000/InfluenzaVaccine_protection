rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(glmnet)

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}

get.rmProteins <- function(dat, NA_cutoff) {
  # Description: identify variable have more missing value (NA) than the NA threshold (NA cutoff)
  
  # Arguments: 
  # dat - a numeric data table with sample in row (sample name is rowname), variable per column
  # selected_sample - name of the selected samples
  
  # Returns: variable name which have more NA than the NA_cutoff
  rmProteins <- which(colSums(is.na(dat))/nrow(dat) > NA_cutoff)
  return(rmProteins)
}

# load data =======================================================================
load("/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~convert_protectees(.x)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "protectee"))) %>%
  mutate(sex = ifelse(sex == "m", 0, 1))

## impute protein data --------------------------------
# check NAs percentage
cutoff <- 0.2
rmProteins <- list()
for (season in names(protein_Dat)) {
  rmProteins[[season]] <- get.rmProteins(protein_Dat[[season]], NA_cutoff = cutoff)
} 
rmProteins # season 2019 and 2020, 14 proteins have a lot of NA --> can use 292/306 protein
rmProtein_names <- unique(rmProteins %>% lapply(function(x) names(x)) %>% unlist())

proteinDat_impute <- protein_Dat %>%
  lapply(function(x) x %>% select(-rmProtein_names) %>%
           mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x)))

## z-score metabolite data --------------------------------
mebo_Dat_zscore <- mebo_Dat %>% 
  lapply(function(x) x %>% mutate_all(~(scale(.) %>% as.vector())))

## prepare input data ------------------------------------------------------
inputDat <- proteinDat_impute %>% purrr::reduce(rbind) %>% 
  cbind(mebo_Dat_zscore %>% purrr::reduce(rbind)) %>%
  as.data.frame %>% rownames_to_column("name")

dat_temp <- metadata_healthy %>% left_join(inputDat)

## prepare train-test --------------------------------
trainSet <- dat_temp %>% filter(season == "2015") %>% 
  rename("reclassify" = "H1N1_reclassify", "ab_T1" = "H1N1_T1")
trainSet %>% count(reclassify)

valiSets <- list()
# for H1N1
for (year in c("2014", "2019", "2020")) {
  valiSets[[paste0("H1N1_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "H1N1_reclassify", "ab_T1" = "H1N1_T1")
}

# for H3N2
valiSets$H3N2_2015 <- dat_temp %>% 
  filter(season == "2015") %>% 
  rename("reclassify" = "H3N2_reclassify", "ab_T1" = "H3N2_T1")

# for B
for (year in c("2014", "2015")) {
  valiSets[[paste0("B_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "B_reclassify", "ab_T1" = "B_T1")
}

# for Byamagata
valiSets$Byamagata<- dat_temp %>% 
  filter(season == "2020") %>% 
  rename("reclassify" = "Byamagata_reclassify", "ab_T1" = "Byamagata_T1")

valiSets %>% lapply(function(x) x%>% count(reclassify))

## input variable ----------------------------------------------------------------
proNames <- colnames(proteinDat_impute$iMED_2014)
meboNames <- colnames(mebo_Dat$ZirFlu_2019)
inputVars <- c("age", "sex", "ab_T1", proNames, meboNames)

## prepare input and validation sets-----------------------------------------
train_inputSet <- trainSet[, inputVars]
train_predOut <- trainSet[, c("reclassify")] %>% 
  as.vector() %>% unlist()

validation_inputSets <- valiSets %>% lapply(function(x) x[, inputVars])

validation_predOuts <- valiSets %>% 
  lapply(function(x) x[, c("reclassify")] %>% as.vector() %>% unlist() %>% as.factor())

# Run models (run in Slum/server) ==========================================
set.seed(123) #set the seed to make your partition reproducible

#selected_metric = 'Accuracy' # for classification
selected_metric =  "Kappa" # similar to classification accuracy but it is useful to normalize the imbalance in classes

## run models to identify the important input variables ---------------------------------
lamda_range <- seq(0, 1.5, by = 0.01) # can reduce from range [0, 10] with step 0.1 to [0, 2.5] step 0.01, to [0, 1.5]
alpha_range <- seq(0, 1, by = 0.01)
crossVal <- trainControl(method="repeatedcv", number = 5, repeats = 10)

nreps <- 10
netFits <- list()
netFit_pred <- list()
for (i in 1:nreps) {
  netFits[[i]] <- train(x = train_inputSet, y = train_predOut,
                        method = "glmnet", metric = selected_metric,
                        tuneGrid=expand.grid(.alpha = alpha_range,
                                             .lambda = lamda_range),
                        trControl = crossVal)
  
  val.pred_temp <- predict(netFits[[i]], newdata = train_inputSet, type = 'raw')
  netFit_pred[["H1N1_2015_train"]][[i]] <- confusionMatrix(val.pred_temp, train_predOut %>% as.factor()) # Training accuracy
  
  for (valiSet in names(validation_inputSets)) {
    val.pred_temp <- predict(netFits[[i]], newdata = validation_inputSets[[valiSet]], type = 'raw')
    netFit_pred[[valiSet]][[i]] <- confusionMatrix(val.pred_temp, validation_predOuts[[valiSet]]) # Training accuracy
  }
  
  print(paste0("Run ", i, " time."))
}

# check lamda range
cowplot::plot_grid(plot(netFits[[1]], metric = "Kappa"), 
                   plot(netFits[[1]], metric = "Accuracy"), ncol = 1)
# save the result -----------------------
save(netFits, netFit_pred, file = "res.elasticModel_allinputs_Zscore_rep10_9.RData")
