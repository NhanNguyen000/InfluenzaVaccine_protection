rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}

get_best_result = function(caret_fit) {
  # Aim: get the bet predict models after training models with the data
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

get.elasticModel <- function(train_inputSet, train_predOut, 
                             validation_inputSets, validation_predOuts) {
  # Aim: run the model 10 times first and 50 time latter
  n <- 50
  netFit_model <- list()
  res_crosVal <- list()
  res_preVal <- list()
  coefs_Tab <- list()
  for (i in 1:n) {
    netFit <- train(x = train_inputSet,
                    y = train_predOut,
                    method = "glmnet", 
                    #metric = 'Accuracy', # for classification
                    metric = "Kappa",
                    tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                                         .lambda = seq(0,1,by=0.01)),
                    trControl = trainControl(method="repeatedcv",
                                             number=5,
                                             repeats=10))
    netFit_model[[i]] <- netFit # save model
    res_crosVal[[i]] <- get_best_result(netFit) # Test accuracy
    
    for (valiSet in names(validation_inputSets)) {
      val.pred_temp <- predict(netFit, newdata = validation_inputSets[[valiSet]], type = 'raw')
      res_preVal[[valiSet]][[i]] <- confusionMatrix(val.pred_temp, validation_predOuts[[valiSet]]) # Training accuracy
    }
    
    
    # Coefficients
    coefs <- coef(object = netFit$finalModel, s = netFit$bestTune$lambda)
    coefs_Tab[[i]] <- coefs[, 1] %>% as.data.frame()
    
    # running time
    print(paste0("Run ", i, " time."))
    
  }
  return(list("netFit_model" = netFit_model, 
              "res_crosVal" = res_crosVal, 
              "res_preVal" = res_preVal,
              "coefs_Tab" = coefs_Tab))
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


## prepare the protein and metabolite data -------------------------
# # iMED cohort 2014 
# metadat_iMED_2014 <- metadata_healthy %>% filter(season == "2014") 
# 
# inputDat_iMED_2014 <- protein_Dat$iMED_2014 %>% as.data.frame %>% rownames_to_column("name") %>%
#   full_join(mebo_Dat$iMED_2014 %>% as.data.frame %>% rownames_to_column("name")) %>%
#   column_to_rownames("name") %>% t() %>% as.data.frame %>% select(metadat_iMED_2014$name)
# # identical(colnames(inputDat_iMED_2014), metadat_iMED_2014$name) # TRUE, the same sample order
# 
# 
# # iMED cohort 2015
# metadat_iMED_2015 <- metadata_healthy %>% filter(season == "2015") 
# 
# inputDat_iMED_2015 <- protein_Dat$iMED_2015 %>% as.data.frame %>% rownames_to_column("name") %>%
#   full_join(mebo_Dat$iMED_2015 %>% as.data.frame %>% rownames_to_column("name")) %>%
#   column_to_rownames("name") %>% t() %>% as.data.frame %>% select(metadat_iMED_2015$name)
# # identical(colnames(inputDat_iMED_2015), metadat_iMED_2015$name) # TRUE, the same sample order
# 
# # ZirFlu cohort 2019
# metadat_ZirFlu_2019 <- metadata_healthy %>% filter(season == "2019") 
# 
# inputDat_ZirFlu_2019 <- protein_Dat$ZirFlu_2019 %>% as.data.frame %>% rownames_to_column("name") %>%
#   full_join(mebo_Dat$ZirFlu_2019 %>% as.data.frame %>% rownames_to_column("name")) %>%
#   column_to_rownames("name") %>% t() %>% as.data.frame %>% select(metadat_ZirFlu_2019$name)
# # identical(colnames(inputDat_ZirFlu_2019), metadat_ZirFlu_2019$name) # TRUE, the same sample order
# 
# # ZirFlu cohort 2020
# metadat_ZirFlu_2020 <- metadata_healthy %>% filter(season == "2020") 
# 
# inputDat_ZirFlu_2020 <- protein_Dat$ZirFlu_2020 %>% as.data.frame %>% rownames_to_column("name") %>%
#   full_join(mebo_Dat$ZirFlu_2020 %>% as.data.frame %>% rownames_to_column("name")) %>%
#   column_to_rownames("name") %>% t() %>% as.data.frame %>% select(metadat_ZirFlu_2020$name)
# # identical(colnames(inputDat_ZirFlu_2020), metadat_ZirFlu_2020$name) # TRUE, the same sample order


# impute protein data --------------------------------
# check NAs percentage
get.rmProteins <- function(dat, NA_cutoff) {
  # Description: identify variable have more missing value (NA) than the NA threshold (NA cutoff)
  
  # Arguments: 
  # dat - a numeric data table with sample in row (sample name is rowname), variable per column
  # selected_sample - name of the selected samples
  
  # Returns: variable name which have more NA than the NA_cutoff
  rmProteins <- which(colSums(is.na(dat))/nrow(dat) > NA_cutoff)
  return(rmProteins)
}

cutoff <- 0.2
rmProteins <- list()
for (season in names(protein_Dat)) {
  rmProteins[[season]] <- get.rmProteins(protein_Dat[[season]], NA_cutoff = cutoff)
} 
rmProteins # season 2019 and 2020, 14 proteins have a lot of NA --> can use 292/306 protein
rmProtein_names <- unique(rmProteins %>% lapply(function(x) names(x)) %>% unlist())

proteinDat_impute <- list()
for (season in names(protein_Dat)) {
  proteinDat_impute[[season]] <- protein_Dat[[season]] %>%
    #dplyr::select(-c(names(rmProteins[[season]]))) %>%
    select(-rmProtein_names) %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))
}

# testing code ------------------
inputDat <- proteinDat_impute %>% purrr::reduce(rbind) %>% 
  cbind(mebo_Dat %>% 
          lapply(function(x) x %>% select(-matches("_2"))) %>% 
          purrr::reduce(rbind)) %>%
  as.data.frame %>% rownames_to_column("name")

dat_temp <- metadata_healthy %>% left_join(inputDat)

dat_temp %>% count(season, H1N1_reclassify) # LL group size >= 3 in all 4 seasons -> can use all 4 seasons
dat_temp %>% count(season, H3N2_reclassify) # LL group size >= 3 in season 2015-> can use season 2015
dat_temp %>% count(season, B_reclassify) # LL group size >= 3 in all 2 seasons 2014 and 2015 -> can use all 2 seasons
dat_temp %>% count(season, Bvictoria_reclassify) # LL group size < 3 in all 2 seasons 2019 and 2020 --> can not use 
dat_temp %>% count(season, Byamagata_reclassify) # LL group size >= 3 in seasons 2020 --> can use season 2020


# Run elastic model --------------------------------
set.seed(123) #set the seed to make your partition reproducible

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


# input variable --------------------------------
proNames <- colnames(proteinDat_impute$iMED_2014)
meboNames <- colnames(mebo_Dat$ZirFlu_2019)

inputVariables <- list()
inputVariables$age_sex_proteins <- c("age", "sex", proNames)
inputVariables$age_sex_abT1_proteins <- c("age", "sex", "ab_T1", proNames)
inputVariables$age_sex_abT1_metabolites <- c("age", "sex", "ab_T1", meboNames)
inputVariables$age_sex_abT1_proteins_metabolites <- c("age", "sex", "ab_T1", proNames, meboNames)
inputVariables$age_sex_proteins_metabolites <- c("age", "sex", proNames, meboNames)


res.elasticModel <- list()
for (inputVal in names(inputVariables)) {
  # prepare input and validation sets
  train_inputSet <- trainSet[, inputVariables[[inputVal]]]
  train_predOut <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()
  
  validation_inputSets <- valiSets %>% 
    lapply(function(x) x[, inputVariables[[inputVal]]])

  validation_predOuts <- valiSets %>% 
    lapply(function(x) x[, c("reclassify")] %>% 
             as.vector() %>% unlist() %>% as.factor())
  
  # run the model
  res.elasticModel[[inputVal]] <- get.elasticModel(train_inputSet, train_predOut, 
                                                   validation_inputSets, validation_predOuts)
}

# save the result -----------------------
save(res.elasticModel, 
     file = "20230522_res.elasticModel.RData")
