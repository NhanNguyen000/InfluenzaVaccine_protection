library(dplyr)
library(magrittr)
library(caret)

get_best_result = function(caret_fit) {
  # Aim: get the bet predict models after training models with the data
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

get.elasticModel <- function(train_inputSet, train_predOut, 
                             validation_inputSet, validation_predOut) {
  # Aim: run the model 50 times
  n <- 50
  res_crosVal <- list()
  res_preVal <- list()
  coefs_Tab <- list()
  for (i in 1:n) {
    netFit <- train(x = train_inputSet,
                    y = train_predOut,
                    method = "glmnet", metric = 'Accuracy', # for classification
                    tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                                         .lambda = seq(0,1,by=0.01)),
                    trControl = trainControl(method="repeatedcv",
                                             number=5,
                                             repeats=10))
    res_crosVal[[i]] <- get_best_result(netFit) # Test accuracy
    
    
    val.pred <- predict(netFit, newdata = validation_inputSet, type = 'raw')
    res_preVal[[i]] <- confusionMatrix(validation_predOut, val.pred) # Training accuracy
    
    # Coefficients
    coefs <- coef(object = netFit$finalModel, s = netFit$bestTune$lambda)
    if (length(unique(train_predOut)) == 2) {
      coefs_Tab[[i]] <- coefs[, 1] %>% as.data.frame()
    } else {
      coefs_Tab[[i]] <- coefs %>% 
        lapply(function(x) x[, 1] %>% as.data.frame) %>% 
        imap(~.x %>% rename_with(function(x) paste(.y))) %>% 
        lapply(function(x) x %>% rownames_to_column("Protein")) %>%
        purrr::reduce(full_join)
    }
    
    # runing time
    print(paste0("Run ", i, " time."))
    
  }
  return(list("res_crosVal" = res_crosVal, 
              "res_preVal" = res_preVal, 
              "coefs_Tab" = coefs_Tab))
}

# input data --------------------------
## metadata for all healthy subjects -------------------------
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(sex = ifelse(sex == "m", 0, 1))

## imputed protein -------------------------
inputDat_protein <- proteinDat_impute %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t() %>%
  as.data.frame() %>% rownames_to_column("OlinkID") %>% 
  left_join(cohorts_dat$proteinAnnot_all) %>% select(-OlinkID, -UniProt) %>%
  column_to_rownames("Assay") %>% drop_na()

## raw metabolite -------------------------
## iMED metabolites
iMED_metabolites <- iMED$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

iMED_metaboliteDat <- iMED$metabolite
if (identical(iMED_metabolites$ionIdx, as.numeric(colnames(iMED$metabolite))) == TRUE) {
  colnames(iMED_metaboliteDat) <- iMED_metabolites$Formula
}

## ZirFlu metabolites 
ZirFlu_metabolites <- ZirFlu$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

ZirFlu_metaboliteDat <- ZirFlu$metabolite_dat
if (identical(ZirFlu_metabolites$ionIdx, as.numeric(colnames(ZirFlu$metabolite_dat))) == TRUE) {
  colnames(ZirFlu_metaboliteDat) <- ZirFlu_metabolites$Formula
}

## overlap metabolites between iMED and ZirFlu 
overlapped_metabolites <- intersect(iMED_metabolites$Formula, ZirFlu_metabolites$Formula)

inputDat_metabolite <- ZirFlu_metaboliteDat[, overlapped_metabolites] %>% 
  as.data.frame %>% rownames_to_column("name") %>%
  full_join(iMED_metaboliteDat[, overlapped_metabolites] %>% 
              as.data.frame %>% rownames_to_column("name")) %>%
  filter(name %in% metadat_healthy$name) %>%
  column_to_rownames("name") %>% t() %>% as.data.frame

## raw metabolite -------------------------
inputDat <- inputDat_protein %>% rbind(inputDat_metabolite)

# Prepare model to run for 2 class (LL vs. Protector, ratio 1:1) --------------------------
dat_temp <- metadat_healthy %>% 
  full_join(inputDat %>% t() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  mutate(reclassify = ifelse(H1N1_reclassify == "LL", "LL", "Protector"))
dat_temp %>% count(reclassify) # LL group = 48 people
dat_temp %>% count(cohort, reclassify) # LL group = 48 people

## extract individual for LL and Protector (LH, HL, HH) groups for iMED cohort ---------------
dat_LL_iMED <- dat_temp %>% filter(cohort == "iMED", reclassify == "LL") # 32 people
dat_Protector_iMED <- dat_temp %>% filter(cohort == "iMED", reclassify != "LL") %>% 
  group_by(H1N1_reclassify) %>% sample_n(size = 11) # 32 people /3 groups = 10.6 ~11 people per group
dat_Protector_iMED %>% count(reclassify)

# merge data together as 1 input table
dat_iMED <- rbind(dat_LL_iMED, dat_Protector_iMED)
dat_iMED %>% count(reclassify) #equal LL (n = 32) and protector (n = 33)

## extract individual for LL and Protector (LH, HL, HH) groups for ZirFlu cohort ---------------
dat_ZirFlu <-  dat_temp %>% filter(cohort == "ZirFlu")
dat_ZirFlu %>% count(reclassify) #equal LL (n = 16) and protector (n = 50)

## Run elastic model --------------------------------
set.seed(123) #set the seed to make your partition reproducible

# split train-test set
trainSet <- dat_iMED
trainSet %>% count(reclassify)

validateSet <- dat_ZirFlu 
validateSet %>% count(reclassify)

# input variable
proteinNames <- rownames(inputDat_protein)
metaboliteNames <- rownames(inputDat_metabolite)

inputVariables <- list()
inputVariables$proteins <- proteinNames
inputVariables$age_proteins <- c("age", proteinNames)
inputVariables$age_sex_proteins <- c("age", "sex", proteinNames)
inputVariables$age_sex_abBaseline_proteins <- c("age", "sex", "H1N1_T1_combine", proteinNames)
inputVariables$age_sex_abBaseline_proteins_metabolites <- c("age", "sex", "H1N1_T1_combine", proteinNames, metaboliteNames)

res.elasticModel <- list()
for (inputVal in names(inputVariables)) {
  # prepare input and validation sets
  train_inputSet <- trainSet[, inputVariables[[inputVal]]]
  train_predOut <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()
  
  validation_inputSet <- validateSet[, inputVariables[[inputVal]]]
  validation_predOut <- validateSet[, c("reclassify")] %>% 
    as.vector() %>% unlist() %>% as.factor()
  
  # run the model
  res.elasticModel[[inputVal]] <- get.elasticModel(train_inputSet, train_predOut, 
                                                   validation_inputSet, validation_predOut)
}

res.elasticModel_2groups_equalRatio_splitCohort <- res.elasticModel
# save the result
save(res.elasticModel_2groups_equalRatio_splitCohort, file = "20230412_res.elasticModel_2groups_equalRatio_splitCohort.RData")


# Prepare model to run for 2 class (LL vs. Protector, ratio 1:2) --------------------------
dat_temp <- metadat_healthy %>% 
  full_join(inputDat %>% t() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  mutate(reclassify = ifelse(H1N1_reclassify == "LL", "LL", "Protector"))
dat_temp %>% count(reclassify) # LL group = 48 people
dat_temp %>% count(cohort, H1N1_reclassify) # LL group = 48 people

## extract individual for LL and Protector (LH, HL, HH) groups for iMED cohort ---------------
dat_LL_iMED <- dat_temp %>% filter(cohort == "iMED", reclassify == "LL") # 32 people
dat_Protector_iMED <- dat_temp %>% filter(cohort == "iMED", reclassify != "LL") %>% 
  group_by(H1N1_reclassify) %>% sample_n(size = 22) # (32*2) people /3 groups = 21.3 ~22 people per group
dat_Protector_iMED %>% count(reclassify)

# merge data together as 1 input table
dat_iMED <- rbind(dat_LL_iMED, dat_Protector_iMED)
dat_iMED %>% count(reclassify) #equal LL (n = 32) and protector (n = 66)

## extract individual for LL and Protector (LH, HL, HH) groups for ZirFlu cohort ---------------
dat_ZirFlu <-  dat_temp %>% filter(cohort == "ZirFlu")
dat_ZirFlu %>% count(reclassify) #equal LL (n = 16) and protector (n = 50)

## Run elastic model --------------------------------
set.seed(123) #set the seed to make your partition reproducible

# split train-test set
trainSet <- dat_iMED
trainSet %>% count(reclassify)

validateSet <- dat_ZirFlu 
validateSet %>% count(reclassify)

# input variable
proteinNames <- rownames(inputDat_protein)
metaboliteNames <- rownames(inputDat_metabolite)

inputVariables <- list()
inputVariables$proteins <- proteinNames
inputVariables$age_proteins <- c("age", proteinNames)
inputVariables$age_sex_proteins <- c("age", "sex", proteinNames)
inputVariables$age_sex_abBaseline_proteins <- c("age", "sex", "H1N1_T1_combine", proteinNames)
inputVariables$age_sex_abBaseline_proteins_metabolites <- c("age", "sex", "H1N1_T1_combine", proteinNames, metaboliteNames)

res.elasticModel <- list()
for (inputVal in names(inputVariables)) {
  # prepare input and validation sets
  train_inputSet <- trainSet[, inputVariables[[inputVal]]]
  train_predOut <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()
  
  validation_inputSet <- validateSet[, inputVariables[[inputVal]]]
  validation_predOut <- validateSet[, c("reclassify")] %>% 
    as.vector() %>% unlist() %>% as.factor()
  
  # run the model
  res.elasticModel[[inputVal]] <- get.elasticModel(train_inputSet, train_predOut, 
                                                   validation_inputSet, validation_predOut)
}
res.elasticModel_2groups_1to2Ratio_splitCohort <- res.elasticModel
# save the result
save(res.elasticModel_2groups_1to2Ratio_splitCohort, file = "20230412_res.elasticModel_2groups_1to2Ratio_splitCohort.RData")

