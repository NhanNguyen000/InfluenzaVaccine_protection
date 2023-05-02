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
# metadata for all healthy subjects
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(sex = ifelse(sex == "m", 0, 1))

# impute protein
inputDat <- proteinDat_impute %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t() %>%
  as.data.frame() %>% rownames_to_column("OlinkID") %>% 
  left_join(cohorts_dat$proteinAnnot_all) %>% select(-OlinkID, -UniProt) %>%
  column_to_rownames("Assay") %>% drop_na()

# Prepare model to run for 4 class --------------------------------------------
dat <- metadat_healthy %>% 
  full_join(inputDat %>% t() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  mutate(reclassify = H1N1_reclassify)
dat %>% count(reclassify)

## Run elastic model --------------------------------
set.seed(123) #set the seed to make your partition reproducible

# split train-test set
trainSet <- dat %>% group_by(reclassify) %>% slice_sample(prop = 0.70)
trainSet %>% count(reclassify)

validateSet <- dat %>% anti_join(trainSet)
validateSet %>% count(reclassify)

# input variable
proteinNames <- rownames(inputDat)

inputVariables <- list()
inputVariables$proteins <- proteinNames
inputVariables$age_proteins <- c("age", proteinNames)
inputVariables$age_sex_proteins <- c("age", "sex", proteinNames)
inputVariables$age_sex_abBaseline_proteins <- c("age", "sex", "H1N1_T1_combine", proteinNames)

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
res.elasticModel_4groups <- res.elasticModel
# save the result
save(res.elasticModel_4groups , file = "20230412_res.elasticModel_4groups.RData")

# Prepare model to run for 2 class (LL vs. Protector, ratio 1:1) --------------------------
dat_temp <- metadat_healthy %>% 
  full_join(inputDat %>% t() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  mutate(reclassify = ifelse(H1N1_reclassify == "LL", "LL", "Protector"))
dat_temp %>% count(reclassify) # LL group = 48 people

# extract individual for LL and Protector (LH, HL, HH) groups
dat_LL <- dat_temp %>% filter(reclassify == "LL") # 48 people
dat_Protector <- dat_temp %>% filter(reclassify != "LL") %>% 
  group_by(H1N1_reclassify) %>% sample_n(size = 16) # 48 people /3 groups = 16 people per group
dat_Protector %>% count(reclassify)

# merge data together as 1 input table
dat <- rbind(dat_LL, dat_Protector)
dat %>% count(reclassify) #equal LL (n = 48) and protector (n = 48)

## Run elastic model --------------------------------
set.seed(123) #set the seed to make your partition reproducible

# split train-test set
trainSet <- dat %>% group_by(reclassify) %>% slice_sample(prop = 0.70)
trainSet %>% count(reclassify)

validateSet <- dat %>% anti_join(trainSet)
validateSet %>% count(reclassify)

# input variable
proteinNames <- rownames(inputDat)

inputVariables <- list()
inputVariables$proteins <- proteinNames
inputVariables$age_proteins <- c("age", proteinNames)
inputVariables$age_sex_proteins <- c("age", "sex", proteinNames)
inputVariables$age_sex_abBaseline_proteins <- c("age", "sex", "H1N1_T1_combine", proteinNames)

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
res.elasticModel_2groups_equalRatio <- res.elasticModel
# save the result
save(res.elasticModel_2groups_equalRatio , file = "20230412_res.elasticModel_2groups_equalRatio.RData")


# Prepare model to run for 2 class (LL vs. Protector, ratio 1:2) --------------------------
dat_temp <- metadat_healthy %>% 
  full_join(inputDat %>% t() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  mutate(reclassify = ifelse(H1N1_reclassify == "LL", "LL", "Protector"))
dat_temp %>% count(reclassify) # LL group = 48 people

# extract individual for LL and Protector (LH, HL, HH) groups
dat_LL <- dat_temp %>% filter(reclassify == "LL") # 48 people
dat_Protector <- dat_temp %>% filter(reclassify != "LL") %>% 
  group_by(H1N1_reclassify) %>% sample_n(size = 32) # (48*2) people /3 groups = 32 people per group
dat_Protector %>% count(reclassify)

# merge data together as 1 input table
dat <- rbind(dat_LL, dat_Protector)
dat %>% count(reclassify) #equal LL (n = 48) and protector (n = 96)

## Run elastic model --------------------------------
set.seed(123) #set the seed to make your partition reproducible

# split train-test set
trainSet <- dat %>% group_by(reclassify) %>% slice_sample(prop = 0.70)
trainSet %>% count(reclassify)

validateSet <- dat %>% anti_join(trainSet)
validateSet %>% count(reclassify)

# input variable
proteinNames <- rownames(inputDat)

inputVariables <- list()
inputVariables$proteins <- proteinNames
inputVariables$age_proteins <- c("age", proteinNames)
inputVariables$age_sex_proteins <- c("age", "sex", proteinNames)
inputVariables$age_sex_abBaseline_proteins <- c("age", "sex", "H1N1_T1_combine", proteinNames)

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
res.elasticModel_2groups_1to2Ratio <- res.elasticModel
# save the result
save(res.elasticModel_2groups_1to2Ratio, file = "20230412_res.elasticModel_2groups_1to2Ratio.RData")

