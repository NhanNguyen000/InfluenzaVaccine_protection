# function & code ------------------------
library(dplyr)
library(magrittr)
library(caret)
library(tidyverse)


get_best_result = function(caret_fit) {
  # Aim: get the bet predict models after training models with the data
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

get.elasticModel <- function(train_inputSet, train_predOut, 
                             validation_inputSet_1, validation_predOut_1,
                             validation_inputSet_2, validation_predOut_2,
                             validation_inputSet_3, validation_predOut_3) {
  # Aim: run the model 10 times first and 50 time latter
  n <- 50
  netFit_model <- list()
  res_crosVal <- list()
  res_preVal_1 <- list()
  res_preVal_2 <- list()
  res_preVal_3 <- list()
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
    netFit_model[[i]] <- netFit # save model
    res_crosVal[[i]] <- get_best_result(netFit) # Test accuracy
    
    
    val.pred_1 <- predict(netFit, newdata = validation_inputSet_1, type = 'raw')
    res_preVal_1[[i]] <- confusionMatrix(val.pred_1, validation_predOut_1) # Training accuracy 1
    
    val.pred_2 <- predict(netFit, newdata = validation_inputSet_2, type = 'raw')
    res_preVal_2[[i]] <- confusionMatrix(val.pred_2, validation_predOut_2) # Training accuracy 2
    
    val.pred_3 <- predict(netFit, newdata = validation_inputSet_3, type = 'raw')
    res_preVal_3[[i]] <- confusionMatrix(val.pred_3, validation_predOut_3) # Training accuracy 3
    
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
  return(list("netFit_model" = netFit_model, 
              "res_crosVal" = res_crosVal, 
              "res_preVal_1" = res_preVal_1, 
              "res_preVal_2" = res_preVal_2, 
              "res_preVal_3" = res_preVal_3, 
              "coefs_Tab" = coefs_Tab))
}


# load data from iMED and ZirFlu cohorts -------------------------------
# load raw data
load("../ZirFlu_NhanNguyen/ZirFlu.RData")
load("../iMED_NhanNguyen/iMED.RData")

# load processed and cleanned data 
load("cohorts_dat.RData")
load("protein_normDat.RData")
load("proteinDat_impute.RData")
load("HAIreclassify.RData")
load("HAIreclassify2.RData")

# prepare input data --------------------------
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

## raw protein + metabolite -------------------------
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

# iMED pilot cohort data ---------------------
# iMED pilot cohort we have metabolites, Ab titers, proteome. 
# Metadata is here: /vol/projects/CIIM/Influenza/iMED/metadata . 
# Proteome: /vol/projects/CIIM/Influenza/iMED/proteomic/analysis/output/cohort1.RData , 
# metabolome: /vol/projects/CIIM/Influenza/iMED/metabolic/analysis/brendan/output/cohort_1/cohort1.Rdata

iMEDpilot_ab <- read.table("/vol/projects/CIIM/Influenza/iMED/metadata/hai_titers_pilot.tsv", header = TRUE)
iMEDpilot_metadat <- iMEDpilot_ab  %>% 
  mutate(ProbandID = as.numeric(substring(pID, 3, 5))) %>% 
  full_join(read.csv("/vol/projects/CIIM/Influenza/iMED/metadata/meta_cohort1.csv", 
                     header = TRUE) %>% select(-X))

# protein data
load("/vol/projects/CIIM/Influenza/iMED/proteomic/analysis/output/cohort1.RData")
protein_meta_T1 <- meta %>% filter(time == "T1")

iMEDpilot_proteinDat_T1 <- df %>% 
  rownames_to_column("name") %>% filter(name %in% protein_meta_T1$name) %>% 
  column_to_rownames("name") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("OlinkID") %>% 
  left_join(annot) %>% 
  filter(Assay %in% rownames(inputDat_protein)) %>%
  column_to_rownames("Assay") %>% select(-OlinkID, -UniProt)
rm(annot, df, meta, protein_meta_T1)

# metabolite data
load("/vol/projects/CIIM/Influenza/iMED/metabolic/analysis/brendan/output/cohort_1/cohort1.Rdata")
mebo_meta_T1 <- meta %>% filter(time == "T1")

iMEDpilot_meboDat_T1 <- df %>% 
  select(mebo_meta_T1$name) %>% 
  rownames_to_column("ionIdx") %>% mutate(ionIdx = as.numeric(ionIdx)) %>%
  left_join(annot2 %>% select(ionIdx, Formula)) %>%
  select(-ionIdx) %>%
  filter(Formula %in% overlapped_metabolites) %>%
  column_to_rownames("Formula")
rm(annot, annot2, df, hmdbnames, mebo_meta_T1)

## prepare iMEDpilot data for the model -----------------
inputDat_iMEDpilot <- iMEDpilot_proteinDat_T1 %>% rbind(iMEDpilot_meboDat_T1)

# iMED pilot data reclassification
iMEDpilot_ab2 <- iMEDpilot_ab %>%
  mutate(ab_H1N1 = ifelse(H1N1_T1 ==0, H1N1_T4, H1N1_T4/H1N1_T1)) %>%
  mutate(ab_H3N2 = ifelse(H3N2_T1 ==0, H3N2_T4, H3N2_T4/H3N2_T1))%>%
  mutate(ab_B = ifelse(B_T1 ==0, B_T4, B_T4/B_T1)) %>%
  mutate(across(where(is.numeric), ~log2(.x))) %>%
  mutate(across(.cols = everything(), ~ ifelse(is.infinite(.x), 0, .x)))

baseline <- log2(40)
abFC <- log2(4)

get.reclassify <- function(dat, baseline_col, abFC_col, baseline_cutoff, abFC_cutoff) {
  # Description: reclassify vaccine response based on the baseline and abFC
  
  # Arguments: 
  # dat: data with patient ID, HAI value per time point, abFC of HAI between time point per column
  # baseline_col, abFC_col: names of the baseline and abFC columns
  # baseline_cutoff, abFC_cutoff: threshold for the classification
  
  # Returns: 
  # 4 reclassified groups: HH - high baseline & high abFC, 
  # HL - high baseline & log abFC, LH - low baseline & high abFC, LL - low baseline & low abFC
  
  reclassify <- paste0(ifelse(dat[, baseline_col] >= baseline_cutoff, "H", "L"),
                       ifelse(dat[, abFC_col] >= abFC_cutoff, "H", "L"))
  return(reclassify)
}
iMED_strains <- c("H1N1", "H3N2", "B")
reclassify_temp <- list()
for (strain in iMED_strains) {
  reclassify_temp[[strain]] <- get.reclassify(dat = iMEDpilot_ab2, 
                                              baseline_col = paste0(strain, "_T1"),
                                              abFC_col = paste0("ab_", strain),
                                              baseline_cutoff = baseline,
                                              abFC_cutoff = abFC)
  
}

iMEDpilot_HAIreclassify <-  iMEDpilot_ab2 %>% 
  cbind(as.data.frame(reclassify_temp) %>% 
          rename_with(~paste0(.x, "_reclassify")))
rm(reclassify_temp)

meta_iMEDpilot <- iMEDpilot_HAIreclassify %>% 
  left_join(iMEDpilot_metadat %>% 
              select(-H1N1_T1, -H1N1_T4, -H3N2_T1, -H3N2_T4, -B_T1, -B_T4) %>%
              filter(time == "T1"))  %>%
  rename("sex" = "gender") %>%
  mutate(sex = ifelse(sex == "male", 0, 1))

meta_iMEDpilot %>% count(H1N1_reclassify) # 3 LLs
meta_iMEDpilot %>% count(H3N2_reclassify) # no LL group
meta_iMEDpilot %>% count(B_reclassify) # 8 LLs

## extract individual for LL and Protector (LH, HL, HH) groups for iMEDpilot cohort ---------------
dat_iMEDpilot <- meta_iMEDpilot %>% 
  full_join(inputDat_iMEDpilot %>% t() %>% as.data.frame %>% rownames_to_column("name"))

dat_iMEDpilot_H1N1 <-  dat_iMEDpilot %>% select(-matches("H3N2|B")) %>% # equal LL (n = 3) and protector (n = 31)
  mutate(reclassify = ifelse(H1N1_reclassify == "LL", "LL", "Protector"))
dat_iMEDpilot_H1N1 %>% count(reclassify)

dat_iMEDpilot_B <-  dat_iMEDpilot %>% select(-matches("H1N1|H3N2")) %>% # equal LL (n = 8) and protector (n = 26)
  mutate(reclassify = ifelse(B_reclassify == "LL", "LL", "Protector"))
dat_iMEDpilot_B %>% count(reclassify)

## Run elastic model --------------------------------
set.seed(123) #set the seed to make your partition reproducible

# split train-test set
trainSet <- dat_iMED %>% rename("ab_baseline" = "H1N1_T1_combine")
trainSet %>% count(reclassify)

validateSet_1 <- dat_ZirFlu %>% rename("ab_baseline" = "H1N1_T1_combine")
validateSet_1 %>% count(reclassify)

validateSet_2 <- dat_iMEDpilot_H1N1 %>% rename("ab_baseline" = "H1N1_T1")
validateSet_2 %>% count(reclassify)

validateSet_3 <- dat_iMEDpilot_B %>% rename("ab_baseline" = "B_T1")
validateSet_3 %>% count(reclassify)

# input variable
proteinNames <- rownames(inputDat_protein)
metaboliteNames <- rownames(inputDat_metabolite)

inputVariables <- list()
inputVariables$proteins <- proteinNames
inputVariables$age_proteins <- c("age", proteinNames)
inputVariables$age_sex_proteins <- c("age", "sex", proteinNames)
inputVariables$age_sex_abBaseline_proteins <- c("age", "sex", "ab_baseline", proteinNames)
inputVariables$age_sex_abBaseline_proteins_metabolites <- c("age", "sex", "ab_baseline", proteinNames, metaboliteNames)

res.elasticModel <- list()
for (inputVal in names(inputVariables)) {
  overlapped_parameters <- intersect(inputVariables[[inputVal]], names(validateSet_2))
  
  # prepare input and validation sets
  train_inputSet <- trainSet[, overlapped_parameters]
  train_predOut <- trainSet[, c("reclassify")] %>% as.vector() %>% unlist()
  
  validation_inputSet_1 <- validateSet_1[, overlapped_parameters]
  validation_predOut_1 <- validateSet_1[, c("reclassify")] %>% 
    as.vector() %>% unlist() %>% as.factor()
  
  # check the ovelap proteins & metabolites
  validation_inputSet_2 <- validateSet_2[, overlapped_parameters]
  validation_predOut_2 <- validateSet_2[, c("reclassify")] %>% 
    as.vector() %>% unlist() %>% as.factor()
  
  validation_inputSet_3 <- validateSet_3[, overlapped_parameters]
  validation_predOut_3 <- validateSet_3[, c("reclassify")] %>% 
    as.vector() %>% unlist() %>% as.factor()
  
  # run the model
  res.elasticModel[[inputVal]] <- get.elasticModel(train_inputSet, train_predOut, 
                                                   validation_inputSet_1, validation_predOut_1,
                                                   validation_inputSet_2, validation_predOut_2,
                                                   validation_inputSet_3, validation_predOut_3)
}

res.elasticModel_2groups_equalRatio_iMEDtraining_ZirFlu_iMEDpilot_valid_reduceParameters <- res.elasticModel

## ISSUE =====
#les parameter in the iMED pilot data (miss 261 metabolites)

# save the result
save(res.elasticModel_2groups_equalRatio_iMEDtraining_ZirFlu_iMEDpilot_valid_reduceParameters, 
     file = "20230503_res.elasticModel_2groups_equalRatio_iMEDtraining_ZirFlu_iMEDpilot_valid_reduceParameters.RData")
