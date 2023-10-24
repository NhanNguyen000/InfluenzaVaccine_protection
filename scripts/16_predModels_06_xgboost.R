rm(list = ls())

library(tidyverse)
library(magrittr)
library(caret)
library(xgboost)

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}

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

# load data =======================================================================
load("/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~convert_protectees(.x)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "protectee"))) %>%
  mutate(sex = ifelse(sex == "m", 0, 1))

# prepare input data ------------------------------------------------------
inputDat <- protein_Dat %>% purrr::reduce(rbind) %>% 
  cbind(mebo_Dat %>% purrr::reduce(rbind)) %>%
  as.data.frame %>% rownames_to_column("name")

dat_temp <- metadata_healthy %>% left_join(inputDat)

dat_temp %>% count(season, H1N1_reclassify) # LL group size >= 3 in all 4 seasons -> can use all 4 seasons
dat_temp %>% count(season, H3N2_reclassify) # LL group size >= 3 in season 2015-> can use season 2015
dat_temp %>% count(season, B_reclassify) # LL group size >= 3 in all 2 seasons 2014 and 2015 -> can use all 2 seasons
dat_temp %>% count(season, Bvictoria_reclassify) # LL group size < 3 in all 2 seasons 2019 and 2020 --> can not use 
dat_temp %>% count(season, Byamagata_reclassify) # LL group size >= 3 in seasons 2020 --> can use season 2020

## prepare train-test --------------------------------
trainSet <- dat_temp %>% filter(season == "2015") %>% 
  rename("reclassify" = "H1N1_reclassify", "ab_d0_log2" = "H1N1_d0_log2")

trainSet %>% count(reclassify)


## prepare validation-test --------------------------------
valiSets <- list()
# for H1N1
for (year in c("2014", "2019", "2020")) {
  valiSets[[paste0("H1N1_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "H1N1_reclassify", "ab_d0_log2" = "H1N1_d0_log2")
}

# for H3N2
valiSets$H3N2_2015 <- dat_temp %>% 
  filter(season == "2015") %>% 
  rename("reclassify" = "H3N2_reclassify", "ab_d0_log2" = "H3N2_d0_log2")

# for B
for (year in c("2014", "2015")) {
  valiSets[[paste0("B_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "B_reclassify", "ab_d0_log2" = "B_d0_log2")
}

# for Byamagata
valiSets$Byamagata<- dat_temp %>% 
  filter(season == "2020") %>% 
  rename("reclassify" = "Byamagata_reclassify", "ab_d0_log2" = "Byamagata_d0_log2")


valiSets %>% lapply(function(x) x%>% count(reclassify))

## input variable --------------------------------
proNames <- colnames(protein_Dat$iMED_2014)
meboNames <- colnames(mebo_Dat$ZirFlu_2019)

inputVars <- c("age", "sex", "ab_d0_log2", proNames, meboNames)

# prepare input and validation sets
train_inputSet <- trainSet[, inputVars] %>% as.matrix()

train_predOut <- trainSet %>% 
  mutate(reclassify_numb = ifelse(reclassify == "LL", 0, 1)) %>%
  select(reclassify, reclassify_numb)

validation_inputSets <- valiSets %>% 
  lapply(function(x) x[, inputVars] %>% as.matrix())

validation_predOuts <- valiSets %>% 
  lapply(function(x)  x %>% 
           mutate(reclassify_numb = ifelse(reclassify == "LL", 0, 1)) %>%
           select(reclassify, reclassify_numb))


# Run models (run in Slum/server) ==========================================
set.seed(123) #set the seed to make your partition reproducible

nreps <- 10
netFits <- list()
netFit_pred <- list()
for (i in 1:nreps) {
  netFits[[i]] <- xgboost(data = train_inputSet, label = train_predOut$reclassify_numb, 
                       max.depth = 2, eta = 1, nthread = 2, nrounds = 2, 
                       objective = "binary:logistic")

  for (valiSet in names(validation_inputSets)) {
    val.pred_temp <- predict(netFits[[i]], newdata = validation_inputSets[[valiSet]]) %>%  
      round() %>% as.data.frame() %>% rename("reclassify_numb" = ".") %>%
      mutate(reclassify = ifelse(reclassify_numb == 0, "LL", "protectee"),
             reclassify = factor(reclassify, levels = c("LL", "protectee")))
    
    netFit_pred[[valiSet]][[i]] <- confusionMatrix(val.pred_temp$reclassify, 
                                                   validation_predOuts[[valiSet]]$reclassify) # Training accuracy
  }
  
  print(paste0("Run ", i, " time."))
}

# save the result -----------------------
save(netFits, netFit_pred, file = "res.predModel_xgboost.RData")
