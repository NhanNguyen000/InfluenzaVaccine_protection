rm(list = ls())

library(tidyverse)

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
load("processedDat/cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
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

# prepare input data --------------------------------------------------------------
inputDat <- proteinDat_impute %>% purrr::reduce(rbind) %>% 
  cbind(mebo_Dat %>% purrr::reduce(rbind)) %>%
  as.data.frame %>% rownames_to_column("name")

dat_temp <- metadata_healthy %>% left_join(inputDat)

dat_temp %>% count(season, H1N1_reclassify) # LL group size >= 3 in all 4 seasons -> can use all 4 seasons
dat_temp %>% count(season, H3N2_reclassify) # LL group size >= 3 in season 2015-> can use season 2015
dat_temp %>% count(season, B_reclassify) # LL group size >= 3 in all 2 seasons 2014 and 2015 -> can use all 2 seasons
dat_temp %>% count(season, Bvictoria_reclassify) # LL group size < 3 in all 2 seasons 2019 and 2020 --> can not use 
dat_temp %>% count(season, Byamagata_reclassify) # LL group size >= 3 in seasons 2020 --> can use season 2020

# save the input (strain and validation set) for the prediction model ----------------

## input variable --------------------------------------------------------------
proNames <- colnames(proteinDat_impute$iMED_2014)
meboNames <- colnames(mebo_Dat$ZirFlu_2019)

## input with H1N1 season 2015 is train set------------------------------------------
trainSet <- dat_temp %>% filter(season == "2015") %>% 
  rename("reclassify" = "H1N1_reclassify", "ab_d0" = "H1N1_d0")
trainSet %>% count(reclassify)

valiSets <- list()
# for H1N1
for (year in c("2014", "2019", "2020")) {
  valiSets[[paste0("H1N1_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "H1N1_reclassify", "ab_d0" = "H1N1_d0")
}

# for H3N2
valiSets$H3N2_2015 <- dat_temp %>% 
  filter(season == "2015") %>% 
  rename("reclassify" = "H3N2_reclassify", "ab_d0" = "H3N2_d0")

# for B
for (year in c("2014", "2015")) {
  valiSets[[paste0("B_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "B_reclassify", "ab_d0" = "B_d0")
}

# for Byamagata
valiSets$Byamagata<- dat_temp %>% 
  filter(season == "2020") %>% 
  rename("reclassify" = "Byamagata_reclassify", "ab_d0" = "Byamagata_d0")

valiSets %>% lapply(function(x) x%>% count(reclassify))

### save the input ---------------------------------------------------------
save(trainSet, valiSets, proNames, meboNames, file = "processedDat/predictInput_H1N1trainSet.RData")

## input with H3N2 season 2015 is train set------------------------------------------
rm(trainSet, valiSets)

trainSet <- dat_temp %>% filter(season == "2015") %>% 
  rename("reclassify" = "H3N2_reclassify", "ab_d0" = "H3N2_d0")
trainSet %>% count(reclassify)

valiSets <- list()
# for H1N1
for (year in c("2014", "2015", "2019", "2020")) {
  valiSets[[paste0("H1N1_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "H1N1_reclassify", "ab_d0" = "H1N1_d0")
}

# for B
for (year in c("2014", "2015")) {
  valiSets[[paste0("B_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "B_reclassify", "ab_d0" = "B_d0")
}

# for Byamagata
valiSets$Byamagata<- dat_temp %>% 
  filter(season == "2020") %>% 
  rename("reclassify" = "Byamagata_reclassify", "ab_d0" = "Byamagata_d0")

valiSets %>% lapply(function(x) x%>% count(reclassify))

### save the input ---------------------------------------------------------
save(trainSet, valiSets, proNames, meboNames, file = "processedDat/predictInput_H3N2trainSet.RData")

## input with B season 2015 is train set------------------------------------------
rm(trainSet, valiSets)

trainSet <- dat_temp %>% filter(season == "2015") %>% 
  rename("reclassify" = "B_reclassify", "ab_d0" = "B_d0")
trainSet %>% count(reclassify)

valiSets <- list()
# for H1N1
for (year in c("2014", "2015", "2019", "2020")) {
  valiSets[[paste0("H1N1_", year)]] <- dat_temp %>% 
    filter(season == year) %>% 
    rename("reclassify" = "H1N1_reclassify", "ab_d0" = "H1N1_d0")
}

# for H3N2
valiSets$H3N2_2015 <- dat_temp %>% 
  filter(season == "2015") %>% 
  rename("reclassify" = "H3N2_reclassify", "ab_d0" = "H3N2_d0")

# for B
valiSets$B_2014 <- dat_temp %>% 
  filter(season == "2014") %>% 
  rename("reclassify" = "B_reclassify", "ab_d0" = "B_d0")

# for Byamagata
valiSets$Byamagata<- dat_temp %>% 
  filter(season == "2020") %>% 
  rename("reclassify" = "Byamagata_reclassify", "ab_d0" = "Byamagata_d0")

valiSets %>% lapply(function(x) x%>% count(reclassify))

### save the input ---------------------------------------------------------
save(trainSet, valiSets, proNames, meboNames, file = "processedDat/predictInput_BtrainSet.RData")
