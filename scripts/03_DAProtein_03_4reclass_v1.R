rm(list = ls())

library(tidyverse)
library(limma)
library(gplots)

get.limmaRes <- function(metaDat, inputDat) {
  # Aim: identify DE proteins/metabolites correct with sex, age, and reclassify (the interested vaccine response reclassification groups) 
  #by running the linear model (limma package)
  
  # input: metadata - interested participant (in "name" column) with sex, age, and reclassify group, 
  #inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  
  # output: res - limma output
  
  inputDat_temp <- inputDat[, metaDat$name]
  
  if (identical(metaDat$name, colnames(inputDat_temp)) == TRUE) {
    res <- lmFit(inputDat_temp,
                 design = model.matrix(~ + sex + age + reclassify, metaDat)) %>% 
      eBayes()
  } else res <- "Error: check input"
  
  return(res)
}

get.limmaRes_perStrain <- function(metadat, inputDat, strain_groups) {
  # Aim: run the linear model (using get.limaRes function) with for multiple strain
  
  # input: metadata - interested participant (in "name" column) with sex, age, and reclassify group, 
  # inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  # strain_groups - strains which the test apply to
  
  # output: res - list of outcome from limmar model per strains
  
  resList <- list()
  for (strain_group in strain_groups) {
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na()
    
    resList[[strain_group]] <- get.limmaRes(metadat_temp, inputDat)
  }
  return(resList)
}

# load data =======================================================================
load("cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

## run lhe limma model -------------------------
iMED_strains <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")
ZirFlu_strains <- c("H1N1_reclassify", "H3N2_reclassify", "Bvictoria_reclassify", "Byamagata_reclassify")

# iMED cohort 2014 
metadat_iMED_2014 <- metadata_healthy %>% filter(season == "2014") 
inputDat_iMED_2014 <- protein_Dat$iMED_2014 %>% t() %>% as.data.frame %>% select(metadat_iMED_2014$name)
# identical(colnames(inputDat_iMED_2014), metadat_iMED_2014$name) # TRUE, the same sample order

res_iMED_2014 <- get.limmaRes_perStrain(metadat =  metadat_iMED_2014,
                                   inputDat = inputDat_iMED_2014, 
                                   strain_groups = iMED_strains)

# iMED cohort 2015
metadat_iMED_2015 <- metadata_healthy %>% filter(season == "2015") 
inputDat_iMED_2015 <- protein_Dat$iMED_2015 %>% t() %>% as.data.frame %>% select(metadat_iMED_2015$name)

res_iMED_2015 <- get.limmaRes_perStrain(metadat =  metadat_iMED_2015,
                                        inputDat = inputDat_iMED_2015, 
                                        strain_groups = iMED_strains)

# ZirFlu cohort 2019
metadat_ZirFlu_2019 <- metadata_healthy %>% filter(season == "2019") 
inputDat_ZirFlu_2019 <- protein_Dat$ZirFlu_2019 %>% t() %>% as.data.frame %>% select(metadat_ZirFlu_2019$name)
# identical(colnames(inputDat_ZirFlu_2019), metadat_ZirFlu_2019$name) # TRUE, the same sample order

res_ZirFlu_2019 <- get.limmaRes_perStrain(metadat =  metadat_ZirFlu_2019,
                                        inputDat = inputDat_ZirFlu_2019, 
                                        strain_groups = ZirFlu_strains)

# ZirFlu cohort 2020
metadat_ZirFlu_2020 <- metadata_healthy %>% filter(season == "2020") 
inputDat_ZirFlu_2020 <- protein_Dat$ZirFlu_2020 %>% t() %>% as.data.frame %>% select(metadat_ZirFlu_2020$name)

res_ZirFlu_2020 <- get.limmaRes_perStrain(metadat =  metadat_ZirFlu_2020,
                                          inputDat = inputDat_ZirFlu_2020, 
                                          strain_groups = ZirFlu_strains)

# save data ------------------------------------------------
resPro_4reclass <- list(  "iMED_2014" = res_iMED_2014, 
                          "iMED_2015" = res_iMED_2015,
                          "ZirFlu_2019" = res_ZirFlu_2019, 
                          "ZirFlu_2020" = res_ZirFlu_2020)

save(resPro_4reclass, file = "resPro_4reclass.RData")
