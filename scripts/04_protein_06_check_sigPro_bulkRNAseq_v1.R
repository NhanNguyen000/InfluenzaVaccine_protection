rm(list = ls())
library(ggpubr)

# load data from iMED transcriptome ---------------
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
# View(transcriptomeList[[1]])
# View(transcriptomeList[[2]])
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1")
iMED_transcripDat <- transcriptomeList[[1]] %>% as.data.frame %>% 
  select(iMED_transcrip_T1$SampleName)


# load data =======================================================================
load("cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}

get.limmaRes <- function(metaDat, inputDat) {
  # Aim: identify DE proteins/metabolites correct with sex, age, and reclassify (the interested vaccine response reclassification groups) 
  #by running the linear model (limma package)
  
  # input: metadata - interested participant (in "SampleName" column) with sex, age, and reclassify group, 
  #inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  
  # output: res - limma output
  
  inputDat_temp <- inputDat[, metaDat$SampleName]
  
  if (identical(metaDat$SampleName, colnames(inputDat_temp)) == TRUE) {
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
      select(SampleName, sex, age, reclassify) %>% drop_na()
    
    resList[[strain_group]] <- get.limmaRes(metadat_temp, inputDat)
  }
  return(resList)
}

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~convert_protectees(.x)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "protectee")))

## run the limma model -------------------------
iMED_strains <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")

# iMED cohort 2015
inputDat_iMED_2015 <- iMED_transcripDat %>% as.data.frame()

metadat_iMED_2015 <- metadata_healthy %>% filter(season == "2015", cohort == "iMED") %>%
  right_join(iMED_transcrip_T1, by = c("probandID" = "patientID")) %>% 
  arrange(factor(SampleName, levels = colnames(inputDat_iMED_2015)))

resPro_2015_bulkRNAseq <- get.limmaRes_perStrain(metadat =  metadat_iMED_2015,
                                        inputDat = inputDat_iMED_2015, 
                                        strain_groups = iMED_strains)

# save data ------------------------------------------------
save(resPro_2015_bulkRNAseq, file = "resPro_2015_bulkRNAseq.RData")
