rm(list = ls())
library(ggpubr)
library(tidyverse)
library(rstatix)

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

# load data =======================================================================

## load data from iMED transcriptome ---------------
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
# View(transcriptomeList[[1]])
# View(transcriptomeList[[2]])
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1") # time T1 in trancriptome = d0
iMED_transcripDat <- transcriptomeList[[1]] %>% as.data.frame %>% 
  select(iMED_transcrip_T1$SampleName)


## metadata for all healthy subjects -------------------------
load("processedDat/cohorts_dat.RData")
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

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
save(resPro_2015_bulkRNAseq, file = "processedDat/resPro_2015_bulkRNAseq_4groups.RData")

# heatmap ------------------------------------------------
load("resPro_2015_bulkRNAseq_4groups.RData")
load("selected_DAPs_padj2015.RData")

pval_dat <- resPro_2015_bulkRNAseq %>% 
  lapply(function(x) x$p.value %>% as.data.frame %>% 
           select(matches("reclassify")) %>% 
           rownames_to_column("valName")) %>% 
  bind_rows(.id = "strain") %>% mutate(strain = gsub("_reclassify", "", strain)) %>%
  pivot_longer(matches("reclassify"), names_to = "reclassify", values_to = "p.value")

tstat_dat <- resPro_2015_bulkRNAseq %>% 
  lapply(function(x) x$t %>% as.data.frame %>% 
           select(matches("reclassify")) %>% 
           rownames_to_column("valName")) %>% 
  bind_rows(.id = "strain") %>% mutate(strain = gsub("_reclassify", "", strain)) %>%
  pivot_longer(matches("reclassify"), names_to = "reclassify", values_to = "relative_diff")

plotDat_DAPs_bulkRNAseq <- tstat_dat %>% 
  left_join(pval_dat) %>% 
  mutate(group = paste0(strain, "_", reclassify), 
         group = gsub("reclassify", "LLvs", group)) %>% 
  filter(valName %in% selected_DAs)

plotDat_DAPs_bulkRNAseq %>%
  ggplot(aes(x = group, y = valName, fill = relative_diff)) + 
  geom_tile() +
 # geom_text(aes(label = ifelse(is.na(p.value), NA, "*"))) +
  scale_fill_gradientn(limits = c(-5, 5), colors = c("blue", "white", "red"), na.value = "grey") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  ggtitle("Bulk RNAseq")