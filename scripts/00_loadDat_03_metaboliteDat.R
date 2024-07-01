rm(list = ls())

library(OlinkAnalyze)
library(openxlsx)
library(tidyverse)
library(dplyr)
# load data =========================================================================
# iMED data
rawDat_iMED_mebo <- read.xlsx('/vol/projects/CIIM/Influenza/iMED/metabolic/raw_data/tables/DATA_CURATED_reformatted.xlsx',
                               sheet = 'ions')
rawDat_iMED_meboAnnot <- read.xlsx('/vol/projects/CIIM/Influenza/iMED/metabolic/raw_data/tables/DATA_CURATED_reformatted.xlsx',
                                     sheet = 'annotation') %>% fill(ionIdx, .direction = "down")

# ZirFlu
rawDat_ZirFlu_mebo <- read.xlsx(xlsxFile = '/vol/projects/CIIM/Influenza/ZirrFlu/metabolic/raw_data/spreadsheets/DATA.xlsx',
                                sheet = 'ions')
rawDat_ZirFlu_meboAnnot <- read.xlsx(xlsxFile = '/vol/projects/CIIM/Influenza/ZirrFlu/metabolic/raw_data/spreadsheets/DATA.xlsx',
                                     sheet = 'annotation') %>% fill(ionIdx, .direction = "down")

# Annotation: filter HMDB endogenous metabolites =========================================================================
# HMDB endogenous metabolites 
hmdb_endogenous <- read_csv("data/20221011_HMDB_endogenousMetabolites")

# iMED
iMED_meboAnnot_temp <- rawDat_iMED_meboAnnot %>% 
  slice(which(Formula %in% hmdb_endogenous$CHEMICAL_FORMULA))
length(unique(iMED_meboAnnot_temp$ionIdx)) # 1345 metabolites ionIdx
length(unique(iMED_meboAnnot_temp$Formula)) # 1400 metabolites formula (1 ionIdx - multiple formulas)

iMED_meboAnnot <- iMED_meboAnnot_temp %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_")) # combine multiple formulas per 1 ionIdx

iMED_meboAnnot$Formula_v2 <- NA # 1 formula - multiple ionIdxs --> add duplicated count to the duplicate formula
for (i in 1:nrow(iMED_meboAnnot)) {
  iMED_meboAnnot$Formula_v2[i] <- ifelse(i %in% which(duplicated(iMED_meboAnnot$Formula)), 
                                         paste0(iMED_meboAnnot$Formula[i], "_2"), 
                                         iMED_meboAnnot$Formula[i])
}

# ZirFlu
ZirFlu_meboAnnot_temp <- rawDat_ZirFlu_meboAnnot %>% 
  slice(which(Formula %in% hmdb_endogenous$CHEMICAL_FORMULA))
length(unique(ZirFlu_meboAnnot_temp$ionIdx)) # 786 metabolites
length(unique(ZirFlu_meboAnnot_temp$Formula)) # 799 metabolites formula (1 ionIdx - multiple formulas)

ZirFlu_meboAnnot <- ZirFlu_meboAnnot_temp %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_")) # combine multiple formulas per 1 ionIdx
which(duplicated(ZirFlu_meboAnnot$Formula)) # interge = 0 --> no case of 1 formula - multiple ionIdxs

# Covert the raw metabolites data to table format =========================================================================
# iMED
iMED_mebo <- rawDat_iMED_mebo %>%
  select(ionIdx, matches("human")) %>% right_join(iMED_meboAnnot) %>%
  select(-ionIdx, -Formula_v2) %>%
  group_by(Formula) %>% summarise_all("mean") %>% # calclulate the average for all ionIdx measures with identical Formula
  column_to_rownames("Formula") %>% 
  t() %>% as.data.frame()

# ZirFlu
ZirFlu_mebo_temp <- rawDat_ZirFlu_mebo  %>%
  select(-all_of(c("ionMz", "ionAverageInt", "ionTopFormula", "ionTopIon", "ionTopName"))) %>%
  right_join(ZirFlu_meboAnnot) %>%
  select(-ionIdx) %>% column_to_rownames("Formula") %>% 
  t() %>% as.data.frame()

old_probenID <- c(339151941, 339156196, 339151948, 339156278, 339151926, 339152850,
                  339156227, 339156287, 339156214, 339156220, 339156213, 339156314,
                  339156322, 339152826, 339152815, 339152794, 339156212, 339151960)
correct_probenID <- c(339151988, 339156157, 339151947, 339156279, 339151924, 339152857,
                      339156180, 339156286, 339156199, 339156204, 339156219, 339156288,
                      339156321, 339152802, 339152806, 339152807, 339156181, 339152000)
change_probenID <- data.frame("old_probenID" = old_probenID, "correct_probenID" = correct_probenID)
for (i in 1:nrow(change_probenID)) {
  rownames(ZirFlu_mebo_temp)[which(rownames(ZirFlu_mebo_temp) == change_probenID$old_probenID[i])] <- change_probenID$correct_probenID[i]
}

ZirFlu_mebo <-  ZirFlu_mebo_temp
rownames(ZirFlu_mebo) <- paste0("0", rownames(ZirFlu_mebo))

# save data ----------------------------------------------------
save(iMED_mebo, iMED_meboAnnot, ZirFlu_mebo, ZirFlu_meboAnnot, 
     file = "processedDat/metaboliteDat.RData")
