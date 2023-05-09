library(OlinkAnalyze)
library(openxlsx)
library(tidyverse)
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
hmdb_endogenous <- read_csv("/vol/projects/CIIM/Influenza/ZirrFlu/ZirFlu_NhanNguyen/reference/20221011_HMDB_endogenousMetabolites")

iMED_meboAnnot <- rawDat_iMED_meboAnnot %>% 
  slice(which(Formula %in% hmdb_endogenous$CHEMICAL_FORMULA))
length(unique(iMED_meboAnnot$ionIdx)) # 1345 metabolites

ZirFlu_meboAnnot <- rawDat_ZirFlu_meboAnnot %>% 
  slice(which(Formula %in% hmdb_endogenous$CHEMICAL_FORMULA))
length(unique(ZirFlu_meboAnnot$ionIdx)) # 786 metabolites

# Covert the raw metabolites data to table format =========================================================================
# iMED
iMED_mebo <- rawDat_iMED_mebo  %>% tibble::column_to_rownames('ionIdx') %>% 
  select(matches("human")) %>%
  t() %>% as.data.frame %>%
  select(unique(iMED_meboAnnot$ionIdx))

# ZirFlu
ZirFlu_mebo_temp <- rawDat_ZirFlu_mebo  %>% tibble::column_to_rownames('ionIdx') %>% 
  select(-all_of(c("ionMz", "ionAverageInt", "ionTopFormula", "ionTopIon", "ionTopName"))) %>%
  t() %>% as.data.frame %>%
  select(unique(ZirFlu_meboAnnot$ionIdx))

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
