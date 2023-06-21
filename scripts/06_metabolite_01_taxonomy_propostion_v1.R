rm(list = ls())
library(tidyverse)
library(openxlsx)

get.tstat <- function(resLimma) {
  # Aim: extract the t-statistic value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_tstat <- resLimma$t %>% as.data.frame %>% 
    dplyr::select(matches("reclassify"))
  
  return(res_tstat)
}

# load the data ----------------
load("cohorts_dat.RData")

# metabolite taxonomy with HMDB ids
mebo_taxo <- read.delim("/vol/projects/CIIM/Influenza/iMED/metabolic/db/hmdb/metabolite_HMDB_taxonomy.csv",
                        quote = "", header = TRUE)

# link the HMDB ids to Formulas
mebo_taxo_fomula <- mebo_taxo %>% rownames_to_column("CompoundID") %>% 
  inner_join(read.xlsx('/vol/projects/CIIM/Influenza/iMED/metabolic/raw_data/tables/DATA_CURATED_reformatted.xlsx',
                       sheet = 'annotation') %>% fill(ionIdx, .direction = "down")) %>% 
  filter(Formula %in% colnames(mebo_Dat$iMED_2014)) %>%
  select(Formula, super_class, class, sub_class) %>% distinct()

length(unique(mebo_taxo_fomula$Formula)) # 499 formulas out of 508 formulas

mebo_taxo_fomula %>% count(super_class)
mebo_taxo_fomula %>% count(class)

# super class 
mebo_taxo_fomula2 <- mebo_taxo_fomula %>% select(Formula, super_class) %>% distinct()
mebo_taxo_fomula2 <- mebo_taxo_fomula %>% select(Formula, class) %>% distinct()


