# load libraries & functions --------------------------
library(tidyverse)

get.selected_sampleDat <- function(dat, selected_sample) {
  # Description: subset the data for only selected samples
  
  # Arguments: 
  # dat - whole data with sample name in rownames, variable in column
  # selected_sample - name of the selected samples
  
  # Returns: data for only selected samples
  
  outcome <- dat %>% rownames_to_column("sample") %>%
    filter(sample %in% selected_sample) %>%
    column_to_rownames("sample") %>% as.data.frame()
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


# get protein data per (cohort & season) at baseline (T1) ---------------------------------
ZirFlu_sampleT1 <- list()
for (year in c("2019", "2020")) {
  sample_temp <- ZirFlu$donorSamples %>% filter(season == year)
  ZirFlu_sampleT1[[year]] <- ZirFlu$donorSamples %>% 
    filter(season == year) %>% filter(time == "T1")
}

iMED_sampleT1 <- iMED$donorSample %>% filter(time == "T1")

proteinDat <- list()
proteinDat[["ZirFlu_2019"]] <- get.selected_sampleDat(dat = protein_normDat$ZirFlu, 
                                                      selected_sample = ZirFlu_sampleT1$`2019`$probenID)
proteinDat[["ZirFlu_2020"]] <- get.selected_sampleDat(dat = protein_normDat$ZirFlu, 
                                                      selected_sample = ZirFlu_sampleT1$`2020`$probenID)
proteinDat[["iMED"]] <- get.selected_sampleDat(dat = protein_normDat$iMED, 
                                               selected_sample = iMED_sampleT1$name)

metadat_sampleT1 <- list()
metadat_sampleT1[["ZirFlu_2019"]] <- ZirFlu_sampleT1$`2019`
metadat_sampleT1[["ZirFlu_2020"]] <- ZirFlu_sampleT1$`2020`
metadat_sampleT1[["iMED"]] <- iMED_sampleT1
rm(iMED_sampleT1, ZirFlu_sampleT1)

# impute protein data --------------------------------
# check NAs percentage
cutoff <- 0.2
rmProteins <- list()
for (cohort in names(proteinDat)) {
  rmProteins[[cohort]] <- get.rmProteins(proteinDat[[cohort]], NA_cutoff = cutoff)
}

proteinDat_impute <- list()
for (cohort in names(proteinDat)) {
  proteinDat_impute[[cohort]] <- proteinDat[[cohort]] %>%
    dplyr::select(-c(names(rmProteins[[cohort]]))) %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))
}

save(proteinDat_impute, metadat_sampleT1, file = "proteinDat_impute.RData")
