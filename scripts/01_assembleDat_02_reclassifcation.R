rm(list = ls())

library(tidyverse)

get.reclassify <- function(dat, baseline_col, abFC_col, baseline_cutoff, abFC_cutoff) {
  # Description: reclassify vaccine response based on the baseline and abFC
  
  # Arguments: 
  # dat: data with patient ID, HAI value per time point, abFC of HAI between time point per column
  # baseline_col, abFC_col: names of the baseline and abFC columns
  # baseline_cutoff, abFC_cutoff: threshold for the classification
  
  # Returns: 
  # 4 reclassified groups: HH - high baseline & high abFC, 
  # HL - high baseline & log abFC, LH - low baseline & high abFC, LL - low baseline & low abFC
  
  reclassify_temp <- paste0(ifelse(dat[, baseline_col] >= baseline_cutoff, "H", "L"),
                       ifelse(dat[, abFC_col] >= abFC_cutoff, "H", "L"))
  reclassify <- ifelse(reclassify_temp == "NANA", NA, reclassify_temp)
  return(reclassify)
}

# load data =======================================================================
load("processedDat/cohorts_dat.RData")

## cutoff ------------------
# (baseline < 40 ) == d0_low. (abFC >= 4) == abFC_high
# (baseline >= 40) == d0_high. (abFC < 4) == abFC_low
baseline <- 40
abFC <- 4

strains <- c("H1N1", "H3N2", "B", "Bvictoria", "Byamagata")
reclassify_temp <- list()
for (strain in strains) {
  reclassify_temp[[strain]] <- get.reclassify(dat = cohorts$HAI_all, 
                                              baseline_col = paste0(strain, "_d0"),
                                              abFC_col = paste0(strain, "_abFC"),
                                              baseline_cutoff = baseline,
                                              abFC_cutoff = abFC)
  
}

cohorts$HAI_all <-  cohorts$HAI_all %>% 
  cbind(as.data.frame(reclassify_temp) %>% 
          rename_with(~paste0(.x, "_reclassify")))
# save data ------------------------------------------------
save(cohorts, protein_Dat, mebo_Dat, file = "processedDat/cohorts_dat.RData")
