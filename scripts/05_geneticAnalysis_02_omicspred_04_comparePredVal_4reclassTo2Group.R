rm(list = ls())

library(tidyverse)

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}

get.ttestRes <- function(dat_temp) {
  # Aim: identify DE molecules between reclassification groups (LL vs. other) by t-test
  
  # input: dat_temp with "model_id", "reclassify" column (column 1 and 2), and the rest of the columns are each separate features 
  # output: res - t-test output
  
  outcome <- data.frame(
    model_id = names(dat_temp)[-c(1,2)],
    t_statistic = NA,
    df = NA,
    p_value = NA,
    stderr = NA,
    alternative = NA,
    method = NA)
  
  for (model_id in names(dat_temp)[-c(1,2)]) {
    
    t_test <- t.test(
      dat_temp  %>% filter(reclassify == "LL") %>% pull(!!sym(model_id)),
      dat_temp  %>% filter(reclassify == "protectee") %>% pull(!!sym(model_id)), 
      var.equal = FALSE)  # Welch's t-test
    
    outcome[outcome$model_id == model_id, 2:7] <- c(
      t_test$statistic, t_test$parameter, t_test$p.value, 
      t_test$stderr, t_test$alternative, t_test$method)
  }
  
  return(outcome)
}

get.ttestRes_perStrain <- function(metadat, pred_Tab, strain_groups) {
  # Aim: run the linear model (using get.limaRes function) with for multiple strain
  
  # input: metadata - interested participant  with reclassify groups, 
  # pred_Tab - data with protein/metabolites (features) in colnames and participant in rownames
  # strain_groups - strains which the test apply to
  
  # output: res - list of outcome from ttest model per strains
  
  resList <- list()
  
  for (strain_group in strain_groups) {
    metadat_temp <- metadat %>% dplyr::select(c("probandID", strain_group)) %>% 
      dplyr::rename("reclassify" := strain_group) %>%
      full_join(pred_Tab %>% rownames_to_column("probandID"))
    
    resList[[strain_group]] <- get.ttestRes(metadat_temp) %>% 
      mutate(padj = p.adjust(p_value, method = "fdr")) 
  }
  return(resList)
}

# load & prepare data ===========================================
load("processedDat/cohorts_dat.RData")
load("processedDat/omicsPred_output.RData")

## metadata for corresponding subjects (who have genetic data) -------------------------
metadat <- cohorts$HAI_all %>% filter(probandID %in% rownames(pred_Tab)) %>%
  left_join(cohorts$donorInfo_all %>% dplyr::select(probandID, season, cohort, sex, age)) %>% 
  mutate_at(vars(contains("reclassify")), ~convert_protectees(.x)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "protectee")))

metadat %>% count(H1N1_reclassify)
metadat %>% count(H3N2_reclassify)
metadat %>% count(B_reclassify)

# run the t-test model (in this case  Welch's t-test) ----------------------------
strains <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")

res_omicsPred <- get.ttestRes_perStrain(metadat =  metadat,
                                        pred_Tab = pred_Tab, 
                                        strain_groups = strains)

# save data ------------------------------------------------------------------------
save(res_omicsPred, file = "processedDat/res_omicsPred_4reclass_2group.RData")


DEs <- res_omicsPred %>% 
  lapply(function(x) x %>% dplyr::select(model_id, p_value, padj))

tstat_dat <- res_omicsPred %>% 
  lapply(function(x) dplyr::select(model_id, t_statistic))
