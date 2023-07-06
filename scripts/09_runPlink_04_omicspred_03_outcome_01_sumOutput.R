rm(list = ls())

library(tidyverse)
# load model information =======================================================================
models <- c("Metabolon", "Nightingale", "Olink", "RNAseq", "Somalogic")
models_info <- list()
for (model_name in models) {
  models_info[[model_name]] <- read.csv2(paste0("data/omicsPred/" , model_name,
                                                "_trait_validation_results_with_OMICSPRED_ID.csv"), sep = "\t")
}

models_Tab <- models_info %>% bind_rows(.id = "model") %>% 
  mutate(type = ifelse(model == "Metabolon"| model == "Nightingale", "metabolite",
                       ifelse(model == "Olink"| model == "Somalogic", "protein", "RNAseq")))

## match metabolite with HMDB database --------------------------
hmdb_endogenous <- read_csv("data/20221011_HMDB_endogenousMetabolites")

models_Tab_mebo <- models_info[c("Metabolon", "Nightingale")] %>% bind_rows(.id = "model") %>% 
  mutate(Name = ifelse(is.na(Biochemical.Name), Biomarker.Name, Biochemical.Name), # merge the name, pathway/group information of 2 platform together
         subClass = ifelse(is.na(Sub.Pathway), Subgroup, Sub.Pathway),
         class = ifelse(is.na(Super.Pathway), Group, Super.Pathway),
         name = tolower(Name)) %>% 
  select(-c("Biomarker.Name", "Biochemical.Name", 
            "Subgroup", "Sub.Pathway", "Group", "Super.Pathway")) %>%
  relocate(c(Name, subClass, class), .before = "X.SNP") %>%
  relocate(c(Metabolon.ID, Trait.ID), .after = name) %>%
  left_join(hmdb_endogenous %>% mutate(name = tolower(NAME))) # annnote to the hmdb database using matched names

models_mebo_formula <- models_Tab_mebo %>% 
  select(model, name, CHEMICAL_FORMULA) %>% drop_na() %>% add_count(CHEMICAL_FORMULA)
models_mebo_formula %>% count(model)
models_mebo_formula %>% filter(n > 1)

# extract model predicted outcome ======================================================
models <- c("Metabolon", "Nightingale", "Olink", "RNAseq", "SomaScan") # SomaScan == Somalogic, OmicsPred use different names in model and summary files
predOut <- list()
for (model_name in models) {
  filesList <- list.files(paste0("processedDat/omicsPred/output/", model_name))
  
  for (fileName in filesList[grep("log", filesList)]) {
    if (file.size(paste0("processedDat/omicsPred/output/", model_name, "/", fileName)) > 0) {
      predOut$log[[substring(fileName, 1, 10)]] <- read.delim(
        paste0("processedDat/omicsPred/output/", model_name, "/", fileName), header = FALSE)}}
   
  
  for (fileName in filesList[grep("score", filesList)]) {
    if (grepl("sscore.vars", fileName)) {
      if (file.size(paste0("processedDat/omicsPred/output/", model_name, "/", fileName)) > 0) {
        predOut$var_rsid[[substring(fileName, 1, 10)]] <- read.delim(
          paste0("processedDat/omicsPred/output/", model_name, "/", fileName), header = FALSE)}
      } else {
        if (file.size(paste0("processedDat/omicsPred/output/", model_name, "/", fileName)) > 0) {
          predOut$score[[substring(fileName, 1, 10)]] <- read.delim(
            paste0("processedDat/omicsPred/output/", model_name, "/", fileName))}}
    }
}

save(predOut, file = "omicsPred_output.RData")

## prediction .log information  ----------------------------------------
pred_infoTab <- pred_info %>% 
  lapply(function(x) if(nrow(x) == 24) {x %>% add_row(V1 = NA, .before = 19)} else x) %>%
  imap(~.x %>% rename_with(function(x) .y)) %>% 
  purrr::reduce(cbind)

pred_model_info <- pred_infoTab[19:20, ] %>% t() %>% as.data.frame() %>% 
  rename("warning" = "19", "number_SNPs_use" = "20") %>%
  mutate(number_SNPs_use = gsub("--score: | variant processed.| variants processed.",
                                "", number_SNPs_use)) %>% 
  rownames_to_column("OMICSPRED.ID") %>% 
  full_join(selected_proPred %>% rename("number_SNPs_total" = "X.SNP")) %>% 
  relocate(number_SNPs_use, .after = "number_SNPs_total")

## prediction outcome  --------------------------------------------------------
predicted_dat <- model_outcome %>% 
  lapply(function(x) x %>% as.data.frame %>% column_to_rownames("X.IID")) %>%
  imap(~.x %>% rename_with(function(x) .y)) %>% 
  purrr::reduce(cbind) %>% 
  t() %>% as.data.frame %>% 
  rownames_to_column("OMICSPRED.ID") %>% full_join(pred_model_info %>% select(OMICSPRED.ID, Gene)) %>%
  column_to_rownames("Gene") %>% select(-OMICSPRED.ID) %>% t() %>% as.data.frame()



