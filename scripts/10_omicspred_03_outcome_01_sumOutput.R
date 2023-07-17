rm(list = ls())

library(tidyverse)
library(stringr)
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

models_Tab_v2 <- models_Tab %>% 
  mutate(Name = ifelse(is.na(Biochemical.Name), Biomarker.Name, Biochemical.Name), # merge the name, pathway/group information of 2 platform together
         subClass = ifelse(is.na(Sub.Pathway), Subgroup, Sub.Pathway),
         class = ifelse(is.na(Super.Pathway), Group, Super.Pathway),
         name = tolower(Name)) %>% 
  select(-c("Biomarker.Name", "Biochemical.Name", 
            "Subgroup", "Sub.Pathway", "Group", "Super.Pathway")) %>%
  relocate(c(Metabolon.ID, Trait.ID), .after = name) %>%
  left_join(hmdb_endogenous %>% mutate(name = tolower(NAME))) # annotate to the hmdb database using matched names

models_mebo_formula <- models_Tab_v2 %>% 
  select(model, name, CHEMICAL_FORMULA) %>% drop_na() %>% add_count(CHEMICAL_FORMULA)
models_mebo_formula %>% count(model)
models_mebo_formula %>% filter(n > 1)

# extract model predicted outcome ======================================================
# models <- c("Metabolon", "Nightingale", "Olink", "RNAseq", "SomaScan") # SomaScan == Somalogic, OmicsPred use different names in model and summary files
# predOut <- list()
# for (model_name in models) {
#   filesList <- list.files(paste0("processedDat/omicsPred/output/", model_name))
#   
#   for (fileName in filesList[grep("log", filesList)]) {
#     if (file.size(paste0("processedDat/omicsPred/output/", model_name, "/", fileName)) > 0) {
#       predOut$log[[substring(fileName, 1, 10)]] <- read.delim(
#         paste0("processedDat/omicsPred/output/", model_name, "/", fileName), header = FALSE)}}
#    
#   
#   for (fileName in filesList[grep("score", filesList)]) {
#     if (grepl("sscore.vars", fileName)) {
#       if (file.size(paste0("processedDat/omicsPred/output/", model_name, "/", fileName)) > 0) {
#         predOut$var_rsid[[substring(fileName, 1, 10)]] <- read.delim(
#           paste0("processedDat/omicsPred/output/", model_name, "/", fileName), header = FALSE)}
#       } else {
#         if (file.size(paste0("processedDat/omicsPred/output/", model_name, "/", fileName)) > 0) {
#           predOut$score[[substring(fileName, 1, 10)]] <- read.delim(
#             paste0("processedDat/omicsPred/output/", model_name, "/", fileName))}}
#     }
# }

#save(predOut, file = "omicsPred_output.RData")

## prediction .log information  ----------------------------------------
load("omicsPred_output.RData")

# logFile_size <- predOut$log %>% lapply(function(x) nrow(x)) %>% unlist() 
# unique(logFile_size) # 4 size: 20, 21 24, 25
# which(logFile_size == 20) #, check files with different file size, such as nrow = 20
# View(predOut$log$OPGS004208) 

pred_logTab <- predOut$log %>% 
  lapply(function(x) # make all logFile size to 25, by adding NA row to files have less than 25 rows
    if(nrow(x) == 24) {x %>% add_row(V1 = NA, .before = 19)
    } else if (nrow(x) == 21) {
      x %>% add_row(V1 = rep(NA, 4), .before = 21)} else if (nrow(x) == 20) {
        x %>% add_row(V1 = rep(NA, 5), .before = 20)} else x) %>%
  imap(~.x %>% rename_with(function(x) .y)) %>% 
  purrr::reduce(cbind)

# pred_logTab[c(1:4, 7:8, 11:18), ] %>% unlist() %>% unique() # return same value across columns, because it is about how the models run
# pred_logTab[24, ] %>% unlist() %>% unique() # return "." (the line from model log file) or NA (filled row in the previous step to match the universe 25 rows on log files)

pred_logInfo <- pred_logTab[19:20, ] %>% # extract the variants (SNP) used information from models'log file
  t() %>% as.data.frame() %>% 
  rename("warning" = "19", "SNPs_use_info" = "20") %>%
  mutate(num_SNPs_use = ifelse(SNPs_use_info == "Error: No valid variants in --score file.", 0, 
                               str_extract(SNPs_use_info, "[[:digit:]]+")) %>% as.numeric) %>%
  mutate(numSNPs_skip = ifelse(is.na(warning), 0, # NA = skip SNPs is not applicable
                               str_extract(warning, "[[:digit:]]+")) %>% as.numeric) %>% 
  rownames_to_column("OMICSPRED.ID")

## model prediction information  --------------------------------------------------------
models_pred_Tab <- models_Tab_v2 %>% 
  select(OMICSPRED.ID, X.SNP, model, type, Gene, Ensembl.ID,
         CHEMICAL_FORMULA, name, Name, subClass, class, 
         matches("R2|Rho|MissingRate")) %>%
  full_join(pred_logInfo) %>% relocate(c(warning, matches("SNPs")), .after = "X.SNP")

# select models to further analsyis ======================================================
selected_models <- models_pred_Tab %>% filter(num_SNPs_use >= 5) # the R2 can be good?

models_pred_Tab %>% count(type, model)
selected_models %>% count(type, model)

models_pred_Tab %>% 
  filter(type == "RNAseq") %>% 
#  filter(type != "RNAseq") %>% 
  filter(X.SNP < 10) %>%
  ggplot(aes(x = X.SNP, y = as.numeric(Internal_R2))) + 
  geom_point() + theme_classic()

models_pred_Tab %>% 
  filter(type == "RNAseq") %>% 
#  filter(type != "RNAseq") %>% 
  filter(num_SNPs_use < 10) %>%
  ggplot(aes(x = num_SNPs_use, y = as.numeric(Internal_R2))) + 
  geom_point() + theme_classic()

## load predicted value ---------------------------------
pred_Tab <- predOut$score %>%
  lapply(function(x) x %>% column_to_rownames("X.IID")) %>% 
  imap(~.x %>% rename_with(function(x) .y)) %>% 
  purrr::reduce(cbind) %>% select(selected_models$OMICSPRED.ID)

#save(pred_Tab, selected_models, file = "omicsPred_predTab.RData")

# OmicPred models overview - barplot ======================================================
selected_models <- models_pred_Tab %>% filter(num_SNPs_use >= 5) # the R2 can be good?
models_Tab %>% count(type)
models_Tab_v2 <- models_Tab %>% full_join(models_pred_Tab) %>%
  mutate(status = ifelse(num_SNPs_use >=5, "used", "non_used")) %>%
  select(model, status) %>% summarise()

models_Tab_v2 %>% ggplot(aes(x = model, y = ))
