rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data =======================================================================
proPred_dat <- read.csv2("data/omicsPred/Olink_trait_validation_results_with_OMICSPRED_ID.csv", sep = "\t")
load("selected_DAPs.RData")

selected_proPred <- proPred_dat %>% 
  filter(Gene %in% selected_DAs) %>% 
  select(OMICSPRED.ID, Gene, X.SNP, Internal_R2, NSPHS_R2, ORCADES_R2)

# extract file ----------------------------------------------------------
model_files <- list.files("processedDat/omicsPred/Olink_output")[
  which(substring(list.files("processedDat/omicsPred/Olink_output"), 1, 10) 
        %in% selected_proPred$OMICSPRED.ID)]

model_info <- list()
model_outcome <- list()
model_rsid <- list()
for (model_id in unique(substring(model_files, 1, 10))) {
  model_info[[model_id]] <- read.delim(paste0("processedDat/omicsPred/Olink_output/", 
                                             model_id, "_model.txt.log"), header = FALSE)
  model_outcome[[model_id]] <- read.delim(paste0("processedDat/omicsPred/Olink_output/", 
                                                model_id, "_model.txt.sscore"))
  model_rsid[[model_id]] <- read.delim(paste0("processedDat/omicsPred/Olink_output/", 
                                             model_id, "_model.txt.sscore.vars"), header = FALSE)
}

# get model running information -----------------------------
model_info_v2 <- model_info %>% 
  lapply(function(x) if(nrow(x) == 24) {x %>% add_row(V1 = NA, .before = 19)} else x) %>%
  imap(~.x %>% rename_with(function(x) .y)) %>% 
  purrr::reduce(cbind)

pred_model_info <- model_info_v2[19:20, ] %>% t() %>% as.data.frame() %>% 
  rename("warning" = "19", "number_SNPs_use" = "20") %>%
  mutate(number_SNPs_use = gsub("--score: | variant processed.| variants processed.",
                                "", number_SNPs_use)) %>% 
  rownames_to_column("OMICSPRED.ID") %>% 
  full_join(selected_proPred %>% rename("number_SNPs_total" = "X.SNP")) %>% 
  relocate(number_SNPs_use, .after = "number_SNPs_total")

# get model outcome --------------------------------------
predicted_dat <- model_outcome %>% 
  lapply(function(x) x %>% as.data.frame %>% column_to_rownames("X.IID")) %>%
  imap(~.x %>% rename_with(function(x) .y)) %>% 
  purrr::reduce(cbind) %>% 
  t() %>% as.data.frame %>% 
  rownames_to_column("OMICSPRED.ID") %>% full_join(pred_model_info %>% select(OMICSPRED.ID, Gene)) %>%
  column_to_rownames("Gene") %>% select(-OMICSPRED.ID) %>% t() %>% as.data.frame()

# check the correlation between prediction and real data -------------------
load("cohorts_dat.RData")

## protein measures from coresspoding subjects
measured_dat <- cohorts$donorSample_all %>% 
  filter(time == "T1") %>% select(probandID, name) %>%
  filter(probandID %in% rownames(predicted_dat)) %>%
  left_join(protein_Dat$iMED_2015 %>% rownames_to_column("name")) %>%
  arrange(match(probandID, rownames(predicted_dat))) %>% 
  column_to_rownames("probandID") %>% select(-name)

identical(rownames(predicted_dat), rownames(measured_dat)) # TRUE = same order of subjects

# calculate the r and p value
for (protein_name in pred_model_info$Gene) {
  out_temp <- cor.test(predicted_dat[, protein_name], measured_dat[, protein_name], 
                       use = "complete.obs")
  pred_model_info$cor[pred_model_info$Gene == protein_name] <- out_temp$estimate
  pred_model_info$pval[pred_model_info$Gene == protein_name] <- out_temp$p.value
}

# make scatter plot
protein <- "CXCL1"
protein <- "IL10RB"
protein <- "OSCAR"
protein <- "IL15RA"
protein <- "LTA"
protein <- "TNFRSF11A"
protein <- "FGF5"
protein <- "CD4"
protein <- "TNFSF10"
protein <- "VEGFA"
protein <- "LIFR"
protein <- "CCL7"
protein <- "CXCL8"
protein <- "OSM"
protein <- "IL6"
protein <- "SIGLEC1"

dat_temp <- cbind("predicted" = predicted_dat[, protein],
                  "measured" = measured_dat[, protein])

dat_temp %>% as.data.frame() %>% 
  ggplot(aes(x = measured, y = predicted)) + 
  geom_point() + 
  stat_smooth(method = "lm", col = "blue") + stat_cor() + 
  theme_classic() + ggtitle(protein)

# check rsid ------------------
protein <- "CXCL1"
write.table(model_rsid[[pred_model_info$OMICSPRED.ID[pred_model_info$Gene == protein]]],
            file = paste0("processedDat/omicsPred/", protein, "_rsid.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

