rm(list = ls())

library(tidyverse)
library(ggpubr)
library(org.Hs.eg.db)
# load data =======================================================================
load("omicsPred_predTab.RData")
load("cohorts_dat.RData")

pred_dat <- pred_Tab %>% 
  t() %>% as.data.frame %>% rownames_to_column("OMICSPRED.ID") %>% 
  full_join(selected_models %>% dplyr::select(OMICSPRED.ID, model, type, Gene, Ensembl.ID))

## protein data ---------------------------------------
pred_proDat <- pred_dat %>% filter(model == "Olink") %>%
  dplyr::select(-OMICSPRED.ID, -model, -type, -Ensembl.ID) %>% 
  column_to_rownames("Gene") %>% t() %>% as.data.frame()

# for FCRL3 (OPGS003959 - RNAseq, OPGS000406 - protein), TNFSF12(OPGS002615 - protein, OPGS006022 - RNAseq) 

predDat_v2 <- pred_dat %>% as.data.frame() %>%
  filter(OMICSPRED.ID %in% c("OPGS003959", "OPGS000406", "OPGS002615", "OPGS006022")) %>%
  mutate(valName = ifelse(OMICSPRED.ID == "OPGS003959", "FCRL3_RNAseq",
                          ifelse(OMICSPRED.ID == "OPGS000406", "FCRL3_protein", 
                                 ifelse(OMICSPRED.ID == "OPGS006022", "TNFSF12_RNAseq", "TNFSF12_protein")))) %>%
  dplyr::select(-OMICSPRED.ID, -model, -type, -Ensembl.ID, -Gene) %>% 
  column_to_rownames("valName") %>% t() %>% as.data.frame()

# protein measures from coresspoding subjects
measured_proDat <- cohorts$donorSample_all %>% 
  filter(time == "d0") %>% dplyr::select(probandID, name) %>%
  filter(probandID %in% rownames(pred_proDat)) %>%
  left_join(protein_Dat$iMED_2015 %>% rownames_to_column("name")) %>%
  arrange(match(probandID, rownames(pred_proDat))) %>% 
  column_to_rownames("probandID") %>% dplyr::select(-name) %>% as.data.frame()

identical(rownames(pred_proDat), rownames(measured_proDat)) # TRUE = same order of subjects


measuredDat_v2 <- cohorts$donorSample_all %>% 
  filter(time == "d0") %>% dplyr::select(probandID, name) %>%
  filter(probandID %in% rownames(predDat_v2)) %>%
  left_join(protein_Dat$iMED_2015 %>% rownames_to_column("name")) %>%
  arrange(match(probandID, rownames(predDat_v2))) %>% 
  column_to_rownames("probandID") %>% dplyr::select(-name) %>% as.data.frame()

identical(rownames(predDat_v2), rownames(measuredDat_v2)) # TRUE = same order of subjects

# calculate the correlation and p value
matched_proteins <- intersect(names(pred_proDat), names(measured_proDat))
pred_proInf <- selected_models %>% 
  filter(model == "Olink", Gene %in% matched_proteins) %>% 
  dplyr::select(-Ensembl.ID, -CHEMICAL_FORMULA, -name, -Name, -subClass, -class)

for (protein in matched_proteins) {
  out_temp <- cor.test(pred_proDat[, protein], measured_proDat[, protein], 
                       use = "complete.obs")
  pred_proInf$cor[pred_proInf$Gene == protein] <- out_temp$estimate
  pred_proInf$pval[pred_proInf$Gene == protein] <- out_temp$p.value
}

pred_proInf_shortView <- pred_proInf %>% 
  dplyr::select(OMICSPRED.ID, Gene, SNPs_use_info, num_SNPs_use, numSNPs_skip,
         Internal_R2, NSPHS_R2, ORCADES_R2, cor, pval)

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


dat_temp <- cbind("predicted" = pred_proDat[, protein],
                  "measured" = measured_proDat[, protein])
dat_temp %>% as.data.frame() %>% 
  ggplot(aes(x = measured, y = predicted)) + 
  geom_point() + 
  stat_smooth(method = "lm", col = "blue") + stat_cor() + 
  theme_classic() + ggtitle(protein)


protein <- "FCRL3"
protein <- "TNFSF12"
dat_temp <- cbind("predicted_Somalogic" = predDat_v2[, paste0(protein, "_protein")],
                  "measured_Olink" = measuredDat_v2[, protein])

dat_temp %>% as.data.frame() %>% 
  ggplot(aes(x = measured_Olink, y = predicted_Somalogic)) + 
  geom_point() + 
  stat_smooth(method = "lm", col = "blue") + stat_cor() + 
  theme_classic() + ggtitle(protein)

## RNAseq data ---------------------------------------
# RNAseq data of 20 donors in season 2015
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1")

# RNAseq prediction based on the genetic data of 159 donors in season 2015
matched_IDs <- intersect(iMED_transcrip_T1$patientID, names(pred_dat)) # 13 donors
pred_rnaDat <- pred_dat %>% 
  filter(model == "RNAseq") %>% 
  dplyr::select(c("OMICSPRED.ID", "model", "type", "Gene", "Ensembl.ID", matched_IDs)) %>% # Note: some row have Ensemble ID but no gene name
  left_join(mapIds(org.Hs.eg.db, keys = pred_dat$Ensembl.ID, # add missing gene name using org.Hs.eg.db with Ensemble ID
                   column = "SYMBOL", keytype =  "ENSEMBL") %>% 
              unlist() %>% as.data.frame() %>% 
              dplyr::rename("Gene2" = ".") %>% rownames_to_column("Ensembl.ID")) %>%
  mutate(Gene = ifelse(Gene == "", Gene2, Gene)) %>% 
  dplyr::select(-OMICSPRED.ID, -model, -type, -Ensembl.ID, -Gene2) %>% 
  drop_na() %>% # remove genes have no gene name and could not convert Ensembl.ID to gene name
  column_to_rownames("Gene") %>% t() %>% as.data.frame()

# RNAseq measures from 13 coresspoding subjects in season 2015
measured_rnaDat <- transcriptomeList[[1]] %>% 
  t() %>% as.data.frame %>% rownames_to_column("SampleName") %>% 
  right_join(iMED_transcrip_T1 %>% dplyr::select(SampleName, patientID)) %>%
  filter(patientID %in% rownames(pred_rnaDat)) %>%
  arrange(match(patientID, rownames(pred_rnaDat))) %>% 
  column_to_rownames("patientID") %>% dplyr::select(-SampleName)

identical(rownames(measured_rnaDat), rownames(pred_rnaDat)) # TRUE = same order of subjects

# calculate the correlation and p value
matched_genes <- intersect(names(pred_rnaDat), names(measured_rnaDat))

pred_rnaInf <- selected_models %>% 
  filter(model == "RNAseq", Gene %in% matched_genes) %>% 
  dplyr::select(-CHEMICAL_FORMULA, -name, -Name, -subClass, -class)

for (gene in matched_genes) {
  out_temp <- cor.test(pred_rnaDat[, gene], measured_rnaDat[, gene], 
                       use = "complete.obs")
  pred_rnaInf$cor[pred_rnaInf$Gene == gene] <- out_temp$estimate
  pred_rnaInf$pval[pred_rnaInf$Gene == gene] <- out_temp$p.value
}

pred_rnaInf_shortView <- pred_rnaInf %>% 
  dplyr::select(OMICSPRED.ID, Gene, SNPs_use_info, num_SNPs_use, numSNPs_skip,
         Internal_R2, INTERVAL_Withheld_R2, cor, pval)

# make scatter plot
gene <- "CD83" # no data in the selcted models
gene <- "NBN"

gene <- "CYP26B1"
gene <- "SLC12A7"

gene <- "RIPOR1"
gene <- "CRYL1"

gene <- "FCRL3"
gene <- "TNFSF12"

dat_temp <- cbind("predicted" = pred_rnaDat[, gene],
                  "measured" = measured_rnaDat[, gene])

dat_temp %>% as.data.frame() %>% 
  ggplot(aes(x = measured, y = predicted)) + 
  geom_point() + 
  stat_smooth(method = "lm", col = "blue") + stat_cor() + 
  theme_classic() + ggtitle(gene)

# sum up the cor(predicted, measured values) for protein Olink and RNAseq data ------------------
pred_inf <- pred_rnaInf %>% full_join(pred_proInf)

pred_inf %>% 
  ggplot(aes(x = cor, y = -log10(pval), col = model)) + 
  geom_jitter(size = 1, width = 0.15, height = 0.15) +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  theme_classic()

save(pred_inf, file = "pred_with_proRna_measures.RData")

# extract rsids per prediction into txt file ----------------------------------------------------------
load("omicsPred_predTab.RData")

model_files <- list.files("processedDat/omicsPred/Olink_output")[
  which(substring(list.files("processedDat/omicsPred/Olink_output"), 1, 10) 
        %in% selected_proPred$OMICSPRED.ID)]

model_rsid <- list()
for (model_id in unique(substring(model_files, 1, 10))) {
  model_rsid[[model_id]] <- read.delim(paste0("processedDat/omicsPred/Olink_output/", 
                                             model_id, "_model.txt.sscore.vars"), header = FALSE)
}

protein <- "CXCL1"
write.table(model_rsid[[pred_proInf$OMICSPRED.ID[pred_model_info$Gene == protein]]],
            file = paste0("processedDat/omicsPred/", protein, "_rsid.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
