rm(list = ls())

library(tidyverse)
library(ggpubr)
library(org.Hs.eg.db)
# load data =======================================================================
load("processedDat/cohorts_dat.RData")
load("processedDat/omicsPred_output.RData")

pred_dat <- pred_Tab %>% 
  t() %>% as.data.frame %>% rownames_to_column("OMICSPRED.ID") %>% 
  full_join(selected_models %>% dplyr::select(OMICSPRED.ID, model, type, Gene, Ensembl.ID))

# protein -------------------------------------------------------------------------
## prediction and measurement data ---------------------------------------------------
pred_proDat_temp <- pred_dat %>% 
  filter(model %in% c("Olink", "Somalogic")) %>%
  dplyr::select(-OMICSPRED.ID, -type, -Ensembl.ID)


measured_proDat <- cohorts$donorSample_all %>% # protein measures from coresspoding subjects
  filter(time == "d0") %>% dplyr::select(probandID, name) %>%
  filter(probandID %in% colnames(pred_proDat_temp)) %>%
  left_join(protein_Dat$iMED_2015 %>% rownames_to_column("name")) %>%
 # arrange(match(probandID, rownames(pred_proDat_temp))) %>% 
  column_to_rownames("probandID") %>% dplyr::select(-name) %>% as.data.frame()

pred_proDat <- pred_proDat_temp %>% # protein prediction from coresspoding subjects
  filter(Gene %in% colnames(measured_proDat)) # %>%  column_to_rownames("Gene") %>% t() %>% as.data.frame()


## scatter plot -----------------------------------------------------------------

#model_name <- "Olink"
model_name <- "Somalogic"
protein <- "CCL7"

dat_temp <- measured_proDat %>% 
  dplyr::select(protein) %>%  tibble::rownames_to_column("name") %>%
  full_join(
    as.data.frame(
      pred_proDat %>% 
        filter(Gene %in% protein, model %in% model_name) %>% 
        mutate(varName =  paste0(Gene, "_", model)) %>% 
        dplyr::select(-Gene, -model) %>% column_to_rownames("varName") %>% 
        t()) %>% 
      tibble::rownames_to_column(var = "name")) %>% 
  column_to_rownames("name")

## save the plot 
png(paste0("output/omicPred_measureVsPredict_protein_" , protein, ".png"))

dat_temp %>% as.data.frame() %>% 
  ggplot(aes_string(x = protein, y = paste0(protein, "_", model_name))) + 
  geom_point(size = 5, alpha = 0.5) +
  geom_smooth(method = "lm", se=FALSE) + stat_cor(size = 8) + 
  xlab(paste0(protein, "_measure")) +
  theme_classic()+
  theme(text = element_text(size = 20))

dev.off()

# RNAseq data -----------------------------------------------------------------------

# RNAseq data of 20 donors in season 2015
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1")

## prediction and measurement data ---------------------------------------------
# RNAseq prediction based on the genetic data of 159 donors in season 2015
matched_IDs <- intersect(iMED_transcrip_T1$patientID, names(pred_dat)) # 13 donors

pred_rnaDat <- pred_dat %>% 
  filter(model == "RNAseq") %>% # RNAseq prediction from 13 out of 20 participants have RNAseq data
  dplyr::select(c("OMICSPRED.ID", "model", "type", "Gene", "Ensembl.ID", matched_IDs)) %>% # Note: some row have Ensemble ID but no gene name
  left_join(mapIds(org.Hs.eg.db, keys = pred_dat$Ensembl.ID, # add missing gene name using org.Hs.eg.db with Ensemble ID
                   column = "SYMBOL", keytype =  "ENSEMBL") %>% 
              unlist() %>% as.data.frame() %>% 
              dplyr::rename("Gene2" = ".") %>% rownames_to_column("Ensembl.ID")) %>%
  mutate(Gene = ifelse(Gene == "", Gene2, Gene)) %>% 
  dplyr::select(-OMICSPRED.ID, -model, -type, -Ensembl.ID, -Gene2) %>% 
  drop_na() %>% # remove genes have no gene name and could not convert Ensembl.ID to gene name
  column_to_rownames("Gene") %>% t() %>% as.data.frame()


measured_rnaDat <- transcriptomeList[[1]] %>% # RNAseq measures from 13 coresspoding subjects
  t() %>% as.data.frame %>% rownames_to_column("SampleName") %>% 
  right_join(iMED_transcrip_T1 %>% dplyr::select(SampleName, patientID)) %>%
  filter(patientID %in% rownames(pred_rnaDat)) %>%
  arrange(match(patientID, rownames(pred_rnaDat))) %>% 
  column_to_rownames("patientID") %>% dplyr::select(-SampleName)

identical(rownames(measured_rnaDat), rownames(pred_rnaDat)) # TRUE = same order of subjects

## scatter plot -----------------------------------------------------------------
gene <- "CD83" # no data in the selcted models

gene <- "TMEM51"
gene <- "IFT52"
gene <- "ARSA"
gene <- "RIPOR1"
gene <- "CRYL1"
gene <- "TMEM204"
gene <- "ACVR2B"
gene <- "RNMT"

dat_temp <- cbind("predicted" = pred_rnaDat[, gene],
                  "measured" = measured_rnaDat[, gene])

## save the plot 
png(paste0("output/omicPred_measureVsPredict_RNAseq_" , gene, ".png"))

dat_temp %>% as.data.frame() %>% 
  ggplot(aes(x = measured, y = predicted)) + 
  geom_point(size = 5, alpha = 0.5) +
  geom_smooth(method = "lm", se=FALSE) + stat_cor(size = 8) + 
  theme_classic()+
  theme(text = element_text(size = 24)) + ggtitle(gene)

dev.off()

## calculate the correlation and p value ----------------------------------------
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


## sum up the cor(predicted, measured values) for protein Olink and RNAseq data ------------------
pred_inf <- pred_rnaInf #%>% full_join(pred_proInf)

pred_inf %>% 
  ggplot(aes(x = cor, y = -log10(pval), col = model)) + 
  geom_jitter(size = 1, width = 0.15, height = 0.15) +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  theme_classic()
