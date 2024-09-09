rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data =======================================================================
load("processedDat/cohorts_dat.RData")

## metadata for all healthy subjects -------------------------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))


## RNAseq data of 20 donors in season 2015 --------------------------------
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1") # time T1 in trancriptome = d0

# prepare data -----------------------------------------------------------------
selected_genes <- c("ORC6", "SHCBP1", "GPT2", "ITFG1", 
                    "VPS35", "MYLK3", "DNAJA2", 
                    #"RP11-329J18.5", "RP11-93O14.3", # these gene are not in the RNAseq dataset
                    "C16orf87", "NETO2")

dat <- transcriptomeList[[1]] %>% 
  t() %>% as.data.frame %>% select(selected_genes) %>%
  rownames_to_column("SampleName") %>% 
  right_join(iMED_transcrip_T1 %>% dplyr::select(SampleName, patientID)) %>%
  left_join(metadata_healthy %>%  # get H3N2 antibody titer data
              select(probandID, H3N2_d28, H3N2_d28_log2),
             by = c("patientID" = "probandID"))

# correlation between gene expression and H3N2 ab titer at day 28 ----------------------------------------
cor_outcome <- cbind(selected_genes, "cor" = NA, "pval" = NA) %>% as.data.frame()

for (gene in cor_outcome$selected_genes) {
  out_temp <- cor.test(dat[, gene], dat$H3N2_d28_log2, 
                       use = "complete.obs")
  cor_outcome$cor[cor_outcome$selected_genes == gene] <- out_temp$estimate
  cor_outcome$pval[cor_outcome$selected_genes == gene] <- out_temp$p.value
}

## plot data ---------------------------------------------------------------------
gene <- "ITFG1"

dat_temp <- cbind("gene_expression" = dat[, gene],
                  "H3N2_ab_d28_log2" = dat$H3N2_d28_log2)

## save the plot 
png(paste0("output/cor_RNAseq_" , gene, "H3N2_abD28.png"))

dat_temp %>% as.data.frame() %>% 
  ggplot(aes(x = gene_expression, y = H3N2_ab_d28_log2)) + 
  geom_point(size = 5, alpha = 0.5) +
  geom_smooth(method = "lm", se=FALSE) + stat_cor(size = 8) + 
  theme_classic()+
  theme(text = element_text(size = 24)) + ggtitle(gene)

dev.off()

