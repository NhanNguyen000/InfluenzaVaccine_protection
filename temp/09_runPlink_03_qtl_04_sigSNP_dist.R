rm(list = ls())
library(tidyverse)

# load data =================================================
load("cohorts_dat.RData")
sigSNPs <- read.table("sendOut/qtl_sigSNPs_H3N2.txt", header = TRUE)
sigSNPs_genotype_raw <- read.table("processedDat/qtl_sigSNPs_H3N2_genotype.txt", header = TRUE)

# convert genotype keep 2 hetero-types -----------------------
sigSNPs_genotype <- sigSNPs_genotype_raw %>% 
  mutate(across(c(10:168), function(x) substring(x, 1, 3))) %>%
  mutate(across(c(10:168), 
                function(x) paste0(ifelse(substring(x, 1, 1) == "0", REF, 
                                          ifelse(substring(x, 1, 1) == "1", ALT, 
                                                 substring(x, 1, 1))),
                                   ifelse(substring(x, 3, 3) == "0", REF, 
                                          ifelse(substring(x, 3, 3) == "1", ALT, 
                                                 substring(x, 1, 1)))))) %>%
  separate(col = ID, into = c("CHR:POS", "rsid"), sep = "%") %>% 
  select(rsid, matches("I.")) %>% select(-2, -3) %>%
  column_to_rownames("rsid") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("probandID") %>% mutate(probandID = gsub("[.]", "-", probandID)) 

# convert genotype merge 2 hetero-types into one -----------------------
sigSNPs_genotype_v2 <- sigSNPs_genotype_raw %>% 
  mutate(across(c(10:168), function(x) paste0(substring(x, 1, 1), substring(x, 3, 3)))) %>%
  mutate(across(c(10:168), 
                function(x) ifelse(x == "00", paste0(REF, REF),
                                   ifelse(x == "11", paste0(ALT, ALT),
                                          ifelse(x == "10" |x == "01", paste0(REF, ALT),
                                                 ifelse(x == "..", "..", NA)))))) %>%
  separate(col = ID, into = c("CHR:POS", "rsid"), sep = "%") %>% 
  select(rsid, matches("I.")) %>% select(-2, -3) %>%
  column_to_rownames("rsid") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("probandID") %>% mutate(probandID = gsub("[.]", "-", probandID)) 

# genotype + abTiter data, and plot  -----------------------
#sigSNPs_genotype_H3N2_T4 <- sigSNPs_genotype %>% left_join(cohorts$HAI_all) # use raw abTiter T4
sigSNPs_genotype_H3N2_T4 <- sigSNPs_genotype_v2 %>% 
  left_join(cohorts$HAI_all) %>% mutate(H3N2_T4_log2 = log2(H3N2_T4)) # use log2(abTiter T4)

plotList <- list()
for (rsid in sigSNPs$rsid) {
  plotList[[rsid]] <- sigSNPs_genotype_H3N2_T4 %>% select(c("probandID", rsid, "H3N2_T4", "H3N2_T4_log2")) %>% 
    #ggplot(aes_string(x = rsid, y = "H3N2_T4")) + 
    ggplot(aes_string(x = rsid, y = "H3N2_T4_log2")) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + 
    theme_bw()
}

pdf(file = "sendOut/sigSNPs_dist_H3N2_T4.pdf", width = 18, height = 12)
gridExtra::grid.arrange(grobs = plotList, ncol = 6)
dev.off()

# use normalized(abTiter T4) - Martijn did ----------------
normalized_abTiter <- read.table("/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/qtl_t4_abtiters/phenotypes.csv") %>% 
  t() %>% as.data.frame %>% rownames_to_column("probandID") %>%
  mutate(probandID = gsub("X", "", probandID)) %>%
  mutate(probandID = ifelse(nchar(probandID) == 1, 
                            paste0("I-000", probandID), 
                            ifelse(nchar(probandID) == 2, paste0("I-00", probandID), 
                                   ifelse(nchar(probandID) == 3, 
                                          paste0("I-0", probandID), paste0("I-", probandID)))))

sigSNPs_genotype_H3N2_T4 <- sigSNPs_genotype_v2 %>% 
  left_join(normalized_abTiter)

plotList <- list()
for (rsid in sigSNPs$rsid) {
  plotList[[rsid]] <- sigSNPs_genotype_H3N2_T4 %>% select(c("probandID", rsid, "H3N2_T4")) %>% 
    ggplot(aes_string(x = rsid, y = "H3N2_T4")) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + 
    theme_bw()
}

pdf(file = "sendOut/sigSNPs_normalizedDist_H3N2_T4.pdf", width = 18, height = 12)
gridExtra::grid.arrange(grobs = plotList, ncol = 6)
dev.off()
