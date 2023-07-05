rm(list = ls())
library(tidyverse)

# load data =================================================
load("cohorts_dat.RData")
sigSNPs <- read.table("sendOut/qtl_sigSNPs_H3N2.txt", header = TRUE)
sigSNPs_genotype_raw <- read.table("processedDat/qtl_sigSNPs_H3N2_genotype.txt", header = TRUE)

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


#sigSNPs_genotype_H3N2_T4 <- sigSNPs_genotype %>% left_join(cohorts$HAI_all)
sigSNPs_genotype_H3N2_T4 <- sigSNPs_genotype_v2 %>% 
  left_join(cohorts$HAI_all) %>% mutate(H3N2_T4_log2 = log2(H3N2_T4))

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
