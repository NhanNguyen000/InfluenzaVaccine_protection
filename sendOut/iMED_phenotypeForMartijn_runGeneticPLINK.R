rm(list = ls())
library(tidyverse)

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}
# genetic data ----------------------------------
library(data.table)
iMED_vcf <- fread("/vol/projects/CIIM/Influenza/iMED/genotype/genotypes_combined_new/qtl_mapping/iMED_vcf.vcf",
                  nrow = 100)
subjectIDs <- names(iMED_vcf)[grep("I-", names(iMED_vcf))]

# metadata for the selected subjects -------------------------
load("/vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/cohorts_dat.RData")

iMED_pheno <- cohorts$HAI_all %>% filter(probandID %in% subjectIDs) %>%
  mutate(probandID = factor(probandID, levels = subjectIDs)) %>% 
  arrange(probandID) %>%
  mutate_at(vars(contains("reclassify")), ~convert_protectees(.x))

all.equal(as.vector(subjectIDs), as.vector(iMED_pheno$probandID)) # TRUE -> the same order as iMED_vcf

# write the data files for GWAS analysis -----------------------------
iMED_pheno_H1N1 <- iMED_pheno %>% select(probandID, H1N1_reclassify) %>%
  rename("#IID" = "probandID") %>%
  mutate(H1N1_reclassify = ifelse(H1N1_reclassify == "LL", 0, 1))
#write.table(iMED_pheno_H1N1, file = "sendOut/iMED_pheno_H1N1.txt", quote = FALSE, row.names = FALSE)

iMED_pheno_H3N2 <- iMED_pheno %>% select(probandID, H3N2_reclassify) %>%
  rename("#IID" = "probandID") %>%
  mutate(H3N2_reclassify = ifelse(H3N2_reclassify == "LL", 0, 1))
#write.table(iMED_pheno_H3N2, file = "sendOut/iMED_pheno_H3N2.txt", quote = FALSE, row.names = FALSE)

iMED_pheno_B <- iMED_pheno %>% select(probandID, B_reclassify) %>%
  rename("#IID" = "probandID") %>%
  mutate(B_reclassify = ifelse(B_reclassify == "LL", 0, 1))
#write.table(iMED_pheno_B, file = "sendOut/iMED_pheno_B.txt", quote = FALSE, row.names = FALSE)

# write the data files for QTL analysis -----------------------------
iMED_H1N1_T4log2 <- iMED_pheno %>% select(probandID, H1N1_T4_log2) %>%
  rename("#IID" = "probandID") 
write.table(iMED_H1N1_T4log2, file = "sendOut/iMED_H1N1_T4log2.txt", quote = FALSE, row.names = FALSE)

iMED_H3N2_T4log2 <- iMED_pheno %>% select(probandID, H3N2_T4_log2) %>%
  rename("#IID" = "probandID") 
write.table(iMED_H3N2_T4log2, file = "sendOut/iMED_H3N2_T4log2.txt", quote = FALSE, row.names = FALSE)

iMED_B_T4log2 <- iMED_pheno %>% select(probandID, B_T4_log2) %>%
  rename("#IID" = "probandID") 
write.table(iMED_B_T4log2, file = "sendOut/iMED_B_T4log2.txt", quote = FALSE, row.names = FALSE)
