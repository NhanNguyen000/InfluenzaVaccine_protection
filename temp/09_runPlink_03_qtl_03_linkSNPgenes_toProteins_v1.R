rm(list = ls())
library(tidyverse)

# load data =================================================
strains <- c("H1N1", "H3N2", "B")

sig_SNIPs <- list()
for (strain in strains) {
  sig_SNIPs[[strain]] <- read.table(
    file = paste0("/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/qtl_t4_abtiters/", 
                  tolower(strain), "_suggestivesig.txt"), header = TRUE)
}

loci <- sig_SNIPs %>% 
  bind_rows(.id = "strain") %>% filter(SNP != ".") %>%
  dplyr::rename('rsid' = SNP, 'pval'= p)

write.table(loci %>% filter(pval < 5e-8), file = "sendOut/qtl_sigSNPs_H3N2.txt", 
            row.names = FALSE, quote = FALSE)

write.table((loci %>% filter(pval < 5e-8))$rsid,
            file = "sendOut/qtl_sigSNPs_H3N2_onlyRsid.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
## info for sig. SNPs --------------------------------------
load("loci_info_qtl.RData")

identical(loci_info$gwas$snps,loci_info$eqtl$snps) # TRUE, use one of them
loci_strains <- loci %>% full_join(loci_info$gwas$snps)

# check if gene related to SNPs also in the DA proteins =================================================
load("selected_DAPs.RData")
intersect(loci_strains$hgnc, selected_DAs) # no overlap gene-related to sig. SNPs and selected DA proteins
