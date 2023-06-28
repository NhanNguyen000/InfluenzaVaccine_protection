rm(list = ls())
library(tidyverse)
library(phenoscanner)

phenoscanner_longList <- function(loci_rsids, catalogue) {
  snps <- c()
  results <- c()
  loci_rsids_chunks <- split(loci_rsids, ceiling(seq_along(loci$rsid)/100))
  
  for (chunk_name in names(loci_rsids_chunks)) {
    snps <- rbind(snps, phenoscanner(loci_rsids_chunks[[chunk_name]], catalogue = catalogue)$snps)
    results <- rbind(results, phenoscanner(loci_rsids_chunks[[chunk_name]], catalogue = catalogue)$results)
  }
  return(list("snps" = snps, "results" = results))
}
# load data =================================================
strains <- c("H1N1", "H3N2", "B")

sig_SNIPs <- list()
for (strain in strains) {
  sig_SNIPs[[strain]] <- read.table(
    file = paste0("/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/qtl_t4_abtiters/", 
                  tolower(strain), "_suggestivesig.txt"), header = TRUE)
}

## info for sig. SNPs --------------------------------------
loci <- sig_SNIPs %>% 
  bind_rows(.id = "strain") %>% filter(SNP != ".") %>%
  dplyr::rename('rsid' = SNP, 'pval'= p)

# loci_info <-list()
# loci_info$gwas <- phenoscanner_longList(loci_rsids = loci$rsid, catalogue = "GWAS")
# loci_info$eqtl <- phenoscanner_longList(loci_rsids = loci$rsid, catalogue = "eQTL")
# save(loci_info, file = "loci_info_qtl.RData")


# SNPs =================================================
load("loci_info_qtl.RData")

identical(loci_info$gwas$snps,loci_info$eqtl$snps) # TRUE, use one of them

loci_info_all <- loci_info$gwas$snps %>% 
  full_join(loci_info$gwas$results) %>%
  full_join(loci_info$eqtl$results) %>%
  group_by(across(c(-exp_gene)))%>%
  summarise_at("exp_gene", function(x) paste(unique(x), collapse = ";"))%>%
  group_by(across(c(-trait)))%>%
  summarise_at("trait", function(x) paste(unique(x), collapse = ";"))

unique(loci_info_all$exp_gene)
unique(loci_info_all$hgnc)
unique(loci_info_all$trait) # find some immune-related traits

selected_traits <- c("Basal metabolic rate",
                     "1c+mDC:%32+; CD1c+ mDC subset (CD32+)") # the most related trait in the result so far

a <- loci %>% rename("beta_cal" = beta) %>%
  right_join(loci_info_all %>% filter(trait %in% selected_traits)) 
