rm(list = ls())
library(tidyverse)
library(phenoscanner)

phenoscanner_longList <- function(loci_rsids, catalogue) {
  snps <- c()
  results <- c()
  loci_rsids_chunks <- split(loci_rsids, ceiling(seq_along(loci_rsids)/100))
  
  for (chunk_name in names(loci_rsids_chunks)) {
    snps <- rbind(snps, phenoscanner(loci_rsids_chunks[[chunk_name]], catalogue = catalogue)$snps)
    results <- rbind(results, phenoscanner(loci_rsids_chunks[[chunk_name]], catalogue = catalogue)$results)
  }
  return(list("snps" = snps, "results" = results))
}

# load data =================================================
times <- c("T1", "T4")

sig_SNIPs <- list()
for (time in times) {
  sig_SNIPs[[time]] <- read.table(
    file = paste0("processedDat/qtl/outcome_", time, ".csv"), header = TRUE) %>% 
    as.data.frame %>% filter(p.value < 10e-5)
}

## info for sig. SNPs --------------------------------------
loci <- sig_SNIPs %>%  bind_rows(.id = "time")
loci %>% group_by(gene) %>% count()

loci_v2 <- loci %>% filter(p.value < 10e-6)
loci_v2 %>% group_by(gene) %>% count()

# loci_v2 <- loci %>% filter(p.value < 10e-7)
# loci_v2 %>% group_by(gene) %>% count()
# 
# loci_v2 <- loci %>% filter(p.value < 10e-8)
# loci_v2 %>% group_by(gene) %>% count()
# 
# loci_v2 <- loci %>% group_by(SNP) %>% add_count() %>%
#   filter(n > 1)

snps <- loci_v2 %>% separate(SNP, into = c("location", "SNP"), sep = "[%]")

snps_info <- phenoscanner_longList(loci_rsids = snps$SNP, catalogue = "eQTL")

# save data =================================================
save(snps_info, file = "snps.RData")

load("snps.RData")
write.table(unique(snps_info$results$exp_gene), file = "processedDat/qtl/snp_associatedGenes.txt", 
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

# result form Martijn =================================================
path <- "/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/qtl_t4_abtiters"
sig_SNIPs_martijn <- list()
for (time in c("h1n1", "h3n2", "b")) {
  sig_SNIPs[[time]] <- read.table(
    file = paste0(path, strain, "_suggestivesig.txt"), header = TRUE) %>% 
    as.data.frame %>% filter(p.value < 10e-5)
}
a <- read.table("/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/qtl_t4_abtiters/out.csv")
b <- read.ftable("processedDat/qtl/outcome_T4.csv")