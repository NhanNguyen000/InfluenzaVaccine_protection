#load result from Martijn analysis
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

# run this code in linux =================================================
cd /vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/new/
awk '$13<10e-5' h1n1.H1N1_reclassify.glm.logistic.hybrid > /vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/gwas/h1n1_suggestiveSig.txt
awk '$13<10e-5' b.B_reclassify.glm.logistic.hybrid > /vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/gwas/b_suggestiveSig.txt

# load data =================================================
gwas <- list()
gwas$H1N1 <- read.table("processedDat/gwas/h1n1_suggestiveSig.txt")
gwas$B <- read.table("processedDat/gwas/b_suggestiveSig.txt")

gwas_colnames <- c("CHROM", "POS", "ID", "REF", "ALT",	"A1", "FIRTH", # based on the "h1n1.H1N1_reclassify.glm.logistic.hybrid" file
                   "TEST",	"OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT",	"P",	"ERRCODE")

names(gwas$H1N1) <- gwas_colnames
names(gwas$B) <- gwas_colnames

## info for sig. SNPs --------------------------------------
loci <- gwas %>%  bind_rows(.id = "strain")
loci %>% group_by(strain) %>% count()

loci_v2 <- loci %>% filter(P < 1e-6)
loci_v2 %>% group_by(strain) %>% count()

loci_v2 <- loci %>% filter(P < 1e-5)
loci_v2 %>% group_by(strain) %>% count()

# load the SNP ids --------------------------------------
rsid_info <- read.table("/vol/projects/CIIM/resources/snpDB/Hg19_snploc.txt")

# H1N1
rs4871494	8	120697077
rs12551619	9	14205169	G	A
# B
