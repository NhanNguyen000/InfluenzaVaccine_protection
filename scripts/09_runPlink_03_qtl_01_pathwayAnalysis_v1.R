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

qtl <- list()
for (strain in strains) {
  qtl[[strain]] <- read.table(
    file = paste0("/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/magma_results/QTL/", 
                  strain, "/magma.gsa.out"), header = TRUE)
}


qtl_genes <- list()
for (strain in strains) {
  qtl_genes[[strain]] <- read.table(
    file = paste0("/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/magma_results/QTL/", 
                  strain, "/magma.genes.out"), header = TRUE)
}

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


# GWAS pathways =================================================
## load data --------------------------------------
qtl_all <- qtl %>% 
  lapply(function(x) x %>% 
           # mutate(Padj = p.adjust(P, method = "fdr")) %>% # no sig. in Padj
           separate(FULL_NAME, c("group", "FULL_NAME"), ":")) %>% 
  bind_rows(.id = "strain") 

qtl_sigPathway <- qtl_all %>% filter(P < 0.05) %>%
  add_count(FULL_NAME) %>% filter(n == 3)

## select interested pathways --------------------------------------
# plot
qtl_sigPathway  %>% filter(P < 0.05) %>% 
  ggplot(aes(x = -log10(P), y = FULL_NAME, color = strain, size = NGENES)) +
  geom_point() +
  theme_bw() + theme(axis.text = element_text(size = 11)) +
  xlim(0, 3)

#vocano plot
qtl_all %>% filter(P < 0.05) %>% 
  ggplot(aes(x = BETA_STD, y = -log10(P))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlim(0, 0.05) + theme_bw()

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
unique(loci_info_all$trait) # find some imune-related traits

selected_traits <- c("Basal metabolic rate",
                     "1c+mDC:%32+; CD1c+ mDC subset (CD32+)") # the most related trait in the result so far

a <- loci %>% rename("beta_cal" = beta) %>%
  right_join(loci_info_all %>% filter(trait %in% selected_traits)) 
