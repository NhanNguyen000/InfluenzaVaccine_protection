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

gwas <- list()
for (strain in strains) {
  gwas[[strain]] <- read.table(
    file = paste0("/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/magma_results/GWAS/", 
                  strain, "/magma.gsa.out"), header = TRUE)
}

gwas_genes <- list()
for (strain in strains) {
  gwas_genes[[strain]] <- read.table(
    file = paste0("/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/magma_results/GWAS/", 
                  strain, "/magma.genes.out"), header = TRUE)
}

sig_SNIPs <- list()
for (strain in strains) {
  sig_SNIPs[[strain]] <- read.table(
    file = paste0("/vol/projects/CIIM/Influenza/iMED/genotype/nhan_ll_vs_other_gwas/new/suggestivesig_", 
                  strain, ".txt"), header = TRUE)
}

## info for sig. SNPs --------------------------------------
loci <- sig_SNIPs %>% 
  bind_rows(.id = "strain") %>% filter(SNP != ".") %>%
  dplyr::rename('rsid' = SNP, 'pval'= p)

# loci_info <-list()
# loci_info$gwas <- phenoscanner_longList(loci_rsids = loci$rsid, catalogue = "GWAS")
# loci_info$eqtl <- phenoscanner_longList(loci_rsids = loci$rsid, catalogue = "eQTL")
# save(loci_info, file = "loci_info_gwas.RData")

# GWAS pathways =================================================
## load data --------------------------------------
gwas_all <- gwas %>% 
  lapply(function(x) x %>% 
          # mutate(Padj = p.adjust(P, method = "fdr")) %>% # no sig. in Padj
           separate(FULL_NAME, c("group", "FULL_NAME"), ":")) %>% 
  bind_rows(.id = "strain") 

a <- gwas_all %>% filter(P < 0.05) %>%
  filter(strain != "H3N2") %>%
  add_count(FULL_NAME) %>% filter(n == 2)

a2 <- gwas_all %>% filter(P < 0.05) %>%
  add_count(FULL_NAME) %>% filter(n == 3)

## select interested pathways --------------------------------------
pathways <- c("go_immunological_synapse_formation", 
              "go_negative_regulation_of_interleukin_5_production", 
              "go_immunoglobulin_secretion",
              "go_inflammasome_complex",
              "go_nlrp3_inflammasome_complex",
              "go_ciliary_neurotrophic_factor_receptor_activity",
              "go_arginine_transmembrane_transporter_activity",
              "go_coenzyme_transmembrane_transporter_activity")

gwas_selected_pathways <- gwas_all %>% filter(FULL_NAME %in% pathways)

# plot
gwas_selected_pathways %>% filter(P < 0.05) %>% 
  ggplot(aes(x = -log10(P), y = FULL_NAME, color = strain, size = NGENES)) +
  geom_point() +
  theme_bw() + theme(axis.text = element_text(size = 11)) +
  xlim(0, 3)

# SNPs =================================================
load("loci_info_gwas.RData")

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

selected_traits <- c( "White blood cell count",
                      "Sum neutrophil eosinophil counts",                                                         
                      "Granulocyte count",                                                                            
                      "Neutrophil count",                                                                             
                      "Sum basophil neutrophil counts",                                                               
                      "Myeloid white cell count",                                                                    
                      "Red cell distribution width",                                                                  
                      "Monocyte count", 
                      "Basal metabolic rate",      
                      "Neutrophil percentage of white cells",                                                         
                      "Lymphocyte percentage of white cells",                                                         
                      "Lymphocyte count",   
                      "Eosinophil percentage of white cells",                                                         
                      "Eosinophil count",                                                                             
                      "Sum eosinophil basophil counts",     
                      "Arachidonic acid 20:4n6",       
                      "IgA deficiency",  
                      "CD8:%SCM; TSCM (CD27+CD28+CD57-CD95+CD127+CD45RA+)",  
                      "Granulocyte percentage of myeloid white cells",                                                
                      "Monocyte percentage of white cells" )

a <- loci %>% 
  right_join(loci_info_all %>% filter(trait %in% selected_traits)) # all traits are come from H3N2, which is not reliable
