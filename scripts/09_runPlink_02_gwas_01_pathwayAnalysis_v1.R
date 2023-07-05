rm(list = ls())
library(tidyverse)

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