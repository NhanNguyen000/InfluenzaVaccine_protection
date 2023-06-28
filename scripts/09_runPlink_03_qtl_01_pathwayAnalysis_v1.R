rm(list = ls())
library(tidyverse)

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

