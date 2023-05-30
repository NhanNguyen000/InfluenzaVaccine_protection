rm(list = ls())

library(tidyverse)
library(fgsea)
library('org.Hs.eg.db')

get.tstat <- function(resLimma) {
  # Aim: extract the t-statistic value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_tstat <- resLimma$t %>% as.data.frame %>% 
    dplyr::select(matches("reclassify"))
  
  return(res_tstat)
}


# load hallmark gene sets from MSigDB --------------------------
hallmarkSets <- gmtPathways("reference/h.all.v2022.1.Hs.entrez.gmt.txt")
hallmarkSets %>% head() %>% lapply(head)

# load data =======================================================================


## convert protein to Entrez gene ID----------------------------------
load("cohorts_dat.RData")
Entrez_IDs <- AnnotationDbi::select(org.Hs.eg.db, 
                     keys = colnames(protein_Dat$iMED_2014),
                     keytype = 'SYMBOL', columns = c('ENTREZID', 'SYMBOL'))

## tstat for GSEA ---------------------------------------------------
load("resPro_4reclass_2group.RData")

# tstat from previous lima analysis
tstat_dat <- resPro_4reclass_2group %>%
  lapply(function(x) x %>% 
           lapply(function(y) get.tstat(y) %>% as.data.frame %>%
                    rownames_to_column("SYMBOL") %>% full_join(Entrez_IDs)))

# rank data set
ranks_dat <- tstat_dat %>%
  lapply(function(x) x %>%
           lapply(function(y) y %>% 
                    dplyr::select(ENTREZID, reclassifyprotectee) %>% 
                    distinct() %>% group_by(ENTREZID) %>% summarise(tstat = mean(reclassifyprotectee)) %>%
                    deframe(.)))
ranks_dat2 <- ranks_dat %>%
  lapply(function(x) x %>% map(discard, is.na) %>% compact())

# GSEA
fgsea_res <- ranks_dat2 %>%
  lapply(function(x) x %>% 
           lapply(function(y) fgsea(pathways=hallmarkSets, stats=y, nPerm = 5000) %>%
                    as_tibble() %>% arrange(desc(NES))))

fgsea_padj <- fgsea_res %>%
  lapply(function(x) x %>% 
           lapply(function(y) y %>% filter(padj < 0.05)))

fgsea_padj_names <- fgsea_padj %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% dplyr::select(pathway)) %>% unlist()) %>%
  unlist() %>% unique() # 10 pathways are unique

# plot the sig. pathway (padj < 0.05) -----------------
gsea_total <- fgsea_res %>% 
  lapply(function(x) x %>% 
           imap(~mutate(.x, strain = .y)) %>% 
           purrr::reduce(full_join)) %>%
  imap(~mutate(.x, season = .y)) %>%
  purrr::reduce(full_join)

gsea_plot <- gsea_total %>% filter(padj < 0.05) %>%
  rename("size" = "Number of protein") %>%
  mutate(strain = gsub("_reclassify", "", strain)) %>%
  separate(season, c("cohort", "season"))

gsea_plot %>%
  ggplot(aes(x = -log10(padj), y = pathway, size = `Number of protein`, color = strain)) + 
  geom_point() + 
  facet_grid(season ~., scales = "free_y", space = "free") + 
  xlim(c(0, 3)) + theme_bw()

