rm(list = ls())

library(tidyverse)
library(fgsea)
library('org.Hs.eg.db')

get.DAs <- function(resLimma) {
  # Aim: extract the DAPs/DAMs with p.value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_DA <- resLimma$p.value %>% as.data.frame %>% 
    dplyr::select(matches("reclassify")) %>% filter(. <0.05)
  
  return(res_DA)
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

Entrez_IDs_pathway <- Entrez_IDs
for(pathway in names(hallmarkSets)) {
  Entrez_IDs_pathway[[pathway]] <- Entrez_IDs_pathway$ENTREZID %in% hallmarkSets$HALLMARK_TNFA_SIGNALING_VIA_NFKB
}



## tstat for GSEA ---------------------------------------------------
load("resPro_4reclass_2group.RData")

# tstat from previous lima analysis
DAs <- resPro_4reclass_2group %>% 
  lapply(function(x) x %>%
           lapply(function(y) get.DAs(y) %>% 
                    as.data.frame %>% rownames_to_column("SYMBOL") %>% 
                    left_join(Entrez_IDs_pathway)))

DAs_pathway <- DAs %>% 
  lapply(function(x) x %>% 
           imap(~mutate(.x, strain = .y)) %>% 
           purrr::reduce(full_join)) %>%
  imap(~mutate(.x, season = .y)) %>%
  purrr::reduce(full_join) %>% 
  slice(which(SYMBOL %in% selected_DAs))

a <- DAs_pathway %>% dplyr::select(where(~ is.logical(.) && any(.)), where(negate(is.logical)))
colSums(DAs_pathway[, c(4:53)])
rowSums(DAs_pathway[, 4:53])
DAs_pathway$SYMBOL[which(rowSums(DAs_pathway[, 4:53]) !=0)]

# Check venn diagram ---------------------------------------------
res.venn <- DAs %>% lapply(function(x) x %>% lapply(function(y) rownames(y)))
selected_DAs <- unique(as.vector(unlist(res.venn)[which(duplicated(unlist(res.venn)))])) # DAs in at least twice acorss strain and season


