# load libraries & functions --------------------------
library(fgsea)
library('org.Hs.eg.db')

# load BTM set from https://github.com/shuzhao-li/BTM  --------------------------
BTM_sets <- gmtPathways("reference/BTM_for_GSEA_20131008.gmt")
BTM_sets %>% head() %>% lapply(head)
