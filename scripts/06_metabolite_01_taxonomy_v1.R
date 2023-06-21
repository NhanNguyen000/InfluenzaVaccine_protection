rm(list = ls())
library(tidyverse)
library(openxlsx)
library(fgsea)

get.tstat <- function(resLimma) {
  # Aim: extract the t-statistic value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_tstat <- resLimma$t %>% as.data.frame %>% 
    dplyr::select(matches("reclassify"))
  
  return(res_tstat)
}

# load the data ----------------
load("cohorts_dat.RData")

# metabolite taxonomy with HMDB ids
mebo_taxo <- read.delim("/vol/projects/CIIM/Influenza/iMED/metabolic/db/hmdb/metabolite_HMDB_taxonomy.csv",
                quote = "", header = TRUE)

# link the HMDB ids to Formulas
mebo_taxo_fomula <- mebo_taxo %>% rownames_to_column("CompoundID") %>% 
  inner_join(read.xlsx('/vol/projects/CIIM/Influenza/iMED/metabolic/raw_data/tables/DATA_CURATED_reformatted.xlsx',
                       sheet = 'annotation') %>% fill(ionIdx, .direction = "down")) %>% 
  filter(Formula %in% colnames(mebo_Dat$iMED_2014)) %>%
  select(Formula, super_class, class, sub_class) %>% distinct()

length(unique(mebo_taxo_fomula$Formula)) # 499 formulas out of 508 formulas

mebo_taxo_fomula %>% count(super_class)
mebo_taxo_fomula %>% count(class)

# make the class 
mebo_taxo_class <- mebo_taxo_fomula %>% select(Formula, class) %>% distinct()
mebo_classSet <- split(mebo_taxo_class$Formula, mebo_taxo_class$class)

## tstat for GSEA ---------------------------------------------------
load("resMebo_4reclass_2group.RData")

# tstat from previous lima analysis
tstat_dat <- resMebo_4reclass_2group %>%
  lapply(function(x) x %>% 
           lapply(function(y) get.tstat(y) %>% as.data.frame %>% rownames_to_column("Formula")))

# rank data set
ranks_dat <- tstat_dat %>%
  lapply(function(x) x %>%
           lapply(function(y) y %>% 
                    dplyr::select(Formula, reclassifyprotectee) %>% deframe(.)))

ranks_dat2 <- ranks_dat %>%
  lapply(function(x) x %>% map(discard, is.na) %>% compact())

# GSEA
fgsea_res <- ranks_dat2 %>%
  lapply(function(x) x %>% 
           lapply(function(y) fgsea(pathways = mebo_classSet, stats=y, nPerm = 5000) %>%
                    as_tibble() %>% arrange(desc(NES))))

fgsea_padj <- fgsea_res %>%
  lapply(function(x) x %>% 
           lapply(function(y) y %>% filter(padj < 0.05)))

fgsea_padj_names <- fgsea_padj %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% dplyr::select(pathway)) %>% unlist()) %>%
  unlist() %>% unique() # 8 pathways are unique

# plot the sig. pathway (padj < 0.05) -----------------
gsea_total <- fgsea_res %>% 
  lapply(function(x) x %>% 
           imap(~mutate(.x, strain = .y)) %>% 
           purrr::reduce(full_join)) %>%
  imap(~mutate(.x, season = .y)) %>%
  purrr::reduce(full_join)

gsea_plot <- gsea_total %>% filter(padj < 0.05) %>%
  rename("Number of Metabolites" = "size") %>%
  mutate(strain = gsub("_reclassify", "", strain)) %>%
  separate(season, c("cohort", "season"))

gsea_plot %>%
  ggplot(aes(x = -log10(padj), y = pathway, size = `Number of Metabolites`, color = strain)) + 
  geom_point() + 
  facet_grid(season ~., scales = "free_y", space = "free") + 
  xlim(c(0, 3)) + theme_bw()

# plot the heatmap --------------------------------
seleted_pathways <- gsea_total %>% filter(padj < 0.05) %>% 
  select(pathway) %>% distinct() %>% unlist()

plotDat <- gsea_total %>% 
  filter(pathway %in% seleted_pathways) %>%
  separate(col = "season", into = c("cohort", "year")) %>%
  mutate(strain = gsub("_reclassify", "", strain)) %>%
  mutate(group = paste0(year, "_", strain))

plotDat %>%
  ggplot(aes(x = group, y = pathway, fill = NES)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(padj < 0.05, "*", NA))) +
  scale_fill_gradientn(limits = c(-3, 3), colors = c("blue", "white", "red"), na.value = "grey") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  ggtitle("Metabolite taxonomy class (padj < 0.05)")
