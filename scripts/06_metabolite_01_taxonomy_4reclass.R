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

get.plotDat_clusterRow <- function(plotDat, colName, valColumn) {
  # Aim: make the heatmap data with dendrogram clustering order, so can have clustered heatmap with ggplot
  # this function comvert the plot data (long format) to wide format 
  # with column name from the given "colName" column, value from given "pathway" colum, and rowname from fixed "pathway" column
  
  # outcome: plotDat_order - plot data long format with order in pathway, so it will show clustering with ggplot heatmap
  
  plotDat_wide <- plotDat %>% select(c("pathway", colName, valColumn)) %>% 
    pivot_wider(names_from = colName, values_from = valColumn) %>% 
    column_to_rownames("pathway")
  
  plotDat_dendrogram <- as.dendrogram(hclust(d = dist(plotDat_wide)))
  plotDat_denOrder <- order.dendrogram(plotDat_dendrogram)
  
  plotDat_order <- plotDat %>%
    mutate(pathway = factor(pathway, levels = pathway[plotDat_denOrder], ordered = TRUE))
  
  return(plotDat_order)
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

## tstat for GSEA for 4 reclass ---------------------------------------------------
load("resMebo_4reclass.RData")

# tstat from previous lima analysis
tstat_dat <- resMebo_4reclass %>%
  lapply(function(x) x %>% 
           lapply(function(y) get.tstat(y) %>% 
                    as.data.frame %>% rownames_to_column("Formula") %>%
                    pivot_longer(cols = starts_with("reclassify")) %>%
                    split(.$name) %>% map(~.x %>% select(-name))))

# rank data set
ranks_dat <- tstat_dat %>%
  lapply(function(x) x %>%
           lapply(function(y) y %>% 
                    lapply(function(z) z %>% 
                             dplyr::select(Formula, value) %>% deframe(.))))
# stop here --------------------------
ranks_dat2 <- ranks_dat %>%
  lapply(function(x) x %>% lapply(function(y) y %>% map(discard, is.na) %>% compact()))

# GSEA
fgsea_res <- ranks_dat2 %>%
  lapply(function(x) x %>% 
           lapply(function(y) y %>%
                    lapply(function(z) fgsea(pathways = mebo_classSet, stats=z, nPerm = 5000) %>%
                    as_tibble() %>% arrange(desc(NES)))))

fgsea_padj <- fgsea_res %>%
  lapply(function(x) x %>% 
           lapply(function(y) y %>% 
                    lapply(function(z) z %>% filter(padj < 0.05))))

fgsea_padj_names <- fgsea_padj %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% 
                    lapply(function(z) z %>% dplyr::select(pathway)) %>% unlist())) %>%
  unlist() %>% unique() # 8 pathways are unique

# plot the sig. pathway (padj < 0.05) for 2015 -----------------
gsea_total <- fgsea_res %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% 
                    imap(~mutate(.x, compare = .y)) %>% purrr::reduce(full_join)) %>%
           imap(~mutate(.x, strain = .y)) %>% purrr::reduce(full_join)) %>%
  imap(~mutate(.x, season = .y)) %>% purrr::reduce(full_join) %>%
  mutate(compare = gsub("reclassify", "LLvs", compare),
         strain = gsub("_reclassify", "", strain)) %>% 
  separate(col = "season", into = c("cohort", "year"))

gsea_plot <- gsea_total %>% filter(padj < 0.05) %>%
  rename("Number of Metabolites" = "size") 

gsea_plot %>% filter(year == "2015") %>%
  ggplot(aes(x = -log10(padj), y = pathway, 
             size = `Number of Metabolites`, color = compare)) + 
  geom_point() + 
  facet_grid(strain ~., scales = "free_y", space = "free") + 
  xlim(c(0, 8)) + 
  theme_bw()

# plot the heatmap --------------------------------
seleted_pathways <- gsea_total %>% filter(padj < 0.05) %>% 
  select(pathway) %>% distinct() %>% unlist()

plotDat_2015 <- gsea_total %>% 
  filter(pathway %in% seleted_pathways, year == "2015") %>%
  mutate(compare = paste0(strain, "_", compare),
         compare = factor(compare, 
                          levels = c("H1N1_LLvsLH", "H1N1_LLvsHL", "H1N1_LLvsHH",
                                     "H3N2_LLvsLH", "H3N2_LLvsHL", "H3N2_LLvsHH",
                                     "B_LLvsLH", "B_LLvsHL", "B_LLvsHH")))

plotDat_2015 %>%
  ggplot(aes(x = compare, y = pathway, fill = NES)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(padj < 0.05, "*", NA))) +
  scale_fill_gradientn(limits = c(-3, 3), colors = c("blue", "white", "red"), na.value = "grey") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  ggtitle("Metabolite taxonomy class (padj < 0.05)")

plotDat_order <- get.plotDat_clusterRow(plotDat_2015, 
                                        colName = "compare", 
                                        valColumn = "NES")
plotDat_order %>%
  ggplot(aes(x = compare, y = pathway, fill = NES)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(padj < 0.05, "*", NA))) +
  scale_fill_gradientn(limits = c(-3, 3), colors = c("blue", "white", "red"), na.value = "grey") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  ggtitle("Metabolite taxonomy class (padj < 0.05)")
