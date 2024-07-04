rm(list = ls())

library(tidyverse)
library(venn)

get.DAs <- function(resLimma) {
  # Aim: extract the DAPs/DAMs with p.value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_DA <- resLimma$p.value %>% as.data.frame %>% 
    select(matches("reclassify"))
  
  return(res_DA)
}


get.tstat <- function(resLimma) {
  # Aim: extract the t-statistic value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_tstat <- resLimma$t %>% as.data.frame %>% 
    select(matches("reclassify"))
  
  return(res_tstat)
}

get.plotDat_clusterRow <- function(plotDat, colName, varColumn) {
  # Aim: make the heatmap data with dendrogram clustering order, so can have clustered heatmap with ggplot
  # this function comvert the plot data (long format) to wide format 
  # with column name from the given "colName" column, value from given "varName" colum, and rowname from fixed "varName" column
  
  # outcome: plotDat_order - plot data long format with order in varName, so it will show clustering with ggplot heatmap
  
  plotDat_wide <- plotDat %>% select(c("varName", colName, varColumn)) %>% 
    pivot_wider(names_from = colName, values_from = varColumn) %>% 
    column_to_rownames("varName")
  
  plotDat_dendrogram <- as.dendrogram(hclust(d = dist(plotDat_wide)))
  plotDat_denOrder <- order.dendrogram(plotDat_dendrogram)
  
  plotDat_order <- plotDat %>%
    mutate(varName = factor(varName, levels = varName[plotDat_denOrder], ordered = TRUE))
  
  return(plotDat_order)
}

# load data --------------------------------------------------------
load("processedDat/resMebo_4reclass_2group.RData")

# significant proteins / metabolites ------------------------------------------
DAs <- resMebo_4reclass_2group %>% 
  lapply(function(x) x %>% lapply(function(y) get.DAs(y)))

DAs_all <- DAs %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% rownames_to_column("varName")) %>%
           bind_rows(.id = "strain")) %>%
  bind_rows(.id = "season") %>%
  mutate(strain = gsub("_reclassify", "", strain)) %>%
  rename("p.value" = "reclassifyprotectee")

# t-statistic values ---------------------------------------------------------
tstat_dat <- resMebo_4reclass_2group %>%
  lapply(function(x) x %>% lapply(function(y) get.tstat(y)))

tstat_all <- tstat_dat %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% rownames_to_column("varName")) %>%
           bind_rows(.id = "strain")) %>%
  bind_rows(.id = "season") %>%
  mutate(strain = gsub("_reclassify", "", strain))  %>%
  rename("tstat" = "reclassifyprotectee")

tstat_longDat <- tstat_all %>%
  full_join(DAs_all) %>%
  separate(season, sep = "_", into = c("cohort", "season")) %>%
  mutate(group = paste0(season, "_", strain))

# prepare data with metabolites that show consistent trend across strain and season ------------------------------- 
selectedSeasons <- c("2014_H1N1", "2015_H1N1", "2019_H1N1", "2020_H1N1",
                     "2014_B", "2015_B", "2015_H3N2", "2020_Byamagata") 

consistVars <- tstat_longDat %>%
  filter(group %in% selectedSeasons) %>%
  mutate(direction = ifelse(tstat < 0, "down", ifelse(tstat > 0, "up", tstat))) %>%
  group_by(varName) %>% add_count(direction) %>%
  filter(n >6) %>% # The up/down trend appear in 6 out of 8 comparison (75% across strain and season) 
  select(varName) %>% unlist() %>% unique()

DAs_consistVars <- DAs %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% 
                    rownames_to_column("varName") %>% 
                    filter(varName %in% consistVars) %>%
                    mutate(padj = p.adjust(reclassifyprotectee, method = "fdr"))))

DAs_consistVars_all <- DAs_consistVars %>% 
  lapply(function(x) x %>% bind_rows(.id = "strain")) %>%
  bind_rows(.id = "season") %>% 
  rename("p.value" = "reclassifyprotectee") %>%
  mutate(strain = gsub("_reclassify", "", strain)) 

tstat_longDat_consistVars <- tstat_all %>%
  inner_join(DAs_consistVars_all) %>%
  separate(season, sep = "_", into = c("cohort", "season")) %>%
  mutate(group = paste0(season, "_", strain))

## heatmap --------------------------------
sigPval_vars <- tstat_longDat %>%
  filter(group %in% selectedSeasons) %>%
  filter(p.value < 0.05) %>% select(varName) %>% unlist() %>% unique()

selected_vars <- intersect(consistVars, sigPval_vars)
#save(selected_vars, sigPval_vars, file = "selected_consist_sigPval_DAMs.RData") # 30 and 180 formulas

plotDat <- tstat_longDat_consistVars %>% 
  filter(varName %in% selected_vars) %>%
  slice(which(group %in% selectedSeasons ))

unique((plotDat %>% filter(padj < 0.05))$varName)

plotDat_order <- get.plotDat_clusterRow(plotDat, 
                                        colName = "group", 
                                        varColumn = "tstat")

heatmapPlot <- plotDat_order %>%
  ggplot(aes(x = group, y = varName, fill = tstat)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(padj < 0.05, "**", NA)), size = 10) +
  geom_text(aes(label = ifelse(p.value < 0.05, "*", NA)), size = 10) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1), text = element_text(size = 24))

# save the plot 
png("output/heatmapMetabolites.png", width = 960, height = 1008)
heatmapPlot
dev.off()

# plotDat_DAPs <- plotDat_order
# save(plotDat_DAPs, file = "plotDat_DAPs.RData")

## select metabolites  --------------------------------
library(openxlsx)
rawDat_iMED_meboAnnot <- read.xlsx('/vol/projects/CIIM/Influenza/iMED/metabolic/raw_data/tables/DATA_CURATED_reformatted.xlsx',
                                   sheet = 'annotation') %>% fill(ionIdx, .direction = "down")

selectedMebos <- rawDat_iMED_meboAnnot %>% filter(Formula %in% selected_vars) # 30 formulas
all_sigMebos <- rawDat_iMED_meboAnnot %>% filter(Formula %in% sigPval_vars) # 180 formulas

write.table(selectedMebos$CompoundID, file = "processedDat/compoundIDs_selectedMebos.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(all_sigMebos$CompoundID, file = "processedDat/compoundIDs_allSigMebos.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)



