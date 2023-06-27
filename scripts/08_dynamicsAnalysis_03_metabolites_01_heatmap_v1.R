rm(list = ls())

library(tidyverse)

get.plotDat_clusterRow <- function(plotDat, colName, valColumn) {
  # Aim: make the heatmap data with dendrogram clustering order, so can have clustered heatmap with ggplot
  # this function comvert the plot data (long format) to wide format 
  # with column name from the given "colName" column, value from given "valName" colum, and rowname from fixed "valName" column
  
  # outcome: plotDat_order - plot data long format with order in valName, so it will show clustering with ggplot heatmap
  
  plotDat_wide <- plotDat %>% select(c("valName", colName, valColumn)) %>% 
    pivot_wider(names_from = colName, values_from = valColumn) %>% 
    column_to_rownames("valName")
  
  plotDat_dendrogram <- as.dendrogram(hclust(d = dist(plotDat_wide)))
  plotDat_denOrder <- order.dendrogram(plotDat_dendrogram)
  
  plotDat_order <- plotDat %>%
    mutate(valName = factor(valName, levels = valName[plotDat_denOrder], ordered = TRUE))
  
  return(plotDat_order)
}
# load data =======================================================================
load("/vol/projects/CIIM/Influenza/iMED/metabolic/metabolite_reclassify_timeDynamics.RData")
load("cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- mebo_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

# sig. metabolites  ================================
allMeboChangingSigni <- list(
  "LL_T2vsT1" = rownames(LL_T2vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "LL_T3vsT1" = rownames(LL_T3vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "LH_T2vsT1" = rownames(LH_T2vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "LH_T3vsT1" = rownames(LH_T3vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "HL_T2vsT1" = rownames(HL_T2vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "HL_T3vsT1" = rownames(HL_T3vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "HH_T2vsT1" = rownames(HH_T2vsT1_metabol %>% filter(adj.P.Val < 0.05)),
  "HH_T3vsT1" = rownames(HH_T3vsT1_metabol %>% filter(adj.P.Val < 0.05)))

mebos <- unique(unlist(allMeboChangingSigni))

plotDat <- inputDat %>% 
  select(probandID, season, responder, time, c(mebos),
         matches("_abFC|_T1|_T4|_reclassify"))

# heatmap based on expression -------------------------------------------------
## all season -------------------------
plotDat_H1N1 <- plotDat %>% 
  pivot_longer(cols = mebos, names_to = "valName", values_to =  "intensity") %>%
  mutate(group = paste0(season, "_", H1N1_reclassify, "_", time)) %>%
  mutate(group = factor(group, 
                        levels = c("2014_LL_T1", "2014_LL_T3", "2014_LL_T4",
                                   "2014_LH_T1", "2014_LH_T3", "2014_LH_T4",
                                   "2014_HL_T1", "2014_HL_T3", "2014_HL_T4",
                                   "2014_HH_T1", "2014_HH_T3", "2014_HH_T4",
                                   "2015_LL_T1", "2015_LL_T3", "2015_LL_T4",
                                   "2015_LH_T1", "2015_LH_T3", "2015_LH_T4",
                                   "2015_HL_T1", "2015_HL_T3", "2015_HL_T4",
                                   "2015_HH_T1", "2015_HH_T3", "2015_HH_T4",
                                   "2019_LL_T1", "2019_LL_T4", "2019_LL_T5",
                                   "2019_LH_T1", "2019_LH_T4", "2019_LH_T5",
                                   "2019_HL_T1", "2019_HL_T4", "2019_HL_T5",
                                   "2019_HH_T1", "2019_HH_T4", "2019_HH_T5",
                                   "2020_LL_T1", "2020_LL_T4", "2020_LL_T5",
                                   "2020_LH_T1", "2020_LH_T4", "2020_LH_T5",
                                   "2020_HL_T1", "2020_HL_T4", "2020_HL_T5")))

plotDat_H1N1 %>%
  ggplot(aes(x = group, y = valName, fill = intensity)) + 
  geom_tile() +
  scale_fill_gradientn(limits = c(10, 23), colors = c("blue", "white", "red"),  na.value = "grey") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  ggtitle("Sig. metabolite in dynamic analysis, show in H1N1 strain")
# heatmap based on t-stat -------------------------------------------------
## seasons 2015 ---------------------------------
allMeboChangingTstat <- list(
  "LL_T2vsT1" = LL_T2vsT1_metabol %>% select(t),
  "LL_T3vsT1" = LL_T3vsT1_metabol %>% select(t),
  "LH_T2vsT1" = LH_T2vsT1_metabol %>% select(t),
  "LH_T3vsT1" = LH_T3vsT1_metabol %>% select(t),
  "HL_T2vsT1" = HL_T2vsT1_metabol %>% select(t),
  "HL_T3vsT1" = HL_T3vsT1_metabol %>% select(t),
  "HH_T2vsT1" = HH_T2vsT1_metabol %>% select(t),
  "HH_T3vsT1" = HH_T3vsT1_metabol %>% select(t)) %>% 
  lapply(function(x) x %>% rownames_to_column("valName")) %>%
  bind_rows(.id = "group")

plotDat <- allMeboChangingTstat %>% filter(valName %in% mebos) %>%
  mutate(group = factor(group, 
                        levels = c("LL_T2vsT1", "LL_T3vsT1", 
                                   "LH_T2vsT1", "LH_T3vsT1", 
                                   "HL_T2vsT1", "HL_T3vsT1", 
                                   "HH_T2vsT1", "HH_T3vsT1")))

plotDat %>%
  ggplot(aes(x = group, y = valName, fill = t)) + 
  geom_tile() +
  scale_fill_gradientn(limits = c(-8, 8), colors = c("blue", "white", "red"),  na.value = "grey") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

plotDat_order <- get.plotDat_clusterRow(plotDat, 
                                        colName = "group", 
                                        valColumn = "t")
plotDat_order %>%
  ggplot(aes(x = group, y = valName, fill = t)) + 
  geom_tile() +
  scale_fill_gradientn(limits = c(-8, 8), colors = c("blue", "white", "red"),  na.value = "grey") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
