rm(list = ls())

library(tidyverse)
library(venn)

get.DAs <- function(resLimma) {
  # Aim: extract the DAPs/DAMs with p.value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_DA_temp <- resLimma$p.value %>% as.data.frame %>% 
    select(matches("reclassify")) %>% mutate_all(~ifelse(. >= 0.05, NA, .)) %>%
    slice(-which(rowSums(is.na(.)) == 3))
  
  res_DA <- list()
  res_DA$LLvsLH <- res_DA_temp %>% select(matches("LH")) %>% drop_na()
  res_DA$LLvsHL <- res_DA_temp %>% select(matches("HL")) %>% drop_na()
  res_DA$LLvsHH <- res_DA_temp %>% select(matches("HH")) %>% drop_na()
  
  return(res_DA)
}

get.tstat <- function(resLimma) {
  # Aim: extract the t-statistic value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_tstat <- resLimma$t %>% as.data.frame %>% 
    select(matches("reclassify"))
  
  return(res_tstat)
}

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

# load data --------------------------------------------------------
load("resPro_4reclass.RData")

# significant proteins / metabolites ---------------------------------------------------------
DAs <- resPro_4reclass %>% 
  lapply(function(x) x %>% lapply(function(y) get.DAs(y)))

tstat_dat <- resPro_4reclass %>%
  lapply(function(x) x %>% lapply(function(y) get.tstat(y)))

# Check venn diagram ---------------------------------------------
res.venn <- DAs %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% lapply(function(z) z %>% rownames(z))))

v.table <- venn(res.venn$iMED_2014$H1N1_reclassify)
v.table <- venn(res.venn$iMED_2015$H3N2_reclassify)
v.table <- venn(res.venn$ZirFlu_2019$Bvictoria_reclassify)
v.table <- venn(res.venn$ZirFlu_2020$Byamagata_reclassify)

# heatmap  ----------------------------------------------------
DAs_all <- DAs %>% 
  lapply(function(x) 
    x %>% lapply(function(y) 
      y %>% lapply(function(z) 
        z %>% rownames_to_column("valName")) %>% purrr::reduce(full_join)) %>%
      bind_rows(.id = "strain")) %>%
  bind_rows(.id = "season") %>%
  mutate(strain = gsub("_reclassify", "", strain)) %>%
  gather(key = "comparison", value = "p.value", matches("reclassify")) %>% drop_na() %>%
  mutate(comparison = gsub("reclassify", "_LLvs", comparison))

tstat_all <- tstat_dat %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% rownames_to_column("valName")) %>%
           bind_rows(.id = "strain")) %>%
  bind_rows(.id = "season") %>%
  mutate(strain = gsub("_reclassify", "", strain))  %>%
  gather(key = "comparison", value = "tstat", matches("reclassify")) %>%
  mutate(comparison = gsub("reclassify", "_LLvs", comparison))

tstat_longDat <- tstat_all %>%
  full_join(DAs_all) %>%
  mutate(group = paste0(season, "_", strain, comparison)) %>%
  separate(season, sep = "_", into = c("cohort", "season"))

# selected DAs for the heat map
selected_DAs <- unique(c(intersect(res.venn$iMED_2015$H1N1_reclassify, res.venn$iMED_2015$H3N2_reclassify),
                         intersect(res.venn$iMED_2015$B_reclassify, res.venn$iMED_2015$H3N2_reclassify),
                         intersect(res.venn$iMED_2015$H1N1_reclassify, res.venn$iMED_2015$B_reclassify))) # DAs in at least 2 strains

selected_DAs <- unique(unlist(res.venn$iMED_2015))
selected_DAs <- unique(unlist(res.venn_padj))
#save(selected_DAs, file = "selected_DAPs_padj2015.RData")
# heatmap
plotDat <- tstat_longDat %>% 
  filter(valName %in% selected_DAs) %>%
  filter(cohort == "iMED") # pick cohort / comparison needed

plotDat_order <- get.plotDat_clusterRow(plotDat, 
                                        colName = "group", 
                                        valColumn = "tstat")

plotDat_order %>%
  ggplot(aes(x = group, y = valName, fill = tstat)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

## heatmap - only season iMED 2015 -------------------------
plotDat <- tstat_longDat %>% 
  filter(valName %in% selected_DAs) %>%
  filter(cohort == "iMED", season == "2015") %>%
  mutate(group = gsub("iMED_2015_", "", group), comparison = gsub("_", "", comparison)) %>%
  mutate(group = factor(group, levels = c("H1N1_LLvsLH", "H1N1_LLvsHL", "H1N1_LLvsHH",
                                          "H3N2_LLvsLH", "H3N2_LLvsHL", "H3N2_LLvsHH",
                                          "B_LLvsLH", "B_LLvsHL", "B_LLvsHH"))) %>%
  rename("relative_diff" = "tstat")

plotDat_order <- get.plotDat_clusterRow(plotDat, 
                                        colName = "group", 
                                        valColumn = "relative_diff")

plotDat_order %>%
  ggplot(aes(x = group, y = valName, fill = relative_diff)) + 
  geom_tile() +
  # geom_text(aes(label = ifelse(is.na(p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

plotDat_order %>% 
  ggplot(aes(x = group, y = valName, fill = relative_diff)) + 
  geom_tile() +
  # geom_text(aes(label = ifelse(is.na(p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw()

# check padj if possible -----------------------
get.DAs_v2 <- function(resLimma) {
  # Aim: extract the DAPs/DAMs with p.value from linnear comparision comparison
  # following the get.limmaRes result
  res_DA <- resLimma$p.value %>% as.data.frame %>% 
    select(matches("reclassify"))
  
  return(res_DA)
}

tstat_longDat_v2 <- tstat_longDat %>% 
 # filter(season == "2015") %>%
  mutate(direction = ifelse(tstat < 0, "down", ifelse(tstat > 0, "up", tstat))) %>%
  group_by(valName) %>% add_count(direction) %>%
  filter(n == 9)
a <- unique(tstat_longDat_v2$valName)


DAs_v3 <- DAs %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% 
                    rownames_to_column("valName") %>% 
                    filter(valName %in% a)))

a <- DAs_v3 %>% lapply(function(x) x %>%
                         lapply(function(y) y %>%
                                  mutate(padj = p.adjust(reclassifyprotectee, method = "fdr"))))
b <- a %>% lapply(function(x) x %>%
                    lapply(function(y) y %>% filter(padj < 0.05))) # CD83 sig. padj at H1N1 and B strain 2015

