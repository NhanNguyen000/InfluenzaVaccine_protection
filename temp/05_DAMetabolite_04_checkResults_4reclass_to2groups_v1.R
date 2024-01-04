rm(list = ls())

library(tidyverse)
library(venn)

get.DAs <- function(resLimma) {
  # Aim: extract the DAPs/DAMs with p.value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_DA <- resLimma$p.value %>% as.data.frame %>% 
    select(matches("reclassify")) %>% filter(. <0.05)
  
  return(res_DA)
}

get.p.value <- function(resLimma) {
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
load("resMebo_4reclass_2group.RData")

# significant proteins / metabolites ---------------------------------------------------------
pvals <- resMebo_4reclass_2group %>% 
  lapply(function(x) x %>% lapply(function(y) get.p.value(y))) 

DAs <- resMebo_4reclass_2group %>% 
  lapply(function(x) x %>% lapply(function(y) get.DAs(y)))

tstat_dat <- resMebo_4reclass_2group %>%
  lapply(function(x) x %>% lapply(function(y) get.tstat(y)))

# Check venn diagram ---------------------------------------------
res.venn <- DAs %>% lapply(function(x) x %>% lapply(function(y) rownames(y)))
v.table <- venn(res.venn$iMED_2014)
v.table <- venn(res.venn$iMED_2015)
v.table <- venn(res.venn$ZirFlu_2019)
v.table <- venn(res.venn$ZirFlu_2020)

a <- res.venn$iMED_2015
names(a) <- c("H1N1_2015", "H3N2_2015", "B_2015")
a$H1N1_2014 <- res.venn$iMED_2014$H1N1_reclassify
a$B_2014 <- res.venn$iMED_2014$B_reclassify
a$H1N1_2019 <- res.venn$ZirFlu_2019$H1N1_reclassify
a$H1N1_2020 <- res.venn$ZirFlu_2020$H1N1_reclassify
a$Byamagata_2020 <- res.venn$ZirFlu_2020$Byamagata_reclassify
v.table <- venn(a)

# check Padj significant ----------------------------------------
padj_2015 <- resMebo_4reclass_2group$iMED_2015 %>% 
  lapply(function(x) x %>% 
           eBayes %>% topTable(n = Inf, sort.by = "n") %>% 
           filter(adj.P.Val < 0.05)) 

res.venn_padj_2015 <- padj_2015 %>% 
  lapply(function(x) x %>% rownames(x))

v.table <- venn(res.venn_padj_2015)

# # the strain in seasons have p-adj significant (<0.05)
resMebo_padj <- list()
resMebo_padj[["2014_H1N1"]] <- resMebo_4reclass_2group$iMED_2014$H1N1_reclassify
resMebo_padj[["2014_B"]] <- resMebo_4reclass_2group$iMED_2014$B_reclassify
resMebo_padj[["2019_H1N1"]] <- resMebo_4reclass_2group$ZirFlu_2019$H1N1_reclassify
resMebo_padj[["2020_Byamagata"]] <- resMebo_4reclass_2group$ZirFlu_2020$Byamagata_reclassify

padj_otherSeasons <- resMebo_padj %>% 
  lapply(function(x) x %>% 
           eBayes %>% topTable(n = Inf, sort.by = "n") %>% 
           filter(adj.P.Val < 0.05))

res.venn_padj_otherSeason <- padj_otherSeasons %>% 
  lapply(function(x) x %>% rownames(x))

v.table <- venn(res.venn_padj_otherSeason)
res.venn_padj_all <- intersect(res.venn_padj_2015 %>% unlist %>% as.vector, res.venn_padj_otherSeason %>% unlist %>% as.vector)

res.venn_adj <- c(res.venn_padj_otherSeason)
# heatmap  ----------------------------------------------------
DAs_all <- DAs %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% rownames_to_column("valName")) %>%
           bind_rows(.id = "strain")) %>%
  bind_rows(.id = "season") %>%
  mutate(strain = gsub("_reclassify", "", strain)) %>%
  rename("p.value" = "reclassifyprotectee")

tstat_all <- tstat_dat %>% 
  lapply(function(x) x %>% 
           lapply(function(y) y %>% rownames_to_column("valName")) %>%
           bind_rows(.id = "strain")) %>%
  bind_rows(.id = "season") %>%
  mutate(strain = gsub("_reclassify", "", strain))  %>%
  rename("tstat" = "reclassifyprotectee")

tstat_longDat <- tstat_all %>%
  full_join(DAs_all) %>%
  separate(season, sep = "_", into = c("cohort", "season")) %>%
  mutate(group = paste0(season, "_", strain))

# selected DAs for the heat map
# selected_DAs <- unique(c(intersect(res.venn$iMED_2015$H1N1_reclassify, res.venn$iMED_2015$H3N2_reclassify),
#                          intersect(res.venn$iMED_2015$B_reclassify, res.venn$iMED_2015$H3N2_reclassify),
#                          intersect(res.venn$iMED_2015$H1N1_reclassify, res.venn$iMED_2015$B_reclassify))) # DAs in at least 2 strains for iMED_2015
# selected_DAs <- unique(unlist(res.venn$iMED_2015))
selected_DAs <- unique(as.vector(unlist(res.venn)[which(duplicated(unlist(res.venn)))])) # DAs in at least twice acorss strain and season
#save(selected_DAs, file = "selected_DAMs.RData")

selected_DAs <- unique(unlist(res.venn_padj_2015))
selected_DAs <- res.venn_padj_all
# heatmap
plotDat <- tstat_longDat %>% filter(valName %in% selected_DAs) %>%
  slice(-which(group %in% c("2014_H3N2", 
                            "2019_Bvictoria", "2019_Byamagata", "2019_H3N2", 
                            "2020_Bvictoria", "2020_H3N2"))) # remove some group if needed

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
