rm(list = ls())

library(tidyverse)
library(patchwork)
library(Seurat)

get.DAs <- function(resLimma) {
  # Aim: extract the DAPs/DAMs with p.value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_DA <- resLimma$p.value %>% as.data.frame %>% 
    select(matches("reclassify")) %>% filter(. <0.05)
  
  return(res_DA)
}

get.tstat <- function(resLimma) {
  # Aim: extract the t-statistic value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_tstat <- resLimma$t %>% as.data.frame %>% 
    select(matches("reclassify"))
  
  return(res_tstat)
}
# DAPs at protein level, heatmap ------------------------------------------------
load("plotDat_DAPs.RData")

plot_DAPs <- plotDat_DAPs %>%
  ggplot(aes(x = group, y = valName, fill = tstat)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(p.value), NA, "*"))) +
  scale_fill_gradientn(limits = c(-5, 5), colors = c("blue", "white", "red")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  ggtitle("Protein")


# DAPs at bulk RNAseq level, heatmap ------------------------------------------------
#load("selected_DAPs.RData")
load("resPro_2015_bulkRNAseq.RData")

DAs_all_bulkRNAseq <- resPro_2015_bulkRNAseq %>% 
  lapply(function(x) get.DAs(x) %>% rownames_to_column("valName"))  %>%
  bind_rows(.id = "strain") %>% mutate(strain = gsub("_reclassify", "", strain)) %>% 
  rename("p.value" = "reclassifyprotectee") %>%
  mutate(season = "2015")  %>%
  mutate(group = paste0(season, "_", strain))

tstat_all_bulkRNAseq <- resPro_2015_bulkRNAseq %>% 
  lapply(function(x) get.tstat(x) %>% rownames_to_column("valName")) %>%
  bind_rows(.id = "strain") %>% mutate(strain = gsub("_reclassify", "", strain)) %>% 
  rename("tstat" = "reclassifyprotectee") %>%
  mutate(season = "2015")  %>%
  mutate(group = paste0(season, "_", strain))

plotDat_DAPs_bulkRNAseq <- tstat_all_bulkRNAseq %>% 
  left_join(DAs_all_bulkRNAseq) %>%
  filter(valName %in% levels(plotDat_DAPs$valName)) %>% 
  mutate(valName = factor(valName, levels = levels(plotDat_DAPs$valName), ordered = TRUE))

plot_DAPs_bulkRNAseq <- plotDat_DAPs_bulkRNAseq %>%
  ggplot(aes(x = group, y = valName, fill = tstat)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(p.value), NA, "*"))) +
  scale_fill_gradientn(limits = c(-5, 5), colors = c("blue", "white", "red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  ggtitle("Bulk RNAseq")

# DAPs at scRNAseq level (scRNAseq from Duesseldorf), heatmap ------------------------------------------------
load("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/analysis/output/data.RData") # Note: marker.list is based on accumulated previous information (no specific papers)
metadata <- read.csv2("/vol/projects/CIIM/Influenza/Duesseldorf/Influenza_Duesseldorf/scRNAseq/data/metadata.csv")

unique(data$vaccination)
unique(data$celltype)

df_subset <- FetchData(object = data, 
                vars = c(levels(plotDat_DAPs$valName), c("sample", "vaccination", 'celltype')), 
                slot = 'data') %>%
     filter(vaccination ==  "before") # subset only scRNAseq before vaccination
names(df_subset)

DAPs_scRNAseq <- df_subset %>%
  select(-sample, -vaccination) %>%
  group_by(celltype) %>% summarise_all("mean") %>% # calculate the average gene expression per cell type
  pivot_longer(!celltype, names_to = "valName", values_to = "geneExp") %>%
  mutate(valName = factor(valName, levels = levels(plotDat_DAPs$valName), ordered = TRUE))


plot_DAPs_scRNAseq <- DAPs_scRNAseq %>%
  ggplot(aes(x = celltype, y = valName, fill = geneExp)) + 
  geom_tile() +
 # geom_text(aes(label = ifelse(is.na(p.value), NA, "*"))) +
  scale_fill_gradientn(limits = c(0, 2), colors = c("white", "red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  ggtitle("ScRNAseq")

# merge heatmaps together -------------------------------------
plot_DAPs + plot_DAPs_bulkRNAseq + plot_DAPs_scRNAseq +
  plot_layout(ncol = 3, widths = c(4, 1.5, 4), guides = "collect") & 
  theme(legend.position = "right")

