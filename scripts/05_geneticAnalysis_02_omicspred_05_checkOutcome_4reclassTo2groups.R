rm(list = ls())

library(tidyverse)
library(gplots)
library(venn)
library(org.Hs.eg.db)

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}

get.plotDat_clusterRow <- function(plotDat, colName, valColumn) {
  # Aim: make the heatmap data with dendrogram clustering order, so can have clustered heatmap with ggplot
  # this function comvert the plot data (long format) to wide format 
  # with column name from the given "colName" column, value from given "valName" colum, and rowname from fixed "valName" column
  
  # outcome: plotDat_order - plot data long format with order in valName, so it will show clustering with ggplot heatmap
  
  plotDat_wide <- plotDat %>% dplyr::select(c("varName", colName, valColumn)) %>% 
    pivot_wider(names_from = colName, values_from = valColumn) %>% 
    column_to_rownames("varName")
  
  plotDat_dendrogram <- as.dendrogram(hclust(d = dist(plotDat_wide)))
  plotDat_denOrder <- order.dendrogram(plotDat_dendrogram)
  
  plotDat_order <- plotDat %>%
    mutate(varName = factor(varName, levels = varName[plotDat_denOrder], ordered = TRUE))
  
  return(plotDat_order)
}

# load & prepare data ===========================================
load("processedDat/cohorts_dat.RData")
load("processedDat/omicsPred_output.RData")

## metadata for corresponding subjects (who have genetic data) -------------------------
metadat <- cohorts$HAI_all %>% filter(probandID %in% rownames(pred_Tab)) %>%
  left_join(cohorts$donorInfo_all %>% dplyr::select(probandID, season, cohort, sex, age)) %>% 
  mutate_at(vars(contains("reclassify")), ~convert_protectees(.x)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "protectee")))

metadat %>% count(H1N1_reclassify)
metadat %>% count(H3N2_reclassify)
metadat %>% count(B_reclassify)

rm(pred_Tab)

# omicsPred outcome  ------------------------------------------------------------
load("processedDat/res_omicsPred_4reclass_2group.RData")

## significant predicted models  ------------------------------------------------
DEs <- res_omicsPred %>% 
  lapply(function(x) x %>% filter(p_value < 0.05) %>%
           dplyr::select(model_id, p_value, padj))

# venn diagram plot
res.venn <- DEs %>% lapply(function(x) x$model_id)

png("output/omicPred_vennPlot_compare4reclassTo2groups.png", width = 512)
v.table <- venn(res.venn, plotsize = 10, ilcs = 1.5, sncs = 1.5, box = FALSE)
dev.off()


## check for padj
DEs_padj <- res_omicsPred %>% 
  lapply(function(x) x %>% filter(padj < 0.05) %>%
           dplyr::select(model_id, p_value, padj))

res.venn_padj <- DEs_padj %>% lapply(function(x) x$model_id)
v.table <- venn(res.venn_padj)

# heatmap  -------------------------------------------------------------------
tstat_longDat <- res_omicsPred %>% 
  bind_rows(.id = "strain") %>% 
  mutate(strain = gsub("_reclassify", "", strain)) 

# selected DEs for the heat map
selected_DEs_temp <- unique(as.vector(unlist(res.venn)[which(duplicated(unlist(res.venn)))])) # DAs in at least twice across strains

selected_DEs <- tstat_longDat %>% # selectDA elements with consistent trend across 3 strain
  filter(model_id %in% selected_DEs_temp) %>% 
  mutate(direction = ifelse(t_statistic >0, 1, -1)) %>%
  group_by(model_id) %>% summarise(regula_trend = sum(direction)) %>% # the regulation direction across 3 strains
  filter(abs(regula_trend) == 3) %>% # select the consistent regulation direction across 3 strains (= +/-3)
  dplyr::select(model_id) %>% unlist()

selected_DEs_info <- selected_models %>% 
  filter(OMICSPRED.ID %in% selected_DEs) %>% # Note: some row have Ensemble ID but no gene name
  left_join(mapIds(org.Hs.eg.db, keys = selected_models$Ensembl.ID, # add missing gene name using org.Hs.eg.db with Ensemble ID
                   column = "SYMBOL", keytype =  "ENSEMBL") %>% 
              unlist() %>% as.data.frame() %>% 
              dplyr::rename("Gene2" = ".") %>% rownames_to_column("Ensembl.ID")) %>%
  mutate(Gene = ifelse(is.na(Gene), Gene2, Gene)) %>% dplyr::select(-Gene2) %>%
  mutate(Name_all = ifelse(is.na(Gene),  # make a column with name of all metabolite, protein, gene
                           ifelse(is.na(CHEMICAL_FORMULA), Name, CHEMICAL_FORMULA), Gene)) %>%
  mutate(Name_all = ifelse(nchar(Name_all) > 0, Name_all, 
                           ifelse(Ensembl.ID == "ENSG00000215146", "RP11-313J2.1", # manual add missing gene name via internet search
                                  ifelse(Ensembl.ID == "ENSG00000226210", "WASH8P",
                                         ifelse(Ensembl.ID == "ENSG00000260404", "ENSG00000260404", # no gene name from Ensembls
                                                ifelse(Ensembl.ID == "ENSG00000242299", "RPS18P5", NA))))))

selected_DEs_info %>% count(type, model) 

## plot heatmap ----------------------------------------------------------------
plotDat <- tstat_longDat %>% 
  filter(model_id %in% selected_DEs) %>%
  left_join(selected_DEs_info %>% dplyr::select(OMICSPRED.ID, Name_all, model),
            by = c("model_id" = "OMICSPRED.ID")) %>%
  dplyr::select(-model_id) %>% rename("Name_all" = "varName") %>% as.data.frame() %>% 
  mutate(t_statistic = as.numeric(t_statistic), 
         p_value = as.numeric(p_value),
         padj = as.numeric(padj))

plotDat_order <- get.plotDat_clusterRow(plotDat, 
                                        colName = "strain", 
                                        valColumn = "t_statistic")

heatmapPlot <- plotDat_order %>%
  ggplot(aes(x = strain, y = varName, fill = t_statistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(padj < 0.05, "**", NA)), size = 10) +
  geom_text(aes(label = ifelse(p_value < 0.05, "*", NA)), size = 10) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(text = element_text(size = 24))

# save the plot 
png("output/omicsPred_heatmapPredictedModel.png", width = 720, height = 1632)
heatmapPlot
dev.off()

# potential genes/proteins for discussions --------------------------------
# Gene (RNAseq models) and proteins (Somalogic model) that are ppadj < 0.05 in H3N2s strains
TMEM51
IFT52
CAMK1 # Somalogic
ARSA
CCL7 # Somalogic
TMEM204
ACVR2B
RNMT

# gene and proteins that are overlapped across 3 strainss
dat <- selected_DEs_info %>% 
  filter(OMICSPRED.ID %in% intersect(intersect(DEs$H1N1_reclassify$model_id, DEs$H3N2_reclassify$model_id), 
               DEs$B_reclassify$model_id))

IFT52 # appear

# DE proteins predicted by olink model
CDH5
IL18BP

