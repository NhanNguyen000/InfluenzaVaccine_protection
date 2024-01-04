rm(list = ls())

library(tidyverse)
library(limma)
library(gplots)
library(venn)
library(org.Hs.eg.db)

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}

get.limmaRes <- function(metaDat, inputDat) {
  # Aim: identify DE proteins/metabolites correct with sex, age, and reclassify (the interested vaccine response reclassification groups) 
  #by running the linear model (limma package)
  
  # input: metadata - interested participant (in "probandID" column) with sex, age, and reclassify group, 
  #inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  
  # output: res - limma output
  
  inputDat_temp <- inputDat[, metaDat$probandID]
  
  if (identical(metaDat$probandID, colnames(inputDat_temp)) == TRUE) {
    res <- lmFit(inputDat_temp,
                 design = model.matrix(~ reclassify, metaDat)) %>% 
      eBayes()
  } else res <- "Error: check input"
  
  return(res)
}

get.limmaRes_perStrain <- function(metadat, inputDat, strain_groups) {
  # Aim: run the linear model (using get.limaRes function) with for multiple strain
  
  # input: metadata - interested participant (in "name" column) with sex, age, and reclassify group, 
  # inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  # strain_groups - strains which the test apply to
  
  # output: res - list of outcome from limmar model per strains
  
  resList <- list()
  for (strain_group in strain_groups) {
    metadat_temp <- metadat %>% dplyr::rename("reclassify" := strain_group) %>%
      dplyr::select(probandID, sex, age, reclassify) %>% drop_na()
    
    resList[[strain_group]] <- get.limmaRes(metadat_temp, inputDat)
  }
  return(resList)
}

get.DAs <- function(resLimma) {
  # Aim: extract the DAPs/DAMs with p.value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_DA <- resLimma$p.value %>% as.data.frame %>% 
    dplyr::select(matches("reclassify")) %>% filter(. <0.05)
  
  return(res_DA)
}

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
  # with column name from the given "colName" column, value from given "valName" colum, and rowname from fixed "valName" column
  
  # outcome: plotDat_order - plot data long format with order in valName, so it will show clustering with ggplot heatmap
  
  plotDat_wide <- plotDat %>% dplyr::select(c("valName", colName, valColumn)) %>% 
    pivot_wider(names_from = colName, values_from = valColumn) %>% 
    column_to_rownames("valName")
  
  plotDat_dendrogram <- as.dendrogram(hclust(d = dist(plotDat_wide)))
  plotDat_denOrder <- order.dendrogram(plotDat_dendrogram)
  
  plotDat_order <- plotDat %>%
    mutate(valName = factor(valName, levels = valName[plotDat_denOrder], ordered = TRUE))
  
  return(plotDat_order)
}

# load & prepare data ===========================================
load("cohorts_dat.RData")
load("omicsPred_predTab.RData")

## metadata for corresponding subjects (who have genetic data) -------------------------
metadat <- cohorts$HAI_all %>% filter(probandID %in% rownames(pred_Tab)) %>%
  left_join(cohorts$donorInfo_all %>% dplyr::select(probandID, season, cohort, sex, age)) %>% 
  mutate_at(vars(contains("reclassify")), ~convert_protectees(.x)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "protectee")))

metadat %>% count(H1N1_reclassify)
metadat %>% count(H3N2_reclassify)
metadat %>% count(B_reclassify)
## prepare omicsPred predicted data for lm() model  -------------------------
omicsPred_dat <- pred_Tab %>% rownames_to_column("probandID") %>% 
  arrange(match(probandID, metadat$probandID)) %>%
  column_to_rownames("probandID") %>% t()

# run the limma model  ----------------------------------------------------
strains <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")
#identical(colnames(omicsPred_dat), metadat$probandID) # TRUE, the same sample order

res_omicsPred <- get.limmaRes_perStrain(metadat =  metadat,
                                        inputDat = omicsPred_dat, 
                                        strain_groups = strains)

# significant proteins / metabolites  ----------------------------------------------------
DAs <- res_omicsPred %>% lapply(function(x) get.DAs(x))
tstat_dat <- res_omicsPred %>% lapply(function(x) get.tstat(x))


# check venn diagram 
res.venn <- DAs %>% lapply(function(x) rownames(x))
v.table <- venn(res.venn)

# heatmap  ----------------------------------------------------
## prepared dat --------------------------------
DAs_all <- DAs %>% 
  lapply(function(x) x %>% rownames_to_column("valName")) %>%
  bind_rows(.id = "strain") %>% 
  mutate(strain = gsub("_reclassify", "", strain)) %>%
  dplyr::rename("p.value" = "reclassifyprotectee")

tstat_all <- tstat_dat %>% 
  lapply(function(x) x %>% rownames_to_column("valName")) %>%
  bind_rows(.id = "strain") %>%
  mutate(strain = gsub("_reclassify", "", strain))  %>%
  dplyr::rename("tstat" = "reclassifyprotectee")

tstat_longDat <- tstat_all %>% full_join(DAs_all)

# selected DAs for the heat map
selected_DAs <- unique(as.vector(unlist(res.venn)[which(duplicated(unlist(res.venn)))])) # DAs in at least twice across strains
selected_DAs <- unique((DAs_all %>% filter(p.value < 0.5))$valName)

selected_DAs_info <- selected_models %>% 
  filter(OMICSPRED.ID %in% selected_DAs) %>% # Note: some row have Ensemble ID but no gene name
  left_join(mapIds(org.Hs.eg.db, keys = selected_models$Ensembl.ID, # add missing gene name using org.Hs.eg.db with Ensemble ID
                   column = "SYMBOL", keytype =  "ENSEMBL") %>% 
              unlist() %>% as.data.frame() %>% 
              dplyr::rename("Gene2" = ".") %>% rownames_to_column("Ensembl.ID")) %>%
  mutate(Gene = ifelse(is.na(Gene), Gene2, Gene)) %>% dplyr::select(-Gene2) %>%
  mutate(Gene = ifelse(type != "RNAseq", Gene,
                       ifelse(Ensembl.ID == "ENSG00000237039", "RPS28P4", # manual add missing gene name via internet search
                              ifelse(Ensembl.ID == "ENSG00000227097", "RPS28P7",
                                     ifelse(Ensembl.ID == "ENSG00000226210", "WASH8P",
                                            ifelse(Ensembl.ID == "ENSG00000228409", "CCT6AP1", 
                                                   ifelse(Ensembl.ID == "ENSG00000152117", "SMPD4BP", Gene))))))) %>%
  mutate(Name_all = ifelse(is.na(Name), Gene, # make a column with name of all metabolite, protein, gene
                          ifelse(is.na(CHEMICAL_FORMULA), Name, 
                                 paste0(Name, "(", CHEMICAL_FORMULA, ")"))))


selected_DAs_info %>% count(type, model)  # have checked, no gene/ proteins overlap with DAPs


write.table(unique(selected_DAs_info$Gene), file = "processedDat/omicsPred/sig_1393genes.txt", 
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

write.table(unique(selected_DAs_info$Gene), file = "processedDat/omicsPred/sig_98genes.txt", 
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

## selectDA elements with consistent trend across 3 strain --------------------
selected_DAs_v2 <- tstat_longDat %>% filter(valName %in% selected_DAs) %>%
  mutate(direction = ifelse(tstat >0, 1, -1)) %>%
  group_by(valName) %>% summarise(regula_trend = sum(direction)) %>% # the regulation direction across 3 strains
  filter(abs(regula_trend) == 3) %>% # select the consistent regulation direction across 3 strains (= +/-3)
  dplyr::select(valName) %>% unlist()

selected_DAs_info %>% filter(OMICSPRED.ID %in% selected_DAs_v2) %>% count(type, model)
selected_DAs_info_v2 <- selected_DAs_info %>% filter(OMICSPRED.ID %in% selected_DAs_v2)

write.table(unique(selected_DAs_info_v2$Gene), file = "processedDat/omicsPred/sig_59genes.txt", 
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)
## selectDA elements (protein, RNAseq) show positive correlation with measured values --------------------
load("pred_with_proRna_measures.RData")
predMeasure_positivCor <- pred_inf %>% filter(cor > 0)

selected_DAs_info_v2.1 <- selected_DAs_info_v2 %>% #as.data.frame() %>% 
  filter(Gene %in% predMeasure_positivCor$Gene) %>% filter(model %in% c("RNAseq", "Olink"))

selected_DAs_info_v2.1 %>% count(type, model)

## plot heatmap --------------------------------
plotDat <- tstat_longDat %>% 
  filter(valName %in% selected_DAs_v2) %>%
 # filter(valName %in% selected_DAs_info_v2.1$OMICSPRED.ID) %>%
  left_join(selected_DAs_info %>% dplyr::select(OMICSPRED.ID, Name_all, model),
            by = c("valName" = "OMICSPRED.ID")) %>%
  mutate(Name_all = ifelse(Name_all != "ICAM5", Name_all, # The ICAM5 proteins in Somalogic has 2 prediction models
                           paste(Name_all, "_", valName))) %>% 
  dplyr::select(-valName) %>% rename("Name_all" = "valName") %>% as.data.frame()

plotDat_order <- get.plotDat_clusterRow(plotDat, 
                                        colName = "strain", 
                                        valColumn = "tstat")

plotDat_order %>%
  ggplot(aes(x = strain, y = valName, fill = tstat)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme( #axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text=element_text(size=8))

# check pathway ---------------------
test <- tstat_longDat %>% filter(valName %in% unique(unlist(res.venn))) %>%
  left_join(selected_models, by = c("valName" = "OMICSPRED.ID")) %>% 
  filter(type %in% c("RNAseq", "protein")) %>%
 # dplyr::select(Gene, type, p.value, strain) %>% drop_na() %>%
  dplyr::select(Gene) %>% drop_na() %>%
  filter(Gene != "") %>% distinct()

write.table(test$Gene, 
            file = "omicsPred_sigGenes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# Adapt code from Xun's pathway analysis ---------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)

#genenames <- unique(selected_DAs_info$Gene) # sig. gene or protein in at least 2 strains --> no pathway result
genenames <- test$Gene # sig. gene or protein in at least 1 strain --> pathway result

pval <- 0.05
qval <- 0.05

er1 <- data.frame(bitr(geneID = genenames, 
                       fromType = 'SYMBOL', toType = c('ENTREZID', 'ENSEMBL'),
                       OrgDb = org.Hs.eg.db, drop = T))

enr.GO.er1 <- enrichGO(gene = er1$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "fdr", #none or fdr
                       pvalueCutoff = pval,
                       qvalueCutoff = qval,
                       readable = TRUE) 

enr.KEGG.er1 <- enrichKEGG(gene = er1$ENTREZID,
                           organism = "hsa", 
                           keyType = "kegg", 
                           pAdjustMethod = "fdr",
                           pvalueCutoff = pval,
                           qvalueCutoff = qval, 
                           use_internal_data = FALSE)

enr.MF.er1 <- enrichGO(gene = er1$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "MF",
                       pAdjustMethod = "fdr", #none or fdr
                       pvalueCutoff = pval,
                       qvalueCutoff = qval,
                       readable = TRUE)

res = list(
  'GO' = enr.GO.er1,
  'KEGG' = enr.KEGG.er1,
  'MF' = enr.MF.er1
)

dotplot(merge_result(res), showCategory = 10) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 1))

# check the gene/protein have mulitple predictiom model
test2 <- tstat_longDat %>% filter(valName %in% unique(unlist(res.venn))) %>%
  left_join(selected_models, by = c("valName" = "OMICSPRED.ID")) %>% 
  filter(type %in% c("RNAseq", "protein")) %>%
  dplyr::select(Gene, type, p.value, strain, valName) %>% drop_na() %>%
  filter(Gene != "") %>% distinct()

a <- test2 %>% add_count(Gene) %>% select(Gene, n) %>% distinct()
a <- test2 %>% add_count(Gene) 

res$GO@result[which(res$GO@result$Description == "leukocyte mediated cytotoxicity"), ]
res$GO@result[which(res$GO@result$Description == "natural killer cell mediated immunity"), ]

# check the protein
load("cohorts_dat.RData")

View(protein_Dat$iMED_2015)
intersect(names(protein_Dat$iMED_2015), test$Gene) # check TNFSF12, FCRL3
intersect(test$Gene, unlist(res.venn$iMED_2015) %>% as.vector()) # got some gene significants
