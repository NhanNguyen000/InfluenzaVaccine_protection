
# load data from iMED transcriptome ---------------
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
View(transcriptomeList[[1]])
View(transcriptomeList[[2]])
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1")
iMED_transcripDat <- transcriptomeList[[1]] %>% as.data.frame %>% 
  select(iMED_trancrip_T1$SampleName)


# metadata from reclassification --------------------
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("young", "old","cirrhosis")))

iMED_transcrip_metadat <- iMED_transcrip_T1 %>% left_join(metadat)
iMED_transcrip_metadat %>% count(H1N1_reclassify)
iMED_transcrip_metadat %>% count(H3N2_reclassify)
iMED_transcrip_metadat %>% count(B_reclassify)

# heatmap of DEPs ------------------------
# LL vs. Other
strain_group <- "H1N1_reclassify"
strain_group <- "H3N2_reclassify" # have no FABP9 and FGF19, tried different names using Genecards
strain_group <- "B_reclassify"

metadat_heatmap <- iMED_transcrip_metadat %>% 
  rename("reclassify" := strain_group) %>% 
  # mutate(class = ifelse(reclassify == "LL", "LL", "Others")) %>% 
  mutate(class = reclassify) %>% filter(class %in% c("LL", "LH")) %>%
  select(class, SampleName)

# proteins <- get.proteinName(OlinkIDs = rownames(resSig_LLvsProtector[[strain_group]]), 
#                             OlinkDat = cohorts_dat$proteinAnnot_all) 
proteins <- get.proteinName(OlinkIDs = rownames(resSig_LLvsLH[[strain_group]]), 
                            OlinkDat = cohorts_dat$proteinAnnot_all) 
proteins <- proteins[-which(proteins %in% c("FABP9", "FGF19"))]

transcripDat_temp <- iMED_transcripDat %>% t() %>% as.data.frame() %>%
  select(proteins) %>% rownames_to_column("SampleName")

heatmapDat <- metadat_heatmap %>% 
  left_join(transcripDat_temp) %>%
  group_by(class) %>%
  summarise_at(vars(-SampleName), funs(mean(., na.rm = TRUE))) %>%
  column_to_rownames("class") %>% t()

heatmapDat %>% Heatmap()
heatmapDat %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))
heatmapDat %>% t() %>% scale() %>% t() %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))
