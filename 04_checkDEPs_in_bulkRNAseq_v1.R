library(ggpubr)

# load data from iMED transcriptome ---------------
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
# View(transcriptomeList[[1]])
# View(transcriptomeList[[2]])
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1")
iMED_transcripDat <- transcriptomeList[[1]] %>% as.data.frame %>% 
  select(iMED_transcrip_T1$SampleName)


# metadata from reclassification --------------------
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("young", "old","cirrhosis")))

# boxplot of DEPs ------------------------
# have no FABP9 and FGF19, tried different names using Genecards

plotDat <- iMED_transcrip_T1 %>% left_join(metadat) %>%
  left_join(iMED_transcripDat %>% t() %>% as.data.frame %>% 
              rownames_to_column("SampleName"))

plotDat %>% count(H1N1_reclassify)
plotDat %>% count(H3N2_reclassify)
plotDat %>% count(B_reclassify)

# compare re-classfiy groups
my_comparisons_v2.1 <- list( c("LL", "LH"), c("LL", "HL"))
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))

strain <- "H1N1_reclassify"
strain <- "H3N2_reclassify"
strain <- "B_reclassify"
protein <- "CD83"
protein <- "SERPINB8"
# "LSP1", "TNFRSF11A", "KRT19", "CSF3", "CXCL17", "SPINK4", "NBN"
protein <- "VEGFA"
protein <- "TNF"
protein <- "BTN3A2"
protein <- "COLEC12"
protein <- "CLEC7A"
protein <- "IL6"
protein <- "REG4"
protein <- "PLAUR"
protein <- "CD160"
protein <- "PGF"
protein <- "GZMA" #(?)
protein <- "NCR1"
protein <- "SPINT2"
protein <- "SH2D1A" # (?)
protein <- "BCL2L11" #?
protein <- "TNFRSF4" # (?)
protein <- "CCL20" # (?)
protein <- "CD200R1"
protein <- "CLEC4G" # (?)
protein <- "FOXO1"
protein <- "ITM2A"
protein <- "NT5C3A" # (?)
protein <- "MAP2K6" # (?)
protein <- "PTX3" # (?)
protein <- "FASLG"
protein <- "CD79B"
protein <- "TGFA"
protein <- "NPPC"
protein <- "HLA-E"
protein <- "SELPLG"
protein <- "FCRL2"
protein <- "CD22"
protein <- "CLEC4C"
protein <- "PARP1"
protein <- "SIGLEC10"
protein <- "SLAMF7"
protein <- "ICAM4"
protein <- "FCAR"

cowplot::plot_grid(
  plotDat %>%
    ggboxplot(x = "H1N1_reclassify", y = protein,
              paletter = "jco", add = "jitter") + facet_wrap(~group)+
    stat_compare_means(comparisons = my_comparisons_v2.1, method = "t.test"),
  plotDat %>%
    ggboxplot(x = "H3N2_reclassify", y = protein,
              paletter = "jco", add = "jitter") + facet_wrap(~group)+
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  plotDat %>%
    ggboxplot(x = "B_reclassify", y = protein,
              paletter = "jco", add = "jitter") + facet_wrap(~group)+
    stat_compare_means(comparisons = my_comparisons_v2.1, method = "t.test"),
  nrow = 1
)

plotDat %>%
  ggboxplot(x = strain, y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~group)+
  stat_compare_means(comparisons = my_comparisons_v2.1, method = "t.test")


plotDat %>%
  ggboxplot(x = strain, y = protein,
            paletter = "jco", add = "jitter") + facet_wrap(~group)+
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

