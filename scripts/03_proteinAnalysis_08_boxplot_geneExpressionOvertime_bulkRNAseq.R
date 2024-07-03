rm(list = ls())

library(tidyverse)

# load data =======================================================================

## load data from iMED transcriptome, iMED cohort, season 2015 ---------------
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
# View(transcriptomeList[[1]])
# View(transcriptomeList[[2]])
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1") # time T1 in trancriptome = d0
iMED_transcripDat <- transcriptomeList[[1]] %>% as.data.frame %>% 
  select(iMED_transcrip_T1$SampleName)


## metadata for all healthy subjects ----------------------------------------
load("processedDat/cohorts_dat.RData")
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# boxplot for selected proteins  overtime ------------------------------------------
## prepare data --------------------------------------------------------------------
iMED_transcrip_overtime <- transcriptomeList[[2]] %>%
  filter(SampleTime %in% c("T1", "T4")) # time T1 in trancriptome = d0, T4 = d28

iMED_transcripDat <- transcriptomeList[[1]] %>% as.data.frame %>% 
  select(iMED_transcrip_overtime$SampleName)

metadat_iMED_2015_overtime <- metadata_healthy %>% select(-time) %>%
  filter(season == "2015", cohort == "iMED") %>%
  right_join(iMED_transcrip_overtime, by = c("probandID" = "patientID")) %>% 
  rename("time" = "SampleTime")


## boxplot of gene expression overtime ---------------------------------------------------
protein <- "CD83" # selected protein 

plotDat_boxplot_overtime <- metadat_iMED_2015_overtime %>% 
  full_join(iMED_transcripDat %>% t() %>% 
              as.data.frame %>% select(protein) %>% rownames_to_column("SampleName")) 

# plot fo H1N1
stat.test <- plotDat_boxplot_overtime %>% 
  group_by(season, H1N1_reclassify) %>%  
  t_test(CD83 ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(H1N1_reclassify == "LL", 0.9, 
                       ifelse(H1N1_reclassify == "LH", 1.9,
                              ifelse(H1N1_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(H1N1_reclassify == "LL", 1.1, 
                       ifelse(H1N1_reclassify == "LH", 2.1,
                              ifelse(H1N1_reclassify == "HL", 3.1, 4.1))))

bxp_H1N1 <- ggboxplot(plotDat_boxplot_overtime, 
                      x = "H1N1_reclassify", y = protein, color = "time",
                      add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  scale_color_manual(values=c( c("#B65008", "#034E91"))) +
  theme(text = element_text(size = 18)) +
  stat_pvalue_manual(stat.test, label = "p.signif", size = 5)

# plot fo H3N2
stat.test <- plotDat_boxplot_overtime %>% 
  group_by(season, H3N2_reclassify) %>%  
  t_test(CD83 ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(H3N2_reclassify == "LL", 0.9, 
                       ifelse(H3N2_reclassify == "LH", 1.9,
                              ifelse(H3N2_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(H3N2_reclassify == "LL", 1.1, 
                       ifelse(H3N2_reclassify == "LH", 2.1,
                              ifelse(H3N2_reclassify == "HL", 3.1, 4.1))))

bxp_H3N2 <- ggboxplot(plotDat_boxplot_overtime,
                      x = "H3N2_reclassify", y = protein, color = "time",
                      add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  scale_color_manual(values=c( c("#B65008", "#034E91"))) +
  theme(text = element_text(size = 18)) +
  stat_pvalue_manual(stat.test, label = "p.signif", size = 5)

# plot fo B
stat.test <- plotDat_boxplot_overtime %>% 
  group_by(season, B_reclassify) %>%  
  t_test(CD83 ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(B_reclassify == "LL", 0.9, 
                       ifelse(B_reclassify == "LH", 1.9,
                              ifelse(B_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(B_reclassify == "LL", 1.1, 
                       ifelse(B_reclassify == "LH", 2.1,
                              ifelse(B_reclassify == "HL", 3.1, 4.1))))

bxp_B <- ggboxplot(plotDat_boxplot_overtime, 
                   x = "B_reclassify", y = protein, color = "time",
                   add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  scale_color_manual(values=c( c("#B65008", "#034E91"))) +
  theme(text = element_text(size = 18)) +
  stat_pvalue_manual(stat.test, label = "p.signif", size = 5)

# all plot together
cowplot::plot_grid(bxp_H1N1, bxp_H3N2, bxp_B, nrow = 1)

## save the plot 
png("output/boxplot_RNAseq_CD83_expressionOvertime.png", width = 816, height = 432)
cowplot::plot_grid(bxp_H1N1, bxp_H3N2, bxp_B, nrow = 1)
dev.off()

