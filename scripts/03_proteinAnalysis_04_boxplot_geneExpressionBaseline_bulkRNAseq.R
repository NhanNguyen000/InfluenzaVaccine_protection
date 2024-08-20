rm(list = ls())

library(tidyverse)

# load data =======================================================================

## load data from iMED transcriptome, iMED cohort, season 2015 ------------------------------
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
# View(transcriptomeList[[1]])
# View(transcriptomeList[[2]])
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1") # time T1 in trancriptome = d0
iMED_transcripDat <- transcriptomeList[[1]] %>% as.data.frame %>% 
  dplyr::select(iMED_transcrip_T1$SampleName) %>% as.data.frame()


## metadata for all healthy subjects -------------------------
load("processedDat/cohorts_dat.RData")
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              dplyr::select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

metadat_iMED_2015 <- metadata_healthy %>% filter(season == "2015", cohort == "iMED") %>%
  right_join(iMED_transcrip_T1, by = c("probandID" = "patientID")) %>% 
  arrange(factor(SampleName, levels = colnames(iMED_transcripDat)))

# boxplot for selected proteins at baseline------------------------------------------
## prepare plot data ------------------------------------------
protein <- "CD83"
protein <- "TMEM51"
protein <- "IFT52"
protein <- "ARSA"
protein <- "RIPOR1"
protein <- "CRYL1"
protein <- "TMEM204"
protein <- "ACVR2B"
protein <- "RNMT"

plotDat_boxplot <- metadat_iMED_2015 %>% 
  full_join(iMED_transcripDat %>% t() %>% 
              as.data.frame %>% dplyr::select(protein) %>% rownames_to_column("SampleName")) %>%
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_d0 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high")))

# check the numebr of paticipants per reclassification group
plotDat_boxplot %>% count(H1N1_reclassify)
plotDat_boxplot %>% count(H3N2_reclassify)
plotDat_boxplot %>% count(B_reclassify)

## make boxplot ------------------------------------------
compare_responder <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
compare_reClass_v2 <- list( c("LL", "LH"), c("LL", "HL"))
compare_reClass <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
compare_abFC <- list(c("R", "NR"))
compare_abD0 <- list(c("low", "high"))

# based on reclassification 
boxplot_reClass <- cowplot::plot_grid(
  plotDat_boxplot %>%
    ggboxplot(x = "H1N1_reclassify", y = protein,
              color = "#898366", add = "jitter", add.params = list(size = 4, alpha = 0.5)) + 
    stat_compare_means(comparisons = compare_reClass_v2, size = 6, method = "t.test",
                       label = "p.signif", hide.ns = TRUE, tip.length = 0, vjust = 0.5)+
    theme(text = element_text(size = 18)),
  plotDat_boxplot %>%
    ggboxplot(x = "H3N2_reclassify", y = protein,
              color = "#898366", add = "jitter", add.params = list(size = 4, alpha = 0.5)) + 
    stat_compare_means(comparisons = compare_reClass, size = 6, method = "t.test",
                       label = "p.signif", hide.ns = TRUE, tip.length = 0, vjust = 0.5)+
    theme(text = element_text(size = 18)),
  plotDat_boxplot %>%
    ggboxplot(x = "B_reclassify", y = protein,
              color = "#898366", add = "jitter", add.params = list(size = 4, alpha = 0.5)) + 
    stat_compare_means(comparisons = compare_reClass_v2, size = 6, method = "t.test",
                       label = "p.signif", hide.ns = TRUE, tip.length = 0, vjust = 0.5)+
    theme(text = element_text(size = 18)),
  nrow = 1
)

# based on previous vaccien classification concept
boxplot_abFC_abD0 <- cowplot::plot_grid(
  # based on abFC: NR vs R
  ggboxplot(plotDat_boxplot, x = "H1N1_abFC", y = protein,
            color = "#898366", add = "jitter", add.params = list(size = 4, alpha = 0.5)) + 
    stat_compare_means(comparisons = compare_abFC, size = 6, method = "t.test")+
    theme(text = element_text(size = 18)),
  
  # based on abT1: low Ab vs high Ab
  ggboxplot(plotDat_boxplot, x = "H1N1_abBaseline", y = protein,
            color = "#898366", add = "jitter", add.params = list(size = 4, alpha = 0.5)) + 
    stat_compare_means(comparisons = compare_abD0, size = 6, method = "t.test")+
    theme(text = element_text(size = 18))
)

# based on responders group: NR, other, TR
ggboxplot(plotDat_boxplot, x = "responder", y = protein,
          color = "#898366", add = "jitter", add.params = list(size = 4, alpha = 0.5)) + 
  stat_compare_means(comparisons = compare_responder, size = 6, method = "t.test")+
  theme(text = element_text(size = 18)) # only have NR and TR people for CD83

## save the plot --------------------------------------------------
#png(paste0("output/boxplotRNAseq_reClass_", protein, ".png"), width = 624, height = 480)
png(paste0("output/boxplotRNAseq_reClass_", protein, ".png"), width = 624, height = 384)
boxplot_reClass 
dev.off()

png(paste0("output/boxplotRNAseq_abFC_abD0_", protein, ".png"), width = 624, height = 480)
boxplot_abFC_abD0
dev.off()
