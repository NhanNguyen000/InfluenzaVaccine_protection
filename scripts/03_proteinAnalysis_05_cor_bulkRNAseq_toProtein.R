rm(list = ls())

library(tidyverse)
library(ggpubr)

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

# prepare plot data ------------------------------------------
protein <- "CD83"

rnaDat <- metadat_iMED_2015 %>% 
  full_join(iMED_transcripDat %>% t() %>% as.data.frame %>%
              rownames_to_column("SampleName")) %>% 
  dplyr::select(probandID, protein)
names(rnaDat)[2] <- paste0(names(rnaDat)[2], "_geneExpression")

proteinDat <- protein_Dat$iMED_2015 %>% 
  rownames_to_column("name") %>% right_join(metadata_healthy)  %>% 
  dplyr::select(probandID, protein)
names(proteinDat)[2] <- paste0(names(proteinDat)[2], "_proteinAbundance")

plotDat <- rnaDat %>% left_join(proteinDat)

# scatter plots -----------------------------------------------------------------
plotDat %>%
  ggplot(aes(x = CD83_geneExpression, y = CD83_proteinAbundance)) +
  geom_point(size =5, alpha = 0.5) +
  geom_smooth(method = "lm", se=FALSE) + stat_cor(size = 8) + 
  theme_classic()+
  theme(text = element_text(size = 24))
