rm(list = ls())

library(tidyverse)
library(ggpubr)
library(rstatix)
# load data =======================================================================
load("processedDat/cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- mebo_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

# boxplot ======================================================================
# select season and metabolite
strainSeasons <- c("H1N1_2014", "H1N1_2015", "H1N1_2019", "H1N1_2020",
                   "B_2014", "B_2015", "H3N2_2015", "Byamagata_2020") # season have at least >= 3 participant per reclassification group

mebo <- "C3H7NO2S"

metadat_boxplot <- inputDat %>% 
  select(probandID, season, responder, time, c(mebo),
         matches("_abFC|_d0|_d28|_reclassify")) %>% filter(time %in% c("d0", "d28")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season))

## boxplots for each timepoint, per group, per season --------------------------
bxp <- metadat_boxplot %>% 
  filter(strainSeason %in% strainSeasons) %>%
  ggboxplot(x = "reclassify", y = mebo, color = "time",
            add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  scale_color_manual(values=c( c("#B65008", "#034E91"))) +
  facet_wrap(~strainSeason, nrow = 2)+
  theme(text = element_text(size = 18))

bxp

# add non sig. paired t_test to the graph ---------------------------------------
metadat_boxplot2 <- metadat_boxplot %>% add_count(probandID, strainSeason) %>% filter(n == 2)

stat.test <- metadat_boxplot2 %>% 
  filter(strainSeason %in% strainSeasons) %>%
  group_by(strainSeason, reclassify) %>%  
  t_test(C3H7NO2S ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(reclassify == "LL", 0.9, 
                       ifelse(reclassify == "LH", 1.9,
                              ifelse(reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(reclassify == "LL", 1.1, 
                       ifelse(reclassify == "LH", 2.1,
                              ifelse(reclassify == "HL", 3.1, 4.1))))

bxp_stat <- bxp +
  stat_pvalue_manual(stat.test, label = "p.signif", size = 5) + 
  theme(legend.position = "top")

bxp_stat


# save the plot ---------------------------------------------------------------
png("output/boxplotMetabolite_reClass_C3H7NO2S_overTime.png", width = 768, height = 520)
bxp_stat
dev.off()


