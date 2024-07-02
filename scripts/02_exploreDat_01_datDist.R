rm(list = ls())
library(tidyverse)
library(ggpubr)

load("processedDat/cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season),
         sex = ifelse(sex == "f", "female", ifelse(sex == "m", "male", sex))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# cohort size, age, sex -------------------------
cohorts_dat <- list()
cohorts_dat$iMED_2014 <- metadata_healthy %>% filter(group == "iMED_2014")
cohorts_dat$iMED_2015 <- metadata_healthy %>% filter(group == "iMED_2015")
cohorts_dat$ZirFlu_2019 <- metadata_healthy %>% filter(group == "ZirFlu_2019")
cohorts_dat$ZirFlu_2020 <- metadata_healthy %>% filter(group == "ZirFlu_2020")

range(cohorts_dat$iMED_2014$age)
range(cohorts_dat$iMED_2015$age)
range(cohorts_dat$ZirFlu_2019$age)
range(cohorts_dat$ZirFlu_2020$age)

cohorts_dat$iMED_2014 %>% count(sex)
cohorts_dat$iMED_2015 %>% count(sex)
cohorts_dat$ZirFlu_2019 %>% count(sex)
cohorts_dat$ZirFlu_2020 %>% count(sex)

# sex with age, antibody titer at baseline (d0), and antibody foldchange (abFC) distribution --------------------------

## age and sex distribution per each season --------------------------------------
plot_ageSexDist <- metadata_healthy %>%
  ggboxplot(x = "sex", y = "age", color = "sex", 
            add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  scale_color_manual(values=c( c("#9B2727", "#034E91"))) +
  facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = list(c("female", "male")), size = 5, method = "t.test") +
  theme(text = element_text(size = 16))

png("output/boxplot_ageSexDistribution.png")
plot_ageSexDist 
dev.off()

## distribution of antibody titer at basline and abFC per sex in each vaccine season --------------------------
plotList <- c()

# for H1N1 and H3N2
overlaped_HAImeasures <- c("H1N1_d0_log2", "H1N1_abFC_log2", "H3N2_d0_log2", "H3N2_abFC_log2")
for (HAImeasure in overlaped_HAImeasures) {
  plotList[[HAImeasure]] <- metadata_healthy %>%
    ggboxplot(x = "sex", y = HAImeasure, color = "sex", 
              add = "jitter", add.params = list(size = 3, alpha = 0.5), xlab = FALSE) + 
    scale_color_manual(values=c( c("#9B2727", "#034E91"))) +
    facet_wrap(~season, nrow = 1) + 
    stat_compare_means(comparisons = list(c("female", "male")), size = 5, method = "t.test") +
    theme(text = element_text(size = 16))
}

# for B strains
Bstrain_HAImeasures <- c("B_d0_log2", "B_abFC_log2", 
                         "Bvictoria_d0_log2", "Bvictoria_abFC_log2", 
                         "Byamagata_d0_log2", "Byamagata_abFC_log2")

for (HAImeasure in Bstrain_HAImeasures) {
  plotList[[HAImeasure]] <- metadata_healthy %>% drop_na(all_of(HAImeasure)) %>% 
    ggboxplot(x = "sex", y = HAImeasure, color = "sex", 
              add = "jitter", add.params = list(size = 3, alpha = 0.5), xlab = FALSE) + 
    scale_color_manual(values=c( c("#9B2727", "#034E91"))) +
    facet_wrap(~season, nrow = 1) + 
    stat_compare_means(comparisons = list(c("female", "male")), method = "t.test") +
    theme(text = element_text(size = 16))
}

# arrange the plots into 1 graph
arrangedPlots <- cowplot::plot_grid(
  cowplot::plot_grid(
    cowplot::plot_grid(plotList$H1N1_d0_log2 + theme(legend.position="none"), 
                       plotList$H1N1_abFC_log2 + theme(legend.position="none"), 
                       nrow = 2),
    cowplot::plot_grid(plotList$H3N2_d0_log2 + theme(legend.position="none"), 
                       plotList$H3N2_abFC_log2 + theme(legend.position="none"), 
                       nrow = 2)
  ),
  cowplot::plot_grid(
    cowplot::plot_grid(plotList$B_d0_log2 + theme(legend.position="none"), 
                       plotList$B_abFC_log2 + theme(legend.position="none"), 
                       nrow = 2),
    cowplot::plot_grid(plotList$Bvictoria_d0_log2 + theme(legend.position="none"),
                       plotList$Bvictoria_abFC_log2 + theme(legend.position="none"),
                       nrow = 2),
    cowplot::plot_grid(plotList$Byamagata_d0_log2 + theme(legend.position="none"), 
                       plotList$Byamagata_abFC_log2 + theme(legend.position="none"), 
                       nrow = 2),
    nrow = 1
  ),
  nrow = 2
)

arrangedPlots

# add legend color for the whole plot
legend_b <- get_legend(
  plotList$H1N1_d0_log2 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top")
)

cowplot::plot_grid(arrangedPlots, legend_b, ncol = 1, rel_heights = c(1, .1))

# save the plot 
png("output/boxplot_ageAbTiter_abFC_distribution.png", width = 1200, height = 1350)
arrangedPlots
#cowplot::plot_grid(arrangedPlots, legend_b, ncol = 1, rel_heights = c(1, .1)) # add legend color
dev.off()

