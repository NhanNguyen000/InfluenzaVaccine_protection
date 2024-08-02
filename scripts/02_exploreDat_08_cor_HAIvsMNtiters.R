rm(list = ls())

library(tidyverse)
library(ggpubr)

# load data ------------------------------------------------------------------------------
load("processedDat/cohorts_dat.RData")

dat_titer <- list()

dat_titer$HAI <- cohorts$HAI_all %>%
  select(probandID, season, matches("_d0|_d28|_abFC")) %>%
  filter(season != "2020") %>%
  pivot_longer(cols = -c(1, 2), names_to = "strain_titer", values_to = "HAI") %>% drop_na()

dat_titer$MN <- cohorts$MN_all %>%
  select(probandID, season, matches("_d0|_d28|_abFC"))  %>%
  pivot_longer(cols = -c(1, 2), names_to = "strain_titer", values_to = "MN") %>% drop_na()

# plot correlation between HAI and MN ------------------------------------------------------
plotDat <- dat_titer %>% purrr::reduce(full_join)

set.seed(123) # make jitter plot reproducible

plotList <- list()
for (strainTiter_name in  unique(plotDat$strain_titer)) {
  plotDat_temp <- plotDat %>%  filter(strain_titer == strainTiter_name)
  
  plotList[[strainTiter_name]] <- plotDat_temp %>% 
    ggplot(aes(x = HAI, y = MN, color = season)) + 
    geom_jitter( size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", se=FALSE) + 
    stat_cor(aes(color = season), size = 8, method =  "spearman") + 
    theme_classic()+
    theme(text = element_text(size = 20),
          legend.position="top") +
    ggtitle(strainTiter_name)

}

## save the plot 
png("output/cor_HAIvsMN_titers.png", width = 1440, height = 1680)

cowplot::plot_grid(
  plotList$H1N1_d0_log2, plotList$H1N1_d28_log2, plotList$H1N1_abFC_log2,
  plotList$H3N2_d0_log2, plotList$H3N2_d28_log2, plotList$H3N2_abFC_log2,
  plotList$B_d0_log2, plotList$B_d28_log2, plotList$B_abFC_log2,
  plotList$Bvictoria_d0_log2, plotList$Bvictoria_d28_log2, plotList$Bvictoria_abFC_log2,
  plotList$Byamagata_d0_log2, plotList$Byamagata_d28_log2, plotList$Byamagata_abFC_log2,
  nrow = 5
)

dev.off()
