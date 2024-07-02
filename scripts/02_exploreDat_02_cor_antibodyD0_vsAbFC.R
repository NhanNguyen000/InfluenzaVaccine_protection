rm(list = ls())
library(tidyverse)
library(ggpubr)
library(webr)
library(patchwork)

load("processedDat/cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# antibody T1 and abFC negative correlation ------------------ ------------------
## prepare data ------------------ ------------------
plotDat <- metadata_healthy %>% 
  mutate(H1N1_d0 = ifelse(H1N1_d0 > 10, as.character(H1N1_d0), "10")) %>% # Convert to make the plot look better
  mutate(H3N2_d0 = ifelse(H3N2_d0 > 10, as.character(H3N2_d0), "10")) %>% 
  mutate(B_d0 = ifelse(B_d0 > 10, as.character(B_d0), "10")) %>% 
  mutate(Bvictoria_d0 = ifelse(Bvictoria_d0 > 10, as.character(Bvictoria_d0), "10")) %>% 
  mutate(Byamagata_d0 = ifelse(Byamagata_d0 > 10, as.character(Byamagata_d0), "10")) %>% 
  mutate_at(vars(ends_with("_d0")), ~factor(.x, levels = c("10", "20", "40", "80", "160", "320", "640", "1280", "2560")))


## make plot ------------------ ------------------
plotList <- list()
plotList$H1N1 <- plotDat %>% #filter(season == "2015") %>%
  ggplot() +
  geom_boxplot(aes(x = H1N1_d0, y = H1N1_abFC_log2), outlier.color = NA) + 
  geom_jitter(aes(x = H1N1_d0, y = H1N1_abFC_log2, color = H1N1_reclassify), 
              width = 0.2, size = 3, alpha = 0.8) + 
  scale_color_manual(values=c( c("#792C74", "#65771E", "#B65008", "#036879"))) +
  geom_smooth(aes(x = H1N1_d0_log2, y = H1N1_abFC_log2), method = "lm", se = FALSE) + # make regress line based on the log2(exact data value)
  stat_cor(aes(x = H1N1_d0_log2, y = H1N1_abFC_log2), # calculate the correlation based on the log2(exact data value)
           label.x = 5, label.y = 10, size = 5) + 
  ylim(0, 12.5) + theme_classic() + 
  theme(legend.position = "top", text = element_text(size = 16))

plotList$H3N2 <- plotDat %>% 
  ggplot() +
  geom_boxplot(aes(x = H3N2_d0, y = H3N2_abFC_log2), outlier.color = NA) + 
  geom_jitter(aes(x = H3N2_d0, y = H3N2_abFC_log2, color = H3N2_reclassify), 
              width = 0.2, size = 3, alpha = 0.8) + 
  scale_color_manual(values=c( c("#792C74", "#65771E", "#B65008", "#036879"))) +
  geom_smooth(aes(x = H3N2_d0_log2, y = H3N2_abFC_log2), method = "lm", se = FALSE) + # make regress line based on the log2(exact data value)
  stat_cor(aes(x = H3N2_d0_log2, y = H3N2_abFC_log2), # calculate the correlation based on the log2(exact data value)
           label.x = 5, label.y = 10, size = 5) + 
  ylim(0, 12.5) + theme_classic() + 
  theme(legend.position = "top", text = element_text(size = 16))

plotList$B <- plotDat %>% drop_na(B_d0) %>%
  ggplot() +
  geom_boxplot(aes(x = B_d0, y = B_abFC_log2), outlier.color = NA) + 
  geom_jitter(aes(x = B_d0, y = B_abFC_log2, color = B_reclassify), 
              width = 0.2, size = 3, alpha = 0.8) + 
  scale_color_manual(values=c( c("#792C74", "#65771E", "#B65008", "#036879"))) +
  geom_smooth(aes(x = B_d0_log2, y = B_abFC_log2), method = "lm", se = FALSE) + # make regress line based on the log2(exact data value)
  stat_cor(aes(x = B_d0_log2, y = B_abFC_log2), # calculate the correlation based on the log2(exact data value)
           label.x = 5, label.y = 10, size = 5) + 
  ylim(0, 12.5) + theme_classic() + 
  theme(legend.position = "top", text = element_text(size = 16))

plotList$Bvictoria <- plotDat %>% drop_na(Bvictoria_d0) %>%
  ggplot() +
  geom_boxplot(aes(x = Bvictoria_d0, y = Bvictoria_abFC_log2), outlier.color = NA) + 
  geom_jitter(aes(x = Bvictoria_d0, y = Bvictoria_abFC_log2, color = Bvictoria_reclassify), 
              width = 0.2, size = 3, alpha = 0.8) + 
  scale_color_manual(values=c( c("#792C74", "#65771E", "#B65008", "#036879"))) +
  geom_smooth(aes(x = Bvictoria_d0_log2, y = Bvictoria_abFC_log2), method = "lm", se = FALSE) + # make regress line based on the log2(exact data value)
  stat_cor(aes(x = Bvictoria_d0_log2, y = Bvictoria_abFC_log2), # calculate the correlation based on the log2(exact data value)
           label.x = 5, label.y = 10, size = 5) + 
  ylim(0, 12.5) + theme_classic() + 
  theme(legend.position = "top", text = element_text(size = 16))

plotList$Byamagata <- plotDat %>% drop_na(Byamagata_d0) %>%
  ggplot() +
  geom_boxplot(aes(x = Byamagata_d0, y = Byamagata_abFC_log2), outlier.color = NA) + 
  geom_jitter(aes(x = Byamagata_d0, y = Byamagata_abFC_log2, color = Byamagata_reclassify), 
              width = 0.2, size = 3, alpha = 0.8) + 
  scale_color_manual(values=c( c("#792C74", "#65771E", "#B65008", "#036879"))) +
  geom_smooth(aes(x = Byamagata_d0_log2, y = Byamagata_abFC_log2), method = "lm", se = FALSE) + # make regress line based on the log2(exact data value)
  stat_cor(aes(x = Byamagata_d0_log2, y = Byamagata_abFC_log2), # calculate the correlation based on the log2(exact data value)
           label.x = 5, label.y = 10, size = 5) + 
  ylim(0, 12.5) + theme_classic() + 
  theme(legend.position = "top", text = element_text(size = 16))

## save the plot ------------------ ------------------
# for H1N1
ggsave("output/cor_antibodyD0_vsAbFC_H1N1.png", 
       plotList$H1N1, device = "png")

# for H3N2 and B strains
plotGraph <- cowplot::plot_grid(
  plotList$H3N2, plotList$B, plotList$Bvictoria, plotList$Byamagata,
  nrow = 2
)
ggsave("output/cor_antibodyD0_vsAbFC_H3N2_Bstrains.png", 
       plotGraph, device = "png")

