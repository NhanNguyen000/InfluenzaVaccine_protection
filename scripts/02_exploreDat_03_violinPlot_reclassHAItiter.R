rm(list = ls())
library(tidyverse)
library(grid)

load("processedDat/cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))


# prepare the data ------------------------------------------------------
H1N1_2015 <- metadata_healthy %>% filter(season == "2015") %>%
  dplyr::select(probandID, cohort, group, 
                H1N1_d0_log2, H1N1_d28_log2, H1N1_abFC, H1N1_reclassify) %>%
  mutate(category = ifelse(H1N1_abFC >= 4, "Responder", "Non-responder")) %>%
  gather(matches("_d"), key = time, value = HAItiter) %>%
  rename_at(vars(ends_with("_reclassify")), ~ "reclassify") %>%
  mutate(time = gsub("_log2", "", time),
         time = gsub("H1N1_", "", time))


# re-classification (responder vs. non-responder) ------------------
plot_reClass <- H1N1_2015  %>%
  ggplot(aes(time, HAItiter)) + 
  geom_violin(width = 0.9) +
  geom_boxplot(width = 0.2, color="gray41") + 
  geom_point(size = 3) +
  geom_line(aes(group = probandID), color="grey70") + 
  facet_wrap(vars(reclassify), nrow = 1) +
  theme_classic() + 
  theme(legend.position = "top", text = element_text(size = 24)) +
  ylab("Log2 (HAI titer of H1N1)")

plot_reClass 

# add color to the strip.background for facet_warp
plot_reClass_addColor <- ggplot_gtable(ggplot_build(plot_reClass))
strip_both <- which(grepl('strip-', plot_reClass_addColor$layout$name))

fills <- c("#9F5590","#5AA9D6", "#5AA9D6", "#3194CCCC") # follow the color in the scheme of re-calssification
#fills <- c("#751915","#65771E", "#B65008", "#036879") # follow colors in other plots in R

k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', plot_reClass_addColor$grobs[[i]]$grobs[[1]]$childrenOrder))
  plot_reClass_addColor$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
 # plot_reClass_addColor$grobs[[i]]$grobs[[1]]$children[[j]]$gp$alpha <- 0.7 # incase, use 4 colors in other plots in R in the "fills" object
  k <- k+1
}

# violin plot to demonstrate the relcassification concept
png("output/violinPlot_HAIabTiter_reClass.png", width = 720, height = 432)
grid.draw(plot_reClass_addColor)
dev.off()

# Current classification (responder vs. non-responder) ------------------
plot_NRvsR <- H1N1_2015  %>%
  ggplot(aes(time, HAItiter)) + 
  geom_violin(width = 0.9) +
  geom_boxplot(width = 0.2, color="gray41") + 
  geom_point(size = 3) +
  geom_line(aes(group = probandID), color="grey70") + 
  facet_wrap(vars(category), nrow = 1) +
  theme_classic() + 
  theme(legend.position = "top", text = element_text(size = 24)) +
  ylab("Log2 (HAI titer of H1N1)")

plot_NRvsR

# add color to the strip.background for facet_warp
plot_NRvsR_addColor <- ggplot_gtable(ggplot_build(plot_NRvsR))
strip_both <- which(grepl('strip-', plot_NRvsR_addColor$layout$name))

fills <- c("#696969CC","#3194CCCC") # follow the color in the scheme of re-calssification
 
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', plot_NRvsR_addColor$grobs[[i]]$grobs[[1]]$childrenOrder))
  plot_NRvsR_addColor$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

# violin plot to demonstrate the current NR vs.R concept
png("output/violinPlot_HAIabTiter_NRvsR.png", height = 432)
grid.draw(plot_NRvsR_addColor)
dev.off()

