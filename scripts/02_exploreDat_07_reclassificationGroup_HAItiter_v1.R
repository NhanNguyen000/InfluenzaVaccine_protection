rm(list = ls())
library(tidyverse)
library(grid)

load("cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))


# Classification ------------------------------------------------------
# iMED season 2015 - H1N1
H1N1_HAI_2015 <- metadata_healthy %>% filter(season == "2015") %>%
  dplyr::select(probandID, cohort, group, 
                H1N1_T1_log2, H1N1_T4_log2, H1N1_abFC, H1N1_reclassify) %>%
  mutate(category = ifelse(H1N1_abFC >= 4, "Response", "Non-response")) %>%
  gather(matches("_T"), key = time, value = HAItiter) %>%
  rename_at(vars(ends_with("_reclassify")), ~ "reclassify") %>%
  mutate(time = gsub("_log2", "", time),)

## Current classification (responder vs. non-responder) ------------------
plotDat_NRvsR <- H1N1_HAI_2015  %>%
  ggplot(aes(time, HAItiter)) + geom_boxplot() + geom_point() +
  geom_line(aes(group = probandID)) + 
  facet_wrap(vars(category), nrow = 1) +
  theme_bw() +
  ylab("Log2 of HAI titer - H1N1")

plotDat_NRvsR

# change the strip.background for facet_warp
plotDat_NRvsR_v2 <- ggplot_gtable(ggplot_build(plotDat_NRvsR))
strip_both <- which(grepl('strip-', plotDat_NRvsR_v2$layout$name))
fills <- c("#696969CC","#3194CCCC")

k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', plotDat_NRvsR_v2$grobs[[i]]$grobs[[1]]$childrenOrder))
  plotDat_NRvsR_v2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(plotDat_NRvsR_v2)

# re-classification (responder vs. non-responder) ------------------
plotDat_4groups <- H1N1_HAI_2015  %>%
  ggplot(aes(time, HAItiter)) + geom_boxplot() + geom_point() +
  geom_line(aes(group = probandID)) + 
  facet_wrap(vars(reclassify), nrow = 1) +
  theme_bw() +
  ylab("Log2 of HAI titer - H1N1")

plotDat_4groups 

# change the strip.background for facet_warp
plotDat_4groups_v2 <- ggplot_gtable(ggplot_build(plotDat_4groups))
strip_both <- which(grepl('strip-', plotDat_4groups_v2$layout$name))
fills <- c("#AF76A4","#5AA9D6", "#5AA9D6", "#3194CCCC")

k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', plotDat_4groups_v2$grobs[[i]]$grobs[[1]]$childrenOrder))
  plotDat_4groups_v2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(plotDat_4groups_v2)

# can test with y-axis order instead of numeric value --------------
a <- plotDat %>%
  mutate(group = ifelse(HAItiter > 10, as.character(HAItiter), "10")) %>%
  mutate(group = factor(group, levels = c("10", "20", "40", "80", "160", "320", "640", "1280", "2560", "5120", "10240")))
