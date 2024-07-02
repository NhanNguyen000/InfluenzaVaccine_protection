rm(list = ls())
library(tidyverse)
library(grid)
library(ggpubr)

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
  mutate(category = ifelse(H1N1_abFC >= 4, "Response", "Non-response")) %>%
  select(-matches("log2|abFC")) %>%
  inner_join(cohorts$MN_all %>% select(probandID, H1N1_d0_log2, H1N1_d28_log2, H1N1_abFC)) %>%
  gather(matches("_log2"), key = time, value = MNtiter) %>%
  rename_at(vars(ends_with("_reclassify")), ~ "reclassify") %>%
  mutate(time = gsub("_log2", "", time))

# Current classification (responder vs. non-responder) ------------------------------------
plot_NRvsR <- H1N1_2015  %>%
  ggplot(aes(time, MNtiter)) + 
  geom_violin(width = 0.9) +
  geom_boxplot(width = 0.2, color="gray41") + 
  geom_point(size = 3) +
  geom_line(aes(group = probandID), color="grey70") + 
  facet_wrap(vars(category), nrow = 1) +
  theme_classic() + 
  theme(legend.position = "top", text = element_text(size = 24)) +
  ylab("Log2 (MN titer of H1N1)")

plot_NRvsR

# Boxplot compare T1 vs T4 per 2 classes: response vs. non-response ------------------------
H1N1_2015 %>%
  ggboxplot(x = "time", y = "MNtiter", facet.by = "category",
            paletter = "jco", add = "jitter") + 
  stat_compare_means(method = "t.test")

H1N1_2015 %>%
  ggboxplot(x = "category", y = "MNtiter", facet.by = "time",
            paletter = "jco", add = "jitter") + 
  stat_compare_means(method = "t.test") +
  ylab("Log2 of MN titer")

H1N1_2015 %>%
  ggboxplot(x = "category", y = "MNtiter", color = "time",
            paletter = "jco", add = "jitter") + 
  stat_compare_means(aes(group = time), method = "t.test")

# Boxplot compare d0 vs d28 per 4 classes: LL, LH, HL, HH ------------------------------------
H1N1_2015 %>%
  ggboxplot(x = "reclassify", y = "MNtiter", facet.by = "time", nrow = 1,
            paletter = "jco", add = "jitter") + 
  stat_compare_means(comparison = list(c("LL", "LH"), c("LL", "HL"), c("LL", "HH")),
                     method = "t.test") +
  ylab("Log2 of MN titer")
