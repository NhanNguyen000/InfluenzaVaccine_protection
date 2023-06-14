rm(list = ls())
library(tidyverse)
library(ggpubr)

load("cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
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

# sex with age, antibody titer, and antibody abFC distribution --------------------------
# option 1 for boxplot
metadata_healthy %>%
  ggboxplot(x = "season", y = "age", color = "sex",
            paletter = "jco", add = "jitter") + 
  stat_compare_means(aes(group = sex), method = "t.test")

# option 2 for boxplot

metadata_healthy %>%
  ggboxplot(x = "sex", y = "age", color = "sex",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = list(c("f", "m")), method = "t.test")

metadata_healthy %>%
  ggboxplot(x = "sex", y = "H1N1_T1_log2", color = "sex",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = list(c("f", "m")), method = "t.test")

metadata_healthy %>%
  ggboxplot(x = "sex", y = "H1N1_abFC_log2", color = "sex",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = list(c("f", "m")), method = "t.test")




