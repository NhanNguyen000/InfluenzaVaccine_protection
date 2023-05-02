# check groups
cohorts_dat$donorInfo_all %>% count(cohort, season, disease, age_group)

# age distribution
plot_dat <- cohorts_dat$donorInfo_all %>%
  mutate(group = paste(cohort, season, disease, sep = "_"))
ggplot(plot_dat, aes(x=age, fill=group)) + geom_density(alpha=.3) + xlim(20, 81)
ggplot(plot_dat %>% filter(season %in% c(NA, "2019")), 
       aes(x=age, fill=group)) + geom_density(alpha=.3)

ggplot(plot_dat, aes(x=age, fill=cohort)) + geom_density(alpha=.3)+ xlim(20, 81)

plot_dat %>% filter(disease == "healthy") %>% 
  ggplot(aes(x=age, fill=cohort)) + geom_density(alpha=.3) + xlim(20, 81)

plot_dat %>% 
  mutate(cohort2 = ifelse(cohort == "iMED", cohort, paste0(cohort, "_", disease))) %>% 
  ggplot(aes(x=age, fill=cohort2)) + geom_density(alpha=.3) + xlim(20, 81)

plot_dat %>% 
  mutate(cohort2 = ifelse(cohort == "iMED", cohort, paste0(cohort, "_", disease))) %>% 
  filter(cohort2 != "ZirFlu_cirrhosis") %>%
  ggplot(aes(x=age, fill=cohort2)) + geom_density(alpha=.3) + xlim(20, 81)
