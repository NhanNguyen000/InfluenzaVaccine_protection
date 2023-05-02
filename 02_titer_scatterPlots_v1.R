library(ggpubr)
library(readxl)

get.log2 <- function(dat) {
  library(tidyverse)
  library(rstatix)
  outcome <- dat %>% mutate(across(where(is.numeric), ~log2(.x))) %>%
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))
  return(outcome)
}


# load MN titer ----------------------------
# iMED data
MNtiter_raw <- read.csv2("/vol/projects/CIIM/Influenza/iMED/metadata/MNtiters_raw.csv") # raw data
MNtiter <- read.csv("/vol/projects/CIIM/Influenza/iMED/metadata/MNtiters.csv") %>% # log2
  rename("H1N1_abFC" = "H1N1", "H3N2_abFC" = "H3N2", "B_abFC" = "B")

names(MNtiter)[-1] <- paste0("MN_", names(MNtiter)[-1])
MN_iMED <- MNtiter %>% rename("patientID" = "X")

# ZirFlu 2019 
MN_rawDat2019 <- read_excel("../metadata/20230105_ZirFlu-2019-2020_MN_Titer_Valerie&Nhan.xlsx",
                            sheet = "2019-2020 MN Titer")
MN_titer_2019 <- MN_rawDat2019 %>% 
  select(1:14) %>% slice(-c(1, 2, 36, 52)) %>%
  fill(1, .direction = "down") %>%
  as.data.frame() %>%
  mutate_at(c(3:14), as.numeric) %>% get.log2()
names(MN_titer_2019) <- c("condition", "patientID", 
                          paste0("H1N1_", c("T1", "T2", "T3")),
                          paste0("H3N2_", c("T1", "T2", "T3")),
                          paste0("Bvictoria_", c("T1", "T2", "T3")),
                          paste0("Byamagata_", c("T1", "T2", "T3")))

MNabFC_rawDat2019 <- read_excel("../metadata/20230105_ZirFlu-2019-2020_MN_Titer_Valerie&Nhan.xlsx",
                                sheet = "2019-2020 MN Fold-increase")
MN_abFC2019 <- MNabFC_rawDat2019 %>% 
  slice(-c(1, 2, 36, 52)) %>%
  fill(1, .direction = "down") %>%
  as.data.frame() %>%
  mutate_at(c(3:6), as.numeric) %>% get.log2()

names(MN_abFC2019) <- c("condition", "patientID",
                        paste0(c("H1N1_", "H3N2_", "Bvictoria_", "Byamagata_"), "abFC"))

MN_2019 <- MN_titer_2019 %>% full_join(MN_abFC2019) %>% 
  mutate(season = "2019",
         condition = str_to_lower(condition)) %>%
  relocate(season)
names(MN_2019)[4:19] <- paste0("MN_", names(MN_2019)[4:19])

# merge MN titer together
MN_metadat <- MN_iMED %>% full_join(MN_2019) %>% 
  left_join(metadat)

# make scatter plot ----------------------
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("young", "old","cirrhosis")))

metadat %>%
  ggplot(aes(x = H1N1_abFC_combine, y = H1N1_T1_combine)) + 
  geom_point(size = 1.5, position = position_jitter(width = 0.25, height = 0.25)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 7) +
  theme_classic() + ggtitle("HAI titer of the H1N1 strain")

MN_metadat %>%
  ggplot(aes(x = MN_H1N1_abFC, y = MN_H1N1_T1)) + 
  geom_point(size = 1.5, position = position_jitter(width = 0.25, height = 0.25)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 5) +
  theme_classic() + ggtitle("MN titer of the H1N1 strain")

MN_metadat %>%
  ggplot(aes(x = H1N1_T1_combine, y = MN_H1N1_T1)) + 
  geom_point(size = 1.5, position = position_jitter(width = 0.25, height = 0.25)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 5) + xlim(0, 10) + ylim(0, 10) + 
  theme_classic() + ggtitle("MN vs HAI baseline titer of the H1N1 strain")

MN_metadat %>%
  ggplot(aes(x = H1N1_abFC_combine, y = MN_H1N1_abFC)) + 
  geom_point(size = 1.5, position = position_jitter(width = 0.25, height = 0.25)) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 5) + xlim(0, 10) + 
  theme_classic() + ggtitle("MN vs HAI abFC titer of the H1N1 strain")

# previous code --------------------------

metadat %>%
  ggplot(aes(x = H1N1_T1_combine, y = H1N1_abFC_combine, color = group)) + 
  geom_point(size = 2, 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.25)) +
  # geom_jitter(size = 2) +
  # geom_smooth(method = lm, se = FALSE) + 
  # stat_cor(label.x = 8) + 
  theme_classic() + theme(legend.position = "top") +
  geom_hline(yintercept = log2(4), linetype = "dashed") +
  geom_vline(xintercept = log2(30), linetype = "dashed")

metadat %>%
  ggplot(aes(x = H1N1_T1_combine, y = H1N1_abFC_combine, color = group)) + 
  geom_point(size = 2, 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.25)) +
  # geom_jitter(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 7.5) +
  theme_classic() + theme(legend.position = "top")

metadat %>%
  ggplot(aes(x = H1N1_T1_combine, y = H1N1_abFC_combine, color = H1N1_reclassify)) + 
  geom_point(size = 2, 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.25)) +
  # geom_jitter(size = 2) +
  # geom_smooth(method = lm, se = FALSE) + 
  # stat_cor(label.x = 8) + 
  theme_classic() + theme(legend.position = "top") +
  geom_hline(yintercept = log2(4), linetype = "dashed") +
  geom_vline(xintercept = log2(40), linetype = "dashed")

metadat %>%
  ggplot(aes(x = H3N2_T1_combine, y = H3N2_abFC_combine, color = group)) + 
  geom_point(size = 2, 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.25)) +
  # geom_jitter(size = 2) +
  # geom_smooth(method = lm, se = FALSE) + 
  # stat_cor(label.x = 8) + 
  theme_classic() + theme(legend.position = "top") +
  geom_hline(yintercept = log2(4), linetype = "dashed") +
  geom_vline(xintercept = log2(30), linetype = "dashed")

metadat %>%
  ggplot(aes(x = H3N2_T1_combine, y = H3N2_abFC_combine, color = group)) + 
  geom_point(size = 2, 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.25)) +
  # geom_jitter(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(label.x = 7) +
  theme_classic() + theme(legend.position = "top")

# boxplot for 4 reclassificatioin --------------------
plot_dat <- metadat %>% 
  select(patientID, H1N1_T1_combine, H1N1_reclassify,
         H1N1_T2, iMED_H1N1_T4) %>%
  mutate(H1N1_T2_combine = ifelse(is.na(H1N1_T2), iMED_H1N1_T4, H1N1_T2)) %>% 
  select(-H1N1_T2, -iMED_H1N1_T4) %>% 
  gather(matches("combine"), key = time, value = HAItiter) %>%
  mutate(time = str_sub(time, 1, 7)) %>% 
  rename("reclassify" = "H1N1_reclassify")

p <- plot_dat %>% 
  ggplot(aes(time, HAItiter)) + geom_boxplot() + geom_point() +
  geom_line(aes(group = patientID)) + theme_bw() +
  facet_wrap(vars(reclassify), nrow = 1)

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("plum", "lightskyblue", "lightskyblue", "lightskyblue")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
library(grid)
grid.draw(g)

