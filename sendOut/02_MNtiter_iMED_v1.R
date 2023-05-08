# load data -----------------------------------------
MNtiter_raw <- read.csv2("/vol/projects/CIIM/Influenza/iMED/metadata/MNtiters_raw.csv") # raw data
MNtiter <- read.csv("/vol/projects/CIIM/Influenza/iMED/metadata/MNtiters.csv") %>% # log2
  rename("ab_H1N1" = "H1N1", "ab_H3N2" = "H3N2", "ab_B" = "B")

names(MNtiter)[-1] <- paste0("MN_", names(MNtiter)[-1])
dat_plot <- MNtiter %>% rename("patientID" = "X") %>% 
  inner_join(metadat, by = c("patientID"))

cols <- c("LL" = "red", "LH" = "blue", "HL" = "darkgreen", "HH" = "orange")

# abFC -----------------------------
dat_plot %>% 
  ggplot(aes(x = ab_H1N1, y = MN_ab_H1N1)) +
  geom_point(aes(col = H1N1_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

dat_plot %>% 
  ggplot(aes(x = ab_H3N2, y = MN_ab_H3N2)) +
  geom_point(aes(col = H3N2_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

dat_plot %>% 
  ggplot(aes(x = ab_B, y = MN_ab_B)) +
  geom_point(aes(col = B_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

# titer at baseline (T1) --------------------
dat_plot %>% 
  ggplot(aes(x = iMED_H1N1_T1, y = MN_H1N1_T1)) +
  geom_point(aes(col = H1N1_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

dat_plot %>% 
  ggplot(aes(x = iMED_H3N2_T1, y = MN_H3N2_T1)) +
  geom_point(aes(col = H3N2_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

dat_plot %>% 
  ggplot(aes(x = iMED_B_T1, y = MN_B_T1)) +
  geom_point(aes(col = B_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)
