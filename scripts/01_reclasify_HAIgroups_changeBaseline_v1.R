# load libraries & functions --------------------------
library(tidyverse)

get.reclassify <- function(dat, baseline_col, abFC_col, baseline_cutoff, abFC_cutoff) {
  # Description: reclassify vaccine response based on the baseline and abFC
  
  # Arguments: 
  # dat: data with patient ID, HAI value per time point, abFC of HAI between time point per column
  # baseline_col, abFC_col: names of the baseline and abFC columns
  # baseline_cutoff, abFC_cutoff: threshold for the classification
  
  # Returns: 
  # 4 reclassified groups: HH - high baseline & high abFC, 
  # HL - high baseline & log abFC, LH - low baseline & high abFC, LL - low baseline & low abFC
  
  reclassify <- paste0(ifelse(dat[, baseline_col] >= baseline_cutoff, "H", "L"),
                       ifelse(dat[, abFC_col] >= abFC_cutoff, "H", "L"))
  return(reclassify)
}

get.boxplot_reclassify_perStrain <- function(dat, strain_list) {
  # Description: boxplot for 4 reclassified vaccine response groups per strain
  
  # Arguments: 
  # dat: data with patient ID, HAI values, abFC of HAI, vaccine reclassify per column
  # strain_list: names of strain in the season
  # baseline_cutoff, abFC_cutoff: threshold for the classification
  
  # Returns: Boxplot of 4 reclassified groups per strain, with each dot is patients, 
  # and dot connection project the change of HAI level in patient 
  # 4 reclassified groups: HH - high baseline & high abFC, HL - high baseline & log abFC, 
  # LH - low baseline & high abFC, LL - low baseline & low abFC
  
  plot_list <- list()
  
  for (strain in strain_list) {
    plotDat <- dat %>% dplyr::select(patientID, matches(strain), -matches("ab")) %>%
      gather(matches("_T"), key = time, value = HAItiter) %>%
      arrange(patientID) %>%
      rename_at(vars(ends_with("_reclassify")), ~ "reclassify")
    
    plot_list[[strain]] <- plotDat %>%
      ggplot(aes(time, HAItiter)) + geom_boxplot() + geom_point() +
      geom_line(aes(group = patientID)) + theme_bw() +
      facet_wrap(vars(reclassify), nrow = 1)
  }
  
  return(plot_list)
}

# cutoff =======================================================================
# (baseline < 40 ) == T1_low. (abFC >= 4) == abFC_high
# (baseline >= 40) == T1_high. (abFC < 4) == abFC_low
baseline <- log2(40)
abFC <- log2(4)

# ZirFlu =====================================================================
ZirFlu_strains <- c("H1N1", "H3N2", "Bvictoria", "Byamagata")
reclassify_temp <- list()
for (strain in ZirFlu_strains) {
  reclassify_temp[[strain]] <- get.reclassify(dat = ZirFlu$HAItiter, 
                                              baseline_col = paste0(strain, "_T1"),
                                              abFC_col = paste0(strain, "_abFC"),
                                              baseline_cutoff = baseline,
                                              abFC_cutoff = abFC)
  
}

ZirFlu_HAIreclassify <-  ZirFlu$HAItiter %>% 
  cbind(as.data.frame(reclassify_temp) %>% 
          rename_with(~paste0(.x, "_reclassify")))
rm(reclassify_temp)

# check groups
ZirFlu_HAIreclassify %>% count(season, H1N1_reclassify)
ZirFlu_HAIreclassify %>% count(season, H3N2_reclassify)
ZirFlu_HAIreclassify %>% count(season, Bvictoria_reclassify)
ZirFlu_HAIreclassify %>% count(season, Byamagata_reclassify)

## box plot 
seasons <- c("2019", "2020")
plot_list <- list()
for (year in seasons) {
  HAI_reclassify_temp <- ZirFlu_HAIreclassify  %>% filter(season == year) #%>% select(-matches("T3"))
  plot_list[[year]] <- get.boxplot_reclassify_perStrain(dat = HAI_reclassify_temp, 
                                                        strain_list = ZirFlu_strains)
}

cowplot::plot_grid(plot(plot_list$`2019`$H1N1),
                   plot(plot_list$`2019`$H3N2),
                   plot(plot_list$`2019`$Bvictoria),
                   plot(plot_list$`2019`$Byamagata),
                   nrow = 4)

cowplot::plot_grid(plot(plot_list$`2020`$H1N1),
                   plot(plot_list$`2020`$H3N2),
                   plot(plot_list$`2020`$Bvictoria),
                   plot(plot_list$`2020`$Byamagata),
                   nrow = 4)

# iMED =====================================================================
iMED_strains <- c("H1N1", "H3N2", "B")
reclassify_temp <- list()
for (strain in iMED_strains) {
  reclassify_temp[[strain]] <- get.reclassify(dat = iMED$HAItiter, 
                                              baseline_col = paste0(strain, "_T1"),
                                              abFC_col = paste0("ab_", strain),
                                              baseline_cutoff = baseline,
                                              abFC_cutoff = abFC)
  
}

iMED_HAIreclassify <-  iMED$HAItiter %>% 
  cbind(as.data.frame(reclassify_temp) %>% 
          rename_with(~paste0(.x, "_reclassify")))
rm(reclassify_temp)

# check groups
iMED_HAIreclassify %>% count(H1N1_reclassify)
iMED_HAIreclassify %>% count(H3N2_reclassify)
iMED_HAIreclassify %>% count(B_reclassify)

# boxplot
iMED_reclass_plot <- get.boxplot_reclassify_perStrain(dat = iMED_HAIreclassify,
                                                      strain_list = iMED_strains)

cowplot::plot_grid(plot(iMED_reclass_plot$H1N1),
                   plot(iMED_reclass_plot$H3N2),
                   plot(iMED_reclass_plot$B),
                   nrow = 3)
# save data ---------------------------------------
HAIreclassify2 <- list()

HAIreclassify2$ZirFlu <- ZirFlu_HAIreclassify
HAIreclassify2$iMED <- iMED_HAIreclassify

HAIreclassify2$all_cohorts <- HAIreclassify2 %>% purrr::reduce(full_join) %>%
  dplyr::select(patientID, season, matches("reclassify")) %>%
  mutate(Bvictoria_reclassify2 = ifelse(is.na(Bvictoria_reclassify), B_reclassify, Bvictoria_reclassify)) %>%
  mutate(Byamagata_reclassify2 = ifelse(is.na(Byamagata_reclassify), B_reclassify, Byamagata_reclassify)) %>% 
  mutate(sameResponse_2Bstrains = ifelse(Bvictoria_reclassify == Byamagata_reclassify, TRUE, FALSE)) 

HAIreclassify2$all_cohorts %>% count(sameResponse_2Bstrains)

save(HAIreclassify2, file = "HAIreclassify2.RData")


# check iMED --------------
iMED$HAItiter %>% count(H1N1_T1 >= log2(40))
iMED$HAItiter %>% count(ab_H1N1 >= log2(4))

