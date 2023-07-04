rm(list = ls())

library(tidyverse)
library(ggvenn)
# sig. protein dynamic for 4 groups: LL, LH, HL, HH ===========================
load("/vol/projects/CIIM/Influenza/iMED/proteomic/protein_reclassify_timeDynamics.RData")

allPro_4groups <- list(
  "LL_T2andT3vsT1" = c(rownames(LL_T2vsT1_prot %>% filter(adj.P.Val < 0.05)),
                       rownames(LL_T3vsT1_prot %>% filter(adj.P.Val < 0.05))),
  "LH_T2andT3vsT1" = c(rownames(LH_T2vsT1_prot %>% filter(adj.P.Val < 0.05)),
                       rownames(LH_T3vsT1_prot %>% filter(adj.P.Val < 0.05))),
  "HL_T2andT3vsT1" = c(rownames(HL_T2vsT1_prot %>% filter(adj.P.Val < 0.05)),
                       rownames(HL_T3vsT1_prot %>% filter(adj.P.Val < 0.05))),
  "HH_T2andT3vsT1" = c(rownames(HH_T2vsT1_prot %>% filter(adj.P.Val < 0.05)),
                       rownames(HH_T3vsT1_prot %>% filter(adj.P.Val < 0.05))))

## venn diagram of sig. proteins for 4 groups ----------------
ggvenn(allPro_4groups,
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 4)

# load data ==============================================================
rm(list = ls())
load("proDynamic_season2015.RData")
load("selected_DAPs.RData")

proSig_tab <- proteins_comparison%>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = group, values_from = value)

intersect(selected_DAs, proSig_tab$valName) # no overlapped proteins

# interested sig. proteins --------------------------------------------------
load("/vol/projects/CIIM/Influenza/iMED/proteomic/protein_reclassify_timeDynamics.RData")

### check overlap between sig. mebo at baseline and during dynamics --------------
interested_pro <- list("baseline_sigMebo" = selected_DAs,
                        "dynamic_sigMebo" = proSig_tab$valName)
ggvenn(interested_pro,
       # fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 4)

pros <- unique(unlist(interested_pro))

# heatmap with t-statistic =====================================================
allProChangingTstat <- list(
  "LL_T2vsT1" = LL_T2vsT1_prot %>% select(t),
  "LL_T3vsT1" = LL_T3vsT1_prot %>% select(t),
  "LH_T2vsT1" = LH_T2vsT1_prot %>% select(t),
  "LH_T3vsT1" = LH_T3vsT1_prot %>% select(t),
  "HL_T2vsT1" = HL_T2vsT1_prot %>% select(t),
  "HL_T3vsT1" = HL_T3vsT1_prot %>% select(t),
  "HH_T2vsT1" = HH_T2vsT1_prot %>% select(t),
  "HH_T3vsT1" = HH_T3vsT1_prot %>% select(t)) %>% 
  lapply(function(x) x %>% rownames_to_column("valName")) %>%
  bind_rows(.id = "group")

plotDat <- allProChangingTstat %>% 
  filter(valName %in% pros) %>%
  mutate(group = factor(group, 
                        levels = c("LL_T2vsT1", "LL_T3vsT1", 
                                   "LH_T2vsT1", "LH_T3vsT1", 
                                   "HL_T2vsT1", "HL_T3vsT1", 
                                   "HH_T2vsT1", "HH_T3vsT1"))) %>% left_join(proteins_comparison)

plotDat %>%
  mutate(group = factor(group, 
                        levels = c("LL_T2vsT1", "LL_T3vsT1", 
                                   "LH_T2vsT1", "LH_T3vsT1", 
                                   "HL_T2vsT1", "HL_T3vsT1", 
                                   "HH_T2vsT1", "HH_T3vsT1"))) %>%
  ggplot(aes(x = group, y = valName, fill = t)) + 
  geom_tile() +
  scale_fill_gradientn(limits = c(-8, 8), colors = c("blue", "white", "red"),  na.value = "grey") +
  geom_text(aes(label = ifelse(is.na(sigPvalue), NA, "*"))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


# Heatmap with protein levels ==================================
## prepare the data --------------------------------------------------
load("cohorts_dat.RData")

## metadata for all healthy subjects 
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

pros_v2 <- intersect(colnames(inputDat), pros) # use protein exist in all seasons
plotDat <- inputDat %>% 
  select(probandID, season, responder, time, c(pros_v2),
         matches("_abFC|_T1|_T4|_reclassify"))

plotDat_H1N1 <- plotDat %>% 
  pivot_longer(cols = pros_v2, names_to = "valName", values_to =  "intensity") %>%
  mutate(group = paste0(season, "_", H1N1_reclassify, "_", time)) %>%
  mutate(group = factor(group, 
                        levels = c("2014_LL_T1", "2014_LL_T3", "2014_LL_T4",
                                   "2014_LH_T1", "2014_LH_T3", "2014_LH_T4",
                                   "2014_HL_T1", "2014_HL_T3", "2014_HL_T4",
                                   "2014_HH_T1", "2014_HH_T3", "2014_HH_T4",
                                   "2015_LL_T1", "2015_LL_T3", "2015_LL_T4",
                                   "2015_LH_T1", "2015_LH_T3", "2015_LH_T4",
                                   "2015_HL_T1", "2015_HL_T3", "2015_HL_T4",
                                   "2015_HH_T1", "2015_HH_T3", "2015_HH_T4",
                                   "2019_LL_T1", "2019_LL_T4", "2019_LL_T5",
                                   "2019_LH_T1", "2019_LH_T4", "2019_LH_T5",
                                   "2019_HL_T1", "2019_HL_T4", "2019_HL_T5",
                                   "2019_HH_T1", "2019_HH_T4", "2019_HH_T5",
                                   "2020_LL_T1", "2020_LL_T4", "2020_LL_T5",
                                   "2020_LH_T1", "2020_LH_T4", "2020_LH_T5",
                                   "2020_HL_T1", "2020_HL_T4", "2020_HL_T5")))
## heatmap with protein level ------------------------------
plotDat_H1N1 %>% slice(grep("2015", group)) %>%
  ggplot(aes(x = group, y = valName, fill = intensity)) + 
  geom_tile() +
  scale_fill_gradientn(limits = c(-9, 9), colors = c("blue", "white", "red"),  na.value = "grey") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  ggtitle("Sig. protein in dynamic analysis, show in H1N1 strain")

## heatmap with the average protein level ------------------------------
plotDat_H1N1_avg <- plotDat_H1N1 %>% 
  select(valName, intensity, group) %>% group_by(group, valName) %>% summarise(avg = mean(intensity, na.rm = TRUE))

plotDat_H1N1_avg %>% slice(grep("2015", group)) %>%
  ggplot(aes(x = group, y = valName, fill = avg)) + 
  geom_tile() +
  scale_fill_gradientn(limits = c(-6, 6), colors = c("blue", "white", "red"),  na.value = "grey") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  ggtitle("Avg sig. metabolite in dynamic analysis, show in H1N1 strain")

### heatmap with scale row (use complex heatmap) ------------------------------
library(ComplexHeatmap)

plotDat_H1N1_wide <- plotDat %>% 
  select(season, time, H1N1_reclassify, pros_v2) %>% filter(season == "2015") %>%
  mutate(group = paste0(season, "_", H1N1_reclassify, "_", time)) %>%
  group_by(group) %>% 
  summarise(across(pros_v2, function(x) mean(x, na.rm = TRUE))) %>%
  column_to_rownames("group")

# scale the value
my.scale = function(x,na.rm=TRUE) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
plotDat_H1N1_wide_scale  <- plotDat_H1N1_wide  %>%
  mutate_all(scale) %>% as.data.frame() %>% t()

plotDat_H1N1_wide %>% t() %>% Heatmap()
plotDat_H1N1_wide_scale %>% 
  Heatmap(cluster_columns = FALSE,
          column_order = c( "2015_LL_T1", "2015_LL_T3", "2015_LL_T4",
                            "2015_LH_T1", "2015_LH_T3", "2015_LH_T4",
                            "2015_HL_T1", "2015_HL_T3", "2015_HL_T4",
                            "2015_HH_T1", "2015_HH_T3", "2015_HH_T4"),
          column_names_rot = 30)

group_annot <- interested_pro %>% 
  lapply(function(x) x %>% as.data.frame %>% rename("valName" = ".")) %>%
  bind_rows(.id = "group") %>% 
  mutate(value = TRUE) %>% pivot_wider(names_from = "group", values_from = "value") %>%
  mutate(group = ifelse(is.na(baseline_sigMebo), "sigDynamic", "sigBase")) %>%
  filter(valName %in% rownames(plotDat_H1N1_wide_scale)) %>%
  arrange(match(valName, rownames(plotDat_H1N1_wide_scale)))

identical(group_annot$valName, rownames(plotDat_H1N1_wide_scale)) # TRUE = the same order

row_ha = rowAnnotation(
  markers = group_annot$group, 
  col = list(markers = c("sigBase_sigDynamic" = "red", "sigBase" = "green", "sigDynamic" = "blue")))

plotDat_H1N1_wide_scale %>% 
  Heatmap(cluster_columns = FALSE, cluster_rows = FALSE,
          column_order = c( "2015_LL_T1", "2015_LL_T3", "2015_LL_T4",
                            "2015_LH_T1", "2015_LH_T3", "2015_LH_T4",
                            "2015_HL_T1", "2015_HL_T3", "2015_HL_T4",
                            "2015_HH_T1", "2015_HH_T3", "2015_HH_T4"),
          column_split = substring(colnames(plotDat_H1N1_wide_scale ), 6, 7),
          right_annotation = row_ha,
          column_names_rot = 30)

# split group for controling
group_temp <- "sigBase"
group_temp <- "sigDynamic"

plotDat_temp <- plotDat_H1N1_wide_scale[(group_annot %>% filter(group == group_temp))$valName, ]
row_ha_temp = rowAnnotation(
  markers = (group_annot %>% filter(group == group_temp ))$group, 
  col = list(markers = c("sigBase_sigDynamic" = "red", "sigBase" = "green", "sigDynamic" = "blue")))

plotDat_temp %>% 
  Heatmap(cluster_columns = FALSE, #cluster_rows = FALSE,
          column_order = c( "2015_LL_T1", "2015_LL_T3", "2015_LL_T4",
                            "2015_LH_T1", "2015_LH_T3", "2015_LH_T4",
                            "2015_HL_T1", "2015_HL_T3", "2015_HL_T4",
                            "2015_HH_T1", "2015_HH_T3", "2015_HH_T4"),
          column_split = substring(colnames(plotDat_H1N1_wide_scale ), 6, 7),
          right_annotation = row_ha_temp,
          column_names_rot = 30,
          row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
