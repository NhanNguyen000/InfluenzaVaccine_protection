# adapt from Saumya code
rm(list = ls())
library(tidyverse)
library(limma)
library(reshape2)

# load data =======================================================================
load("/vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/cohorts_dat.RData")

metadat <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all) %>%
  filter(condition == "Healthy")

proteinDat <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% 
  column_to_rownames("name") %>% t() %>% as.data.frame %>% 
  select(metadat$name)

all.equal(colnames(proteinDat), metadat$name) # TRUE = the order of samples are the same

# linear fit() model with time ------------------------------------
## strain with season have all 4 groups: LL, LH, HL, HH -----------------------------------
groups <- c("2014_H1N1", "2015_H1N1", "2019_H1N1", #"2020_H1N1",
            "2014_B", "2015_B", "2015_H3N2", "2020_Byamagata")

outcome <- list()
for (group in groups) {
  # prepare inputs ------------------------------------
  metadat_temp <- metadat %>% filter(season == substring(group, 1, 4))
  proteinDat_temp <- proteinDat %>% select(metadat_temp$name)
  
  strain <- substring(group, 6, nchar(group))
  
  # prepare design table ------------------------------------
  metaFileMerged <- metadat_temp %>%
    mutate(reclassifyTime = paste0(metadat_temp[, c(paste0(strain, "_reclassify"))], "_", time))
  
  modTime <- model.matrix(~0 + reclassifyTime + age + sex, data = metaFileMerged)
  
  corfit <- duplicateCorrelation(proteinDat_temp, modTime, block = metaFileMerged$probandID)
  print(paste0(group, " - corfit$consensus.correlation ", corfit$consensus.correlation))
  
  # per timepoint ------------------------------------
  fit_perTime <- lmFit(proteinDat_temp, modTime, 
                       block = metaFileMerged$probandID, 
                       correlation = corfit$consensus.correlation) %>% eBayes()
  
  # overtime ------------------------------------
  cont.matrix = makeContrasts(HH_d7vsd0 = reclassifyTimeHH_d7 - reclassifyTimeHH_d0,
                              HH_d28vsd0 = reclassifyTimeHH_d28 - reclassifyTimeHH_d0,
                              LH_d7vsd0 = reclassifyTimeLH_d7 - reclassifyTimeLH_d0,
                              LH_d28vsd0 = reclassifyTimeLH_d28 - reclassifyTimeLH_d0,
                              HL_d7vsd0 = reclassifyTimeHL_d7 - reclassifyTimeHL_d0,
                              HL_d28vsd0 = reclassifyTimeHL_d28 - reclassifyTimeHL_d0,
                              LL_d7vsd0 = reclassifyTimeLL_d7 - reclassifyTimeLL_d0,
                              LL_d28vsd0 = reclassifyTimeLL_d28 - reclassifyTimeLL_d0, levels = modTime)
  
  fit_overTime = contrasts.fit(fit_perTime, cont.matrix) %>% eBayes()
  
  # outcome ------------------------------------
  outcome[[group]]$fit_perTime <- summary(decideTests(fit_perTime))
  outcome[[group]]$fit_overTime <- summary(decideTests(fit_overTime))
  outcome[[group]]$compare_overTime <- list(
    HH_d7vsd0 = topTable(fit_overTime, coef = 1, number = Inf),
    HH_d28vsd0 = topTable(fit_overTime, coef = 2, number = Inf),
    LH_d7vsd0 = topTable(fit_overTime, coef = 3, number = Inf),
    LH_d28vsd0 = topTable(fit_overTime, coef = 4, number = Inf),
    HL_d7vsd0 = topTable(fit_overTime, coef = 5, number = Inf),
    HL_d28vsd0 = topTable(fit_overTime, coef = 6, number = Inf),
    LL_d7vsd0 = topTable(fit_overTime, coef = 7, number = Inf),
    LL_d28vsd0 = topTable(fit_overTime, coef = 8, number = Inf))
}

## strain with season have all 3 groups: LL, LH, HL-----------------------------------
groups <- c("2020_H1N1")

for (group in groups) {
  # prepare inputs ------------------------------------
  metadat_temp <- metadat %>% filter(season == substring(group, 1, 4))
  proteinDat_temp <- proteinDat %>% select(metadat_temp$name)
  
  strain <- substring(group, 6, nchar(group))
  
  # prepare design table ------------------------------------
  metaFileMerged <- metadat_temp %>%
    mutate(reclassifyTime = paste0(metadat_temp[, c(paste0(strain, "_reclassify"))], "_", time))
  
  modTime <- model.matrix(~0 + reclassifyTime + age + sex, data = metaFileMerged)
  
  corfit <- duplicateCorrelation(proteinDat_temp, modTime, block = metaFileMerged$probandID)
  print(paste0(group, " - corfit$consensus.correlation ", corfit$consensus.correlation))
  
  # per timepoint ------------------------------------
  fit_perTime <- lmFit(proteinDat_temp, modTime, 
                       block = metaFileMerged$probandID, 
                       correlation = corfit$consensus.correlation) %>% eBayes()
  
  # overtime ------------------------------------
  cont.matrix = makeContrasts(LH_d7vsd0 = reclassifyTimeLH_d7 - reclassifyTimeLH_d0,
                              LH_d28vsd0 = reclassifyTimeLH_d28 - reclassifyTimeLH_d0,
                              HL_d7vsd0 = reclassifyTimeHL_d7 - reclassifyTimeHL_d0,
                              HL_d28vsd0 = reclassifyTimeHL_d28 - reclassifyTimeHL_d0,
                              LL_d7vsd0 = reclassifyTimeLL_d7 - reclassifyTimeLL_d0,
                              LL_d28vsd0 = reclassifyTimeLL_d28 - reclassifyTimeLL_d0, levels = modTime)
  
  fit_overTime = contrasts.fit(fit_perTime, cont.matrix) %>% eBayes()
  
  # outcome ------------------------------------
  outcome[[group]]$fit_perTime <- summary(decideTests(fit_perTime))
  outcome[[group]]$fit_overTime <- summary(decideTests(fit_overTime))
  outcome[[group]]$compare_overTime <- list(
    LH_d7vsd0 = topTable(fit_overTime, coef = 1, number = Inf),
    LH_d28vsd0 = topTable(fit_overTime, coef = 2, number = Inf),
    HL_d7vsd0 = topTable(fit_overTime, coef = 3, number = Inf),
    HL_d28vsd0 = topTable(fit_overTime, coef = 4, number = Inf),
    LL_d7vsd0 = topTable(fit_overTime, coef = 5, number = Inf),
    LL_d28vsd0 = topTable(fit_overTime, coef = 6, number = Inf))
}



# plot the number of  DE overtime -----------------------------------
baseline <- data.frame(time = rep("d0", 8), value = rep(0, 8), 
                       variable = c(rep("Down", 4), rep("Up", 4)),
                       Category = rep(c("LL", "LH", "HL", "HH"), 2))

plotDat <- outcome %>% 
  lapply(function(x) x$fit_overTime %>% as.data.frame() %>%
           rename("variable" = "Var1", "value" = "Freq") %>% 
           filter(variable != "NotSig") %>%
           separate(Var2, into = c("Category", "time")) %>%
           mutate(time = gsub("vsd0", "", time),
                  value = ifelse(variable == "Down", -value, value)) %>%
           rbind(baseline) %>%
           mutate(CategoryVariable = paste0(Category,"_",variable),
                  time = factor(time, levels = c("d0", "d7", "d28"))))

## plot the number of  DE overtime per group -----------------------------------
plotDat$`2015_H1N1` %>% 
  ggplot(aes(x=time, y=value, group=CategoryVariable))+
  geom_line(position="identity",aes(col=Category),size = 1)+
  geom_point(aes(col=Category),size=3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites")

plotDat$`2015_H3N2` %>% 
  ggplot(aes(x=time, y=value, group=CategoryVariable))+
  geom_line(position="identity",aes(col=Category),size = 1)+
  geom_point(aes(col=Category),size=3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites")

plotDat$`2015_B` %>% 
  ggplot(aes(x=time, y=value, group=CategoryVariable))+
  geom_line(position="identity",aes(col=Category),size = 1)+
  geom_point(aes(col=Category),size=3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites")

## plot the number of dynamic DE per group together -----------------
plotDat_v2 <- plotDat %>% bind_rows(.id = "group")

plotDat_v2 %>% 
  filter(group %in% c("2014_H1N1", "2015_H1N1", "2014_B", "2015_B", "2015_H3N2")) %>%
  ggplot(aes(x = time, y = value, group = CategoryVariable))+
  geom_line(position = "identity", aes(col = Category),size = 1)+
  geom_point(aes(col = Category), size = 3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites") +
  facet_wrap(~group)

plotDat_v2 %>% 
  filter(group %in% c("2019_H1N1", "2020_H1N1", "2020_Byamagata")) %>%
  ggplot(aes(x = time, y = value, group = CategoryVariable))+
  geom_line(position = "identity", aes(col = Category),size = 1)+
  geom_point(aes(col = Category), size = 3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites") +
  facet_wrap(~group, ncol = 1)

# get the DE overtime -----------------------------------
sigOutcome_temp <- outcome %>% 
  lapply(function(x) x$compare_overTime %>% 
           lapply(function(y) y %>% filter(adj.P.Val < 0.05)))

sigOutcome <- outcome  %>% 
  lapply(function(x) x$compare_overTime %>% 
           lapply(function(y) y %>% rownames_to_column("var")) %>% 
           bind_rows(.id = "compare")) %>% 
  bind_rows(.id = "group") %>% filter(adj.P.Val < 0.05)

length(unique(sigOutcome$var)) # check pathways they belong to, and which pathway is important. Show the change per pathway

save(outcome, sigOutcome, file = "sigMebo_dynamic.RData")