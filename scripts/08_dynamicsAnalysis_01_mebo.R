# adapt from Saumya code
rm(list = ls())
library(tidyverse)
library(limma)
library(reshape2)

# load data =======================================================================
load("/vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/cohorts_dat.RData")
metabols <- mebo_Dat$iMED_2015 %>% t()

metadat <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all) %>%
  filter(season == "2015", condition == "Healthy") %>%
  arrange(match(name, colnames(metabols))) 

all.equal(colnames(metabols), metadat$name) # TRUE = the order of samples are the same

# linear fit() model with time ------------------------------------
strains <- c("H1N1", "H3N2", "B")

outcome <- list()
for (strain in strains) {
  metaFileMerged <- metadat %>%
    mutate(reclassifyTime = paste0(metadat[, c(paste0(strain, "_reclassify"))], "_", time))
  
  modTime <- model.matrix(~0 + reclassifyTime + age + sex, data = metaFileMerged)
  
  corfit <- duplicateCorrelation(metabols,modTime,block = metaFileMerged$probandID)
  print(paste0(strain, " - corfit$consensus.correlation ", corfit$consensus.correlation))
  
  # per timepoint ------------------------------------
  fit_perTime <- lmFit(metabols, modTime, 
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
  outcome[[strain]]$fit_perTime <- summary(decideTests(fit_perTime))
  outcome[[strain]]$fit_overTime <- summary(decideTests(fit_overTime))
  outcome[[strain]]$compare_overTime <- list(
    HH_d7vsd0 = topTable(fit_overTime, coef = 1, number = Inf),
    HH_d28vsd0 = topTable(fit_overTime, coef = 2, number = Inf),
    LH_d7vsd0 = topTable(fit_overTime, coef = 3, number = Inf),
    LH_d28vsd0 = topTable(fit_overTime, coef = 4, number = Inf),
    HL_d7vsd0 = topTable(fit_overTime, coef = 5, number = Inf),
    HL_d28vsd0 = topTable(fit_overTime, coef = 6, number = Inf),
    LL_d7vsd0 = topTable(fit_overTime, coef = 7, number = Inf),
    LL_d28vsd0 = topTable(fit_overTime, coef = 8, number = Inf))
}

# plot the DE overtime -----------------------------------
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
           mutate(CategoryVariable = paste0(Category,"_",variable)))

plotDat$H1N1 %>% 
  ggplot(aes(x=time, y=value, group=CategoryVariable))+
  geom_line(position="identity",aes(col=Category),size = 1)+
  geom_point(aes(col=Category),size=3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites")

plotDat$H3N2 %>% 
  ggplot(aes(x=time, y=value, group=CategoryVariable))+
  geom_line(position="identity",aes(col=Category),size = 1)+
  geom_point(aes(col=Category),size=3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites")

plotDat$B %>% 
  ggplot(aes(x=time, y=value, group=CategoryVariable))+
  geom_line(position="identity",aes(col=Category),size = 1)+
  geom_point(aes(col=Category),size=3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites")
