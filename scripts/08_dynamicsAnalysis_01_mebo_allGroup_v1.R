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

metabols <- mebo_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% 
  column_to_rownames("name") %>% t() %>% as.data.frame %>% 
  select(metadat$name)

all.equal(colnames(metabols), metadat$name) # TRUE = the order of samples are the same

# linear fit() model with time ------------------------------------
years <- c("2014", "2015", "2019", "2020")
outcome <- list()
for (year in years) {
  # prepare inputs ------------------------------------
  metadat_temp <- metadat %>% filter(season == year)
  metabols_temp <- metabols %>% select(metadat_temp$name)
  
  # prepare design table ------------------------------------
  modTime <- model.matrix(~0 + time + age + sex, data = metadat_temp)
  
  corfit <- duplicateCorrelation(metabols_temp, modTime, block = metadat_temp$probandID)
  print(paste0(year, " - corfit$consensus.correlation ", corfit$consensus.correlation))
  
  # per timepoint ------------------------------------
  fit_perTime <- lmFit(metabols_temp, modTime, 
                       block = metadat_temp$probandID, 
                       correlation = corfit$consensus.correlation) %>% eBayes()
  
  # overtime ------------------------------------
  cont.matrix = makeContrasts(d7vsd0 = timed7 - timed0,
                              d28vsd0 = timed28 - timed0, levels = modTime)
  
  fit_overTime = contrasts.fit(fit_perTime, cont.matrix) %>% eBayes()
  
  # outcome ------------------------------------
  outcome[[year]]$fit_perTime <- summary(decideTests(fit_perTime))
  outcome[[year]]$fit_overTime <- summary(decideTests(fit_overTime))
  outcome[[year]]$compare_overTime <- list(
    d7vsd0 = topTable(fit_overTime, coef = 1, number = Inf),
    d28vsd0 = topTable(fit_overTime, coef = 2, number = Inf))
  
}


# plot the number of  DE overtime -----------------------------------
baseline <- data.frame(time = rep("d0", 8), value = rep(0, 8), 
                       variable = c(rep("Down", 4), rep("Up", 4)))

plotDat <- outcome %>% 
  lapply(function(x) x$fit_overTime %>% as.data.frame() %>%
           rename("variable" = "Var1", "time" = "Var2", "value" = "Freq") %>% 
           filter(variable != "NotSig") %>%
           mutate(time = gsub("vsd0", "", time),
                  value = ifelse(variable == "Down", -value, value)) %>%
           rbind(baseline) %>%
           mutate(time = factor(time, levels = c("d0", "d7", "d28"))))

plotDat_v2 <- plotDat %>% bind_rows(.id = "season")

## plot the number of  DE overtime per group -----------------------------------
plotDat$`2014` %>% 
  ggplot(aes(x=time, y=value))+
  geom_line(position="identity", size = 1)+
  geom_point(size=3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites")

plotDat_v2 %>% 
  ggplot(aes(x=time, y=value, group = season))+
  geom_line(position="identity",aes(col=season),size = 1)+
  geom_point(aes(col=season),size=3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites")

# get the DE overtime -----------------------------------
sigOutcome_temp <- outcome %>% 
  lapply(function(x) x$compare_overTime %>% 
           lapply(function(y) y %>% filter(adj.P.Val < 0.05)))

sigOutcome <- outcome  %>% 
  lapply(function(x) x$compare_overTime %>% 
           lapply(function(y) y %>% rownames_to_column("var")) %>% 
           bind_rows(.id = "compare")) %>% 
  bind_rows(.id = "season") %>% filter(adj.P.Val < 0.05)

length(unique(sigOutcome$var)) # check pathways they belong to, and which pathway is important. Show the change per pathway

# not rund this code yet: save(outcome, sigOutcome, file = "sigMebo_dynamic.RData")

b1 <- sigOutcome %>% filter(compare %in% c("d7vsd0")) %>%
  group_by(var) %>% add_count()
unique((b1 %>% filter(n == 2))$var)

b2 <- sigOutcome %>% filter(compare %in% c("d28vsd0")) %>%
  group_by(var) %>% add_count()
unique((b2 %>% filter(n == 2))$var)
