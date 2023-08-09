rm(list = ls())
library(tidyverse)

# load the data ----------------
ora_pathways <- read.table("processedDat/ORA_results_DAM_4reclass_to2groups.tab", sep = "\t", header = TRUE)
ora_pathways <- read.table("processedDat/ORA_results_DAM_4groups.tab", sep = "\t", header = TRUE)

# plot the sig. pathway (padj < 0.05) -----------------
ora_pathways %>% filter(q.value < 0.05) %>%
  ggplot(aes(x = -log10(q.value), y = pathway, size = `effective_size`)) + 
  geom_point() + xlim(c(0, 3.5)) + theme_classic()

# check the pathway ---------------------
fatty_acids <- ora_pathways %>% 
  filter(q.value < 0.05) %>% slice(grep("acid", pathway)) %>%
  select(members_input_overlap) %>% unlist() %>% strsplit("; ") %>% 
  unlist() %>% unique()
