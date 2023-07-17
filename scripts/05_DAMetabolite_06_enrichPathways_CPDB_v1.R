rm(list = ls())
library(tidyverse)

# load the data ----------------
ora_pathways <- read.table("processedDat/ORA_results.tab", sep = "\t", header = TRUE)

# plot the sig. pathway (padj < 0.05) -----------------
ora_pathways %>% filter(q.value < 0.05) %>%
  ggplot(aes(x = -log10(q.value), y = pathway, size = `effective_size`)) + 
  geom_point() + xlim(c(0, 3.5)) + theme_classic()
