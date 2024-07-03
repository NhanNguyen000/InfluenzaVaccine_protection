rm(list = ls())

library(tidyverse)
library(venn)

# load human atlas: gene epression in cell type -----------
human_atlas <- read.table("data/rna_single_cell_type_tissue.tsv", 
                header = TRUE, sep = "\t", fill = NA)

# load DAPs ----------------------
load("selected_DAPs.RData")
selected_DAs

#selected_tissues <- c("pbmc", "lymph node", "bone marrow")
selected_tissues <- "pbmc"

dat <- human_atlas %>% 
  filter(Gene.name %in% selected_DAs,
         Tissue %in% selected_tissues)
dat2 <- dat %>% select(Gene.name, Tissue, Cell.type) %>% distinct()

dat3 <- dat %>% 
  select(Gene.name, Cell.type, nTPM) %>% distinct() %>%
  group_by(Gene.name, Cell.type) %>% summarise(mean_nTPM = mean(nTPM)) %>%
  mutate(log2_mean_nTPM = ifelse(mean_nTPM == 0, 0, log2(mean_nTPM)))

unique(dat3$Cell.type)

dat3 %>%
  ggplot(aes(x = Cell.type, y = Gene.name, fill = log2_mean_nTPM)) + 
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
