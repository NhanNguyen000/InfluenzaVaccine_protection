rm(list = ls())
library(tidyverse)
library(openxlsx)

# load the data ----------------
load("selected_DAMs.RData") # 45 formulas
#load("selected_DAMs_padj2015.RData") # 146 formulas

# metabolite taxonomy with HMDB ids
mebo_taxo <- read.delim("/vol/projects/CIIM/Influenza/iMED/metabolic/db/hmdb/metabolite_HMDB_taxonomy.csv",
                        quote = "", header = TRUE)

# link the HMDB ids to Formulas
mebo_taxo_fomula <- mebo_taxo %>% rownames_to_column("CompoundID") %>% 
  inner_join(read.xlsx('/vol/projects/CIIM/Influenza/iMED/metabolic/raw_data/tables/DATA_CURATED_reformatted.xlsx',
                       sheet = 'annotation') %>% fill(ionIdx, .direction = "down")) %>% 
  filter(Formula %in% selected_DAs) %>%
  select(Formula, super_class, class, sub_class) %>% distinct()

length(unique(mebo_taxo_fomula$Formula)) # get all 45 formulas out of 45 formulas

mebo_taxo_fomula %>% count(super_class)
mebo_taxo_fomula %>% count(class)

# modify the class
mebo_taxo_fomula2 <- mebo_taxo_fomula %>% add_count(class) %>%
  mutate(metabolite_class = ifelse(n >= 5, class, "other"))

# pie chart for class proportion in DAMs---------------------
plotDat <- mebo_taxo_fomula2 %>% select(metabolite_class) %>% 
  add_count(metabolite_class) %>% arrange(metabolite_class) %>% distinct() %>% 
  mutate(prop = round(n / sum(n) *100, digits = 2))

plotDat %>% 
  ggplot(aes(x = "", y = n, fill = metabolite_class)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  geom_text(aes(x = 1.3, label = prop), position = position_stack(vjust = 0.5),
            color = "black", size = 3) +
  coord_polar("y") + theme_void()

  
library(webr)
PieDonut(mebo_taxo_fomula2, aes(super_class, class), 
         r0 = 0, r2 = 1.1, start = -120,
         explode = 3, explodeDonut = TRUE,
         showPieName = FALSE, 
         title = "Distribution of gender per season")

