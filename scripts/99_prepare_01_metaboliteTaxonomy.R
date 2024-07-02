rm(list = ls())

library(xml2)
library(tidyverse)

# load data --------------------------------------------------------
hmdb_xml <- as_list(read_xml("data/hmdb_metabolites.xml"))
names(hmdb_xml$hmdb[[1]])

save(hmdb_xml, file = "processedDat/hmdb_xml.RData")

xml_structure(hmdb_xml)
xml_structure(hmdb_xml$hmdb)

xml_df = tibble::as_tibble(hmdb_xml$hmdb) %>%
  unnest_longer()


a <- hmdb_xml$hmdb %>% lapply(function(x) x[c("chemical_formula", "taxonomy")])
k <- a %>% lapply(function(x) x$chemical_formula) %>% unlist()
names(a) <- k

a2 <- a %>% lapply(function(x) x$taxonomy)
names(a2$C7H11N3O2)
b <- a2 %>% 
  lapply(function(x) x[c("description", "kingdom", "super_class", "class", "sub_class", "molecular_framework")])

b2 <- b %>% lapply(function(x) x %>% unlist() %>% as.data.frame)

xml_df <- as_tibble(b) %>%
  unnest_longer()

b$C7H11N3O2$description %>%
  unnest_longer()
g <- b$C7H11N3O2 %>% unlist() %>% as.data.frame()







a <- hmdb_xml$hmdb %>% lapply(function(x) x[c("chemical_formula", "taxonomy")])
k <- a %>% lapply(function(x) x$chemical_formula) %>% unlist()
names(a) <- k

a2 <- a %>% lapply(function(x) x$taxonomy)
names(a2$C7H11N3O2)
b <- a2 %>% 
  lapply(function(x) x[c("description", "kingdom", "super_class", "class", "sub_class", "molecular_framework")])

b2 <- b %>% lapply(function(x) x %>% unlist() %>% as.data.frame)



q2 <- hmdb_xml$hmdb %>% 
  lapply(function(x) 
    x$taxonomy[c("description", "kingdom", "super_class", "class", "sub_class", "molecular_framework")])
names(q2) <-  hmdb_xml$hmdb %>% lapply(function(x) x$chemical_formula) %>% unlist()
q3 <- q2 %>% 
  lapply(function(x) x %>% unlist() %>% as.data.frame %>% rownames_to_column("name"))
q4 <- q3 %>% imap(~.x %>% rename_with(function(x) paste(x), 1))

w3 <- q2 %>% 
  lapply(function(x) x %>% unlist() %>% as.data.frame) %>% 
  imap(~.x %>% rename_with(function(x) paste(.x), 1))


