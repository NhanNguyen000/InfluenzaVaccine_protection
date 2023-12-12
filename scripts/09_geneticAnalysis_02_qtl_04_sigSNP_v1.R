rm(list = ls())
library(tidyverse)

# load data =================================================
times <- c("T1", "T4")

sig_SNIPs <- list()
for (time in times) {
  sig_SNIPs[[time]] <- read.table(
    file = paste0("processedDat/qtl/outcome_", time, ".csv"), header = TRUE) %>% 
    as.data.frame %>% filter(p.value < 10e-5)
}

## info for sig. SNPs --------------------------------------
loci <- sig_SNIPs %>%  bind_rows(.id = "time")

loci_v2 <- loci %>% group_by(SNP) %>% add_count() %>%
  filter(n > 1)
