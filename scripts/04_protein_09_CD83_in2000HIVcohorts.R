# Aim: using he 2000HIV cohort (patient with HIV) tp check the relationship of 
# CD83 with other CD83-related proteins at protein and RNAseq levels
# CD83 with cytokine level


rm(list = ls())

library(tidyverse)
library(rstatix)
library(ggpubr)

plot.cor <- function(dat, var1, var2) {
  dat %>%
    ggplot(aes_string(x = var1, y = var2))  +
    geom_point() +
    geom_smooth(method = "lm", se=FALSE) + stat_cor() + 
    theme_classic()
}

plot.corSpearman <- function(dat, var1, var2) {
  dat %>%
    ggplot(aes_string(x = var1, y = var2))  +
    geom_point() +
    geom_smooth(method = "lm", se=FALSE) + stat_cor(method = "spearman") + 
    theme_classic()
}
# load data =======================================================================
load("/vol/projects/BIIM/cohorts_CMV/CMV_2000HIV_NhanNguyen/processedDat/cohortDat.RData") # load the data form 2000HIV cohort

CD83_relatedPro <- c("CD83", "CD1A", "CD1B", "CD1C", "CD1D", "CD1E", "CD40", "CD80", "CD86", "TNF", "CCR7") # based on stringDB

# cor(CD83, other proteins) at protein level ---------------------------
matched_proteins <- intersect(CD83_relatedPro, names(cohortDat$allSample$protein)) # 6 proteins
HIVcohorts_protein <- cohortDat$allSample$protein %>% select(matched_proteins)

outcome <- matrix(NA, nrow = length(matched_proteins), ncol = 2) %>% as.data.frame
names(outcome) <- c("p.value", "cor.estimate")
rownames(outcome) <- matched_proteins

for (i in 1:length(matched_proteins)) {
  outcome$p.value[i] <- cor.test(HIVcohorts_protein$CD83, HIVcohorts_protein [[matched_proteins[i]]])$p.value
  outcome$cor.estimate[i] <- cor.test(HIVcohorts_protein$CD83, HIVcohorts_protein [[matched_proteins[i]]])$estimate
}
outcome$padj <- p.adjust(outcome$p.value, method = "fdr")
outcome %>% filter(p.value < 0.05)
outcome %>% filter(padj < 0.05)

## scatter plots ----------------------------------
plot_CD83_CD1C <- plot.cor(HIVcohorts_protein, "CD83", "CD1C")
plot_CD83_CD40 <- plot.cor(HIVcohorts_protein, "CD83", "CD40")
plot_CD83_CD80 <- plot.cor(HIVcohorts_protein, "CD83", "CD80")
plot_CD83_CD86 <- plot.cor(HIVcohorts_protein, "CD83", "CD86")
plot_CD83_TNF <- plot.cor(HIVcohorts_protein, "CD83", "TNF")

cowplot::plot_grid(plot_CD83_CD1C, plot_CD83_CD40, plot_CD83_CD80, plot_CD83_TNF, nrow =1)


# cor(CD83, other genes) at transcriptome level -------------------------
matched_rnas <- intersect(CD83_relatedPro, names(cohortDat$allSample$rna)) # 10 genes
HIVcohorts_rna <- cohortDat$allSample$rna %>% select(matched_rnas)

outcome <- matrix(NA, nrow = length(matched_rnas), ncol = 2) %>% as.data.frame
names(outcome) <- c("p.value", "cor.estimate")
rownames(outcome) <- matched_rnas

for (i in 1:length(matched_rnas)) {
  outcome$p.value[i] <- cor.test(HIVcohorts_rna$CD83, HIVcohorts_rna[[matched_rnas[i]]])$p.value
  outcome$cor.estimate[i] <- cor.test(HIVcohorts_rna$CD83, HIVcohorts_rna[[matched_rnas[i]]])$estimate
}
outcome$padj <- p.adjust(outcome$p.value, method = "fdr")
outcome %>% filter(p.value < 0.05)
outcome %>% filter(padj < 0.05)

## scatter plots ----------------------------------
plot_CD83_CD1C <- plot.cor(HIVcohorts_rna, "CD83", "CD1C")
plot_CD83_CD40 <- plot.cor(HIVcohorts_rna, "CD83", "CD40")
plot_CD83_CD80 <- plot.cor(HIVcohorts_rna, "CD83", "CD80")
plot_CD83_CD86 <- plot.cor(HIVcohorts_rna, "CD83", "CD86")
plot_CD83_TNF <- plot.cor(HIVcohorts_rna, "CD83", "TNF")
plot_CD83_CCR7 <- plot.cor(HIVcohorts_rna, "CD83", "CCR7")

cowplot::plot_grid(plot_CD83_CD1C, plot_CD83_CD40, plot_CD83_CD80, 
                   plot_CD83_TNF, plot_CD83_CCR7, nrow =1)

# cor(CD83, cytokine) -------------------------
HIVcohorts_cytokineCD83 <- cohortDat$allSample$cytokine %>% 
  rownames_to_column("patient") %>%
  full_join(cohortDat$allSample$protein %>% select("CD83") %>% 
               rownames_to_column("patient")) %>%
  column_to_rownames("patient")

cytokines <- names(cohortDat$allSample$cytokine)

outcome <- matrix(NA, nrow = length(cytokines), ncol = 2) %>% as.data.frame
names(outcome) <- c("p.value", "cor.estimate")
rownames(outcome) <- cytokines

for (i in 1:length(cytokines)) {
  outcome$p.value[i] <- cor.test(HIVcohorts_cytokineCD83$CD83, 
                                 HIVcohorts_cytokineCD83[[cytokines[i]]], method = "spearman")$p.value
  outcome$cor.estimate[i] <- cor.test(HIVcohorts_cytokineCD83$CD83, 
                                      HIVcohorts_cytokineCD83[[cytokines[i]]], method = "spearman")$estimate
}
outcome$padj <- p.adjust(outcome$p.value, method = "fdr")
outcome %>% filter(p.value < 0.05)
outcome %>% filter(padj < 0.05)

# scatter plots ----------------------------------
plotList <- list()
for(sigCytos in rownames(outcome %>% filter(padj < 0.05))) {
  plotList[[sigCytos]] <- plot.corSpearman(HIVcohorts_cytokineCD83, "CD83", sigCytos)
}

cowplot::plot_grid(plotlist = plotList,
                   nrow = 3)
