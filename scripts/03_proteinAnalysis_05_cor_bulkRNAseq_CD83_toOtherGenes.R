rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data =======================================================================

## load data from iMED transcriptome ---------------
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
# View(transcriptomeList[[1]])
# View(transcriptomeList[[2]])
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1") # time T1 in trancriptome = d0
iMED_transcripDat <- transcriptomeList[[1]] %>% as.data.frame %>% 
  select(iMED_transcrip_T1$SampleName) 

bulkRNAseq_baseline <- iMED_transcripDat %>% 
  select(iMED_transcrip_T1$SampleName) %>% # time T1 in trancriptome = d0
  t() %>% as.data.frame

# cor(CD83, other proteins) at transcriptome level -----------------------------------
CD83_relatedPro <- c("CD1A", "CD1B", "CD1C", "CD1D", "CD1E", "CD40", "CD80", "CD86", "TNF", "CCR7") # based on stringDB

outcome <- matrix(NA, nrow = length(CD83_relatedPro), ncol = 2) %>% as.data.frame
names(outcome) <- c("p.value", "cor.estimate")
rownames(outcome) <- CD83_relatedPro

for (i in 1:length(CD83_relatedPro)) {
  outcome$p.value[i] <- cor.test(bulkRNAseq_baseline$CD83, bulkRNAseq_baseline[[CD83_relatedPro[i]]])$p.value
  outcome$cor.estimate[i] <- cor.test(bulkRNAseq_baseline$CD83, bulkRNAseq_baseline[[CD83_relatedPro[i]]])$estimate
}
outcome$padj <- p.adjust(outcome$p.value, method = "fdr")
outcome %>% filter(p.value < 0.05)
outcome %>% filter(padj < 0.05)

# scatter plots -----------------------------------------------------------------
plot_CD83_CD40 <- bulkRNAseq_baseline %>%
  ggplot(aes(x = CD83, y = CD40)) +
  geom_point(size =5, alpha = 0.5) +
  geom_smooth(method = "lm", se=FALSE) + stat_cor(size = 8) + 
  theme_classic()+
  theme(text = element_text(size = 24))

plot_CD83_CD80 <- bulkRNAseq_baseline %>%
  ggplot(aes(x = CD83, y = CD80)) +
  geom_point(size = 5, alpha = 0.5) +
  geom_smooth(method = "lm", se=FALSE) + stat_cor(size = 8) + 
  theme_classic()+
  theme(text = element_text(size = 24))


plot_CD83_CCR7 <- bulkRNAseq_baseline %>%
  ggplot(aes(x = CD83, y = CCR7)) +
  geom_point(size = 5, alpha = 0.5) +
  geom_smooth(method = "lm", se=FALSE) + stat_cor(size = 8) + 
  theme_classic()+
  theme(text = element_text(size = 24))

cowplot::plot_grid(plot_CD83_CD40, plot_CD83_CD80, plot_CD83_CCR7, nrow =1)

## save the plot 
png("output/cor_RNAseq_CD83_otherGenes.png", width = 1200, height = 480)
cowplot::plot_grid(plot_CD83_CD40, plot_CD83_CD80, plot_CD83_CCR7, nrow =1)
dev.off()

# correlation heatmap -----------------------------------------------------------------
library(rstatix)

corDat <- bulkRNAseq_baseline %>% select(c("CD83", CD83_relatedPro))
cor.mat <- corDat %>% cor_mat()
cor.mat %>% cor_get_pval()
cor.mat %>%
  cor_reorder() %>%
 # pull_lower_triangle() %>%
  cor_plot(label = TRUE, insignificant = "blank", 
           palette = get_palette(c("blue", "white", "red"), 20))

# with selected genes
corDat_v2 <- bulkRNAseq_baseline %>% select(c("CD83", "CD40", "CD80", "CCR7"))
cor.mat <- corDat_v2 %>% cor_mat()
cor.mat %>% cor_get_pval()
cor.mat %>% cor_get_pval()

png("output/cor_RNAseq_CD83_otherGenes_heatmap.png")

cor.mat %>%
  cor_reorder() %>%
  cor_plot(method = "color", type = "upper",
           label = TRUE, insignificant = "blank", 
           palette = get_palette(c("#3C5488FF", "white", "#E64B35FF"), 10),
           font.label = list(size = 1.8, color = "black", tl.cex = 5))

dev.off()
