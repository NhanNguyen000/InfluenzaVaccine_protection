rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(ggman)
library(qqman)


# load data, qtl T1 =================================================
df <- data.table::fread('processedDat/qtl/outcome_T1.csv') %>% set_colnames(c('SNP', 'gene', 'beta', 't', 'P'))
df <- cbind(df, colsplit(df$SNP, '%', names = c('inf', 'RSid')))
df <- cbind(df, colsplit(df$inf, ':', names = c('CHR', 'BP', 'REF', 'ALT')))

df$CHR <- as.numeric(df$CHR)

df %<>% dplyr::filter(!is.na(CHR))

png('processedDat/qtl/manhattans_T1.png', width = 10, height = 6, units = 'in', res=300)
ggpubr::ggarrange(
  ggmanhattan(df %>% dplyr::filter(gene == 'B_T1'), sparsify = F, significance=NA) +
    labs(title = 'Post-vaccination B') +
    theme_classic() +
    theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
    geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
    geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
    ylim(c(0, 10)),
  
  ggmanhattan(df %>% dplyr::filter(gene == 'H1N1_T1'), sparsify = F, significance=NA) +
    labs(title = 'Post-vaccination H1N1') +
    theme_classic() +
    theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
    geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
    geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
    ylim(c(0, 10)),
  
  ggmanhattan(df %>% dplyr::filter(gene == 'H3N2_T1'), sparsify = F, significance=NA) +
    labs(title = 'Post-vaccination H3N2') +
    theme_classic() +
    theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
    geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
    geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
    ylim(c(0, 10)), ncol = 1
)
dev.off()


# load data, qtl T4 =================================================
df <- data.table::fread('processedDat/qtl/outcome_T4.csv') %>% set_colnames(c('SNP', 'gene', 'beta', 't', 'P'))
df <- cbind(df, colsplit(df$SNP, '%', names = c('inf', 'RSid')))
df <- cbind(df, colsplit(df$inf, ':', names = c('CHR', 'BP', 'REF', 'ALT')))

df$CHR <- as.numeric(df$CHR)

df %<>% dplyr::filter(!is.na(CHR))

png('processedDat/qtl/manhattans_T4.png', width = 10, height = 6, units = 'in', res=300)
ggpubr::ggarrange(
  ggmanhattan(df %>% dplyr::filter(gene == 'B_T4'), sparsify = F, significance=NA) +
    labs(title = 'Post-vaccination B') +
    theme_classic() +
    theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
    geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
    geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
    ylim(c(0, 10)),
  
  ggmanhattan(df %>% dplyr::filter(gene == 'H1N1_T4'), sparsify = F, significance=NA) +
    labs(title = 'Post-vaccination H1N1') +
    theme_classic() +
    theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
    geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
    geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
    ylim(c(0, 10)),
  
  ggmanhattan(df %>% dplyr::filter(gene == 'H3N2_T4'), sparsify = F, significance=NA) +
    labs(title = 'Post-vaccination H3N2') +
    theme_classic() +
    theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
    geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
    geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
    ylim(c(0, 10)), ncol = 1
)
dev.off()

# 
# 
# # Regional association plots for H3N2 gwsig ones
# # rs112240980   46890948
# # rs4751339  132975733
# window <- 1e6/2 #1MB window
# 
# # rs112240980
# pos <- 46890948
# sub <- df %>% filter(CHR == 16) %>% filter(BP > pos - window) %>% filter(BP < pos + window)
# 
# png('topsnp_chr16.png', width = 2, height = 4, units = 'in', res=300)
# ggplot(sub) + 
#   geom_point(aes(x=BP, y=-log10(P))) +
#   geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
#   geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
#   theme_classic() +
#   labs(x = 'Position', y = '-log10(P)')
# dev.off()
# 
# 
# pos <- 132975733
# sub <- df %>% filter(CHR == 10) %>% filter(BP > pos - window) %>% filter(BP < pos + window)
# 
# png('topsnp_chr10.png', width = 2, height = 4, units = 'in', res=300)
# ggplot(sub) + 
#   geom_point(aes(x=BP, y=-log10(P))) +
#   geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
#   geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
#   theme_classic() +
#   labs(x = 'Position', y = '-log10(P)')
# dev.off()
# 
# 
# # Boxplots
# snp = 'rs112240980'
# 
# data.table::fread(cmd = 'cat genotype.csv | grep rs112240980') %>% tibble::column_to_rownames('V1') %>% t() -> dosages
# phenotypes <- read.table('phenotypes.csv')%>% t() %>% as.data.frame() %>% select(H3N2_T4)
# df <- cbind(dosages, phenotypes)
# 
# df$genotype <- 'TC'
# df$genotype <- ifelse(df$`16:46890948:T:C%rs112240980` > 1.5, 'CC', df$genotype)
# df$genotype <- ifelse(df$`16:46890948:T:C%rs112240980` < 0.5, 'TT', df$genotype)
# 
# pdf('boxplot_chr16.pdf', width = 2, height = 4)
# df %>% select(genotype, H3N2_T4) %>% melt() %>% 
#   ggplot() +
#   geom_boxplot(aes(x=genotype, y = value)) +
#   geom_point(aes(x=genotype, y=value), position = position_jitter(.1)) +
#   labs(x = 'Genotype', y = 'H3N2 (T4)', title = 'rs112240980\np=3.4e-08')
# dev.off()
