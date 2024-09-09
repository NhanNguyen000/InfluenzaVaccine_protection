library(dplyr)
library(reshape2)
library(ggplot2)

setwd('/vol/projects/CIIM/Influenza/ZirrFlu/genetics_abqtl/out/mapping/')

data.table::fread('H3N2_T4.txt') -> df
df <- cbind(df, colsplit(df$SNP, '%', c('info', 'rsID')))
df <- cbind(df, colsplit(df$info, ':', c('CHR', 'BP', 'REF', 'ALT')))
df$P <- df$`p-value`


png('h3n2_manhattan.png', width = 8, height = 3, units = 'in', res=300)
ggman::ggmanhattan(df, sparsify = T, significance=NA) +
  labs(title = 'Post-vaccination H3N2 titers') +
  theme_classic() +
  theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
  geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
  geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
  ylim(c(0, 10))
dev.off()



# B
data.table::fread('B_T4.txt') -> df
df <- cbind(df, colsplit(df$SNP, '%', c('info', 'rsID')))
df <- cbind(df, colsplit(df$info, ':', c('CHR', 'BP', 'REF', 'ALT')))
df$P <- df$`p-value`


png('B_manhattan.png', width = 8, height = 3, units = 'in', res=300)
ggman::ggmanhattan(df, sparsify = T, significance=NA) +
  labs(title = 'Post-vaccination B titers') +
  theme_classic() +
  theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
  geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
  geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
  ylim(c(0, 10))
dev.off()



# H1N1
data.table::fread('H1N1_T4.txt') -> df
df <- cbind(df, colsplit(df$SNP, '%', c('info', 'rsID')))
df <- cbind(df, colsplit(df$info, ':', c('CHR', 'BP', 'REF', 'ALT')))
df$P <- df$`p-value`


png('h1n1_manhattan.png', width = 8, height = 3, units = 'in', res=300)
ggman::ggmanhattan(df, sparsify = T, significance=NA) +
  labs(title = 'Post-vaccination H1N1 titers') +
  theme_classic() +
  theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
  geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
  geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
  ylim(c(0, 10))
dev.off()
