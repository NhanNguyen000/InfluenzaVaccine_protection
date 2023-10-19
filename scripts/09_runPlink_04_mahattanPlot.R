# Mahattan plot based on Martijn code --------------------------------
library(dplyr)
library(ggman)
library(reshape2)

setwd('/vol/projects/CIIM/Influenza/ZirrFlu/genetics_abqtl/out/mapping/') # folder with the last updated result uing QTL Javi pipeline

# B strain
data.table::fread('B_T4.txt') -> df
df <- cbind(df, colsplit(df$SNP, '%', c('info', 'rsID')))
df <- cbind(df, colsplit(df$info, ':', c('CHR', 'BP', 'REF', 'ALT')))
df$pvalue <- df$`p-value`


png('b_manhattan.png', width = 8, height = 3, units = 'in', res=300)
ggman(df, sparsify = T, significance=NA) +
  labs(title = 'Post-vaccination B titers') +
  theme_classic() +
  theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
  geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
  geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
  ylim(c(0, 10))
dev.off()

# H1N1 strain
data.table::fread('H1N1_T4.txt') -> df
df <- cbind(df, colsplit(df$SNP, '%', c('info', 'rsID')))
df <- cbind(df, colsplit(df$info, ':', c('CHR', 'BP', 'REF', 'ALT')))
df$pvalue <- df$`p-value`


png('h1n1_manhattan.png', width = 8, height = 3, units = 'in', res=300)
ggman(df, sparsify = T, significance=NA) +
  labs(title = 'Post-vaccination H1N1 titers') +
  theme_classic() +
  theme(plot.title = element_text(hjust = .5), legend.position = 'none') +
  geom_hline(yintercept = -log10(1e-5), col = 'lightblue', lty=2) +
  geom_hline(yintercept = -log10(5e-8), col = 'red', lty=2) +
  ylim(c(0, 10))
dev.off()

# get the top sig. SNP ----------------------------------
strains <- c("H1N1", "H3N2", "B")

snps <- list()
for (strain in strains) {
  snps[[strain]] <- read.table(
    paste0("/vol/projects/CIIM/Influenza/ZirrFlu/genetics_abqtl/out/mapping/", strain, "_T4.txt"), 
    header = TRUE)
}

sigSnps <- snps %>% lapply(function(x) x %>% filter(p.value < 5e-7))
# For H3N2, I got 4 top sig. SNP, and Maritjn only report the 1st, 2nd, and 4th


intersect(sigSnps$H1N1$SNP, sigSnps$H3N2$SNP) # no overlap
intersect(sigSnps$H1N1$SNP, sigSnps$B$SNP) # no overlap
intersect(sigSnps$H3N2$SNP, sigSnps$B$SNP) # no overlap

