rm(list = ls())
library(tidyverse)
library(caret)
library(reshape2)

# load data ---------------------------------------------------
load("cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# protein data --------------------------------------------------
proteinDat_metadata <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  bind_rows(.id = "group") %>% select(group, name) %>% 
  inner_join(metadata_healthy)

proteinDat <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  bind_rows(.id = "group") %>% filter(name %in% proteinDat_metadata$name) %>%
  select(-group) %>% column_to_rownames("name")

protein_pca <- predict(preProcess(x=proteinDat, method = "knnImpute"), # imput the NA value
                       proteinDat) %>% prcomp()
summary(protein_pca)$importance[1:3, 1:5]

# correlations PCs - variables (adapt from Manoj codes) --------------------------------------------------
## calculate the correlation PCs - variable -----------------------
phenoDat <- metadata_healthy %>% 
  select(name, age, sex, season, responder,
         matches("abFC_log2|T1_log2")) %>% 
  column_to_rownames("name")
head(phenoDat)

n_princomps <- 30 # PCs 1-30
princomps <- protein_pca$x[, 1:n_princomps] 
corr.pvalues <- matrix(nrow = ncol(princomps), 
                       ncol = ncol(phenoDat), 
                       dimnames = list(colnames(princomps), colnames(phenoDat)))

for (var1 in colnames(phenoDat)) {
  for ( pc in colnames(princomps)) {
    print(pc)
    print(var1)
    #res <- cor.test(princomps[, pc], pheno2[, var1])
    #pval <- -1 * log10(res$p.value)
    #pval <- res$p.value
    #corr.pvalues[pc, var1] <- pval
    
    a<- lm ( princomps[, pc] ~ phenoDat[, var1] )
    b<- data.frame(summary(a)$coefficients[2,])
    pvali<-b[4,]
    corr.pvalues[pc, var1] <- pvali
  }
}

## prepare the plot data ----------------------------------------------
df.res <- melt(corr.pvalues)
colnames(df.res)=c("PC","cov","pval")
table(df.res)

#df.res$significance <- cut(df.res$pval, breaks=c(-Inf, 0.01, 0.05, Inf), label=c( "P<0.01", "P<0.05", "NS")) 
df.res$significance <- cut(df.res$pval, breaks=c(-Inf, 1e-10,1e-5, 0.01, 0.05, Inf), label=c("P<1e-10","P<1e-5", "P<0.01", "P<0.05", "NS")) 
PC <- colnames(pcs$x[, 1:n_princomps])
cov <- colnames(phenoDat)
df.res$PC <- factor(df.res$PC,levels=PC)
df.res$cov <- factor(df.res$cov,levels = cov)

### calculate the accumulated variance ---------------------------------
pca_var <- protein_pca$sdev^2
pca_var_per <- (pca_var / sum(pca_var) * 100)
pca_var_per[1:30]
acumullated_variance <- sum(pca_var_per[1:30])

## plot the figure ------------------------
ggplot(data = df.res, aes(x = PC, y = cov, fill = significance)) +
  geom_tile() + scale_fill_manual(values = c("#E31A1C","#FC4E2A","#FEB24C","#FFF3B2","white")) +
  ylab("Covariates") + xlab("Principal component")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))+
  theme(legend.position = "bottom") +
  ggtitle(paste0("Acumullated variance: ", round(acumullated_variance, 2), "%"))



