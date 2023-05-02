#library(edgeR)
library(limma)
group
x1 <- as.factor(sample(c("case", "control"), 10, replace=TRUE))
x1
model.matrix(~x1)
design <- model.matrix(~ $er)
fit <- lmFit(vdx, design)
fit <- eBayes(fit)
topTable(fit)

metadat <- cohorts_dat$donorSample_all %>% filter(time == "T1") %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  mutate(sex = substring(sex, 1, 1)) %>%
  filter(disease == "healthy") %>%
  filter(category != "Other")

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat$name) %>% 
  arrange(factor(name, levels = metadat$name)) %>%
  column_to_rownames("name") %>% t()

identical(metadat$name, colnames(inputDat))

designTable <- model.matrix(~ + sex + age + category, metadat)
a <- lmFit(inputDat, 
      design = designTable) %>% eBayes() 

protein_limma.t.stat <- a$t[, "categoryTR"]

# run fgsea package ----------------
library(fgsea)
protein_limma.t.stat_geneID <- as.data.frame(protein_limma.t.stat) %>% 
  rownames_to_column("OlinkID") %>% left_join(protein_Entrez)

ranked_proteins <- protein_limma.t.stat_geneID$protein_limma.t.stat
names(ranked_proteins) <- protein_limma.t.stat_geneID$ENTREZID

fgseaRes_TRvsNR <- fgsea(hallmarkSets, ranked_proteins ) 
fgseaRes_TRvsNR_inflam <- fgseaRes_TRvsNR %>% filter(pathway %in% hallmarkSubsets_v2)
leadingProteins_TRvsNR <- fgseaRes_TRvsNR_inflam$leadingEdge %>% unlist() %>% unique()

leadingProteins_TRvsNR_dat <- protein_Entrez %>% 
  filter(ENTREZID %in% leadingProteins_TRvsNR) # 27 proteins
