rm(list = ls())

library(tidyverse)
library(OlinkAnalyze)
library(openxlsx)
library(readxl)

# load data =========================================================================
## protein data --------------------------------------------------------------------
# iMED data
rawDat_iMEDprotein <- read_NPX('/vol/projects/CIIM/Influenza/iMED/proteomic/raw_data/20210222_Li_NPX_2021-05-18.csv') 

# ZirFlu
rawDat_ZirFluprotein <- read_NPX(filename = "/vol/projects/CIIM/Influenza/ZirrFlu/proteomic/raw_data/20212645_Li_NPX_2022-02-02.csv")

## protein plate and cohort informations ----------------------------------------
ZirFlu_proteinPlate <- read_excel("/vol/projects/CIIM/Influenza/ZirrFlu/proteomic/raw_data/ZirrFlu plates final.xlsx",
                                  sheet = "Phenotypes")

ZirFlu_iMED_brideSamples <- ZirFlu_proteinPlate %>% 
  dplyr::select(patientID, probenID) %>% distinct() %>%
  dplyr::slice(which(grepl("human", patientID))) # 10 bridge samples with number and FR numbers

# iMED_metaCohort2 <- read.csv('/vol/projects/CIIM/Influenza/iMED/metadata/meta_cohort2.csv', row.names = 1)
# ZirFlu_iMED_brideSamples_v2 <- iMED_metaCohort2 %>% filter(name %in% ZirFlu_proteinPlate$patientID) # 5 overlapped samples appear in the used iMED samples (n = 200)
# intersect(ZirFlu_iMED_brideSamples_v2$name, ZirFlu_iMED_brideSamples$patientID) # 5 samples overlap

# rename the sample in ZirFlu to match the iMED sample name
rawDat_ZirFluprotein_rename <- rawDat_ZirFluprotein %>% 
  mutate(SampleID = str_replace_all(SampleID, "0198829641", "human-000434")) %>%
  mutate(SampleID = str_replace_all(SampleID, "0198829891", "human-000136")) %>%
  mutate(SampleID = str_replace_all(SampleID, "0183988146", "human-000082")) %>%
  mutate(SampleID = str_replace_all(SampleID, "0198795849", "human-000039")) %>%
  mutate(SampleID = str_replace_all(SampleID, "0193651919", "human-000516")) %>%
  mutate(SampleID = str_replace_all(SampleID, "FR05324258", "human-000058")) %>%
  mutate(SampleID = str_replace_all(SampleID, "FR05323560", "human-000629")) %>%
  mutate(SampleID = str_replace_all(SampleID, "FR05325648", "human-000059")) %>%
  mutate(SampleID = str_replace_all(SampleID, "FR03528024", "human-000045")) %>%
  mutate(SampleID = str_replace_all(SampleID, "FR05324060", "human-000070"))

# Find overlapping/bridge samples
overlap_samples <- intersect(rawDat_ZirFluprotein_rename$SampleID, rawDat_iMEDprotein$SampleID) %>% 
  data.frame() %>% 
  filter(!str_detect(., 'CONTROL_SAMPLE')) %>% #Remove control samples
  pull(.) # 10 overlapped samples between 2 batches

# normalization --------------------------------------------------------
# Perform Bridging normalization, iMED normalization - "Intensity", ZirFlu = "Plate control" and "Intensity"
norm_OlinkDat <- olink_normalization(df1 = rawDat_ZirFluprotein_rename, 
                                     df2 = rawDat_iMEDprotein, 
                                     overlapping_samples_df1 = overlap_samples,
                                     df1_project_nr = 'ZirFlu',
                                     df2_project_nr = 'iMED',
                                     reference_project = 'iMED')
# Olink norlamized data
protein_normOlink <- list()
protein_normOlink[["combined_2cohorts"]] <- norm_OlinkDat

protein_normOlink[["normed_iMED"]] <- norm_OlinkDat %>% 
  filter(Project == "iMED") %>% dplyr::select(-c(Project, Adj_factor))

protein_normOlink[["normed_ZirFlu"]] <- norm_OlinkDat %>% 
  filter(Project == "ZirFlu") %>% dplyr::select(-c(Project, Adj_factor))

# check the annotation 
iMED_proteinAnnot <- rawDat_iMEDprotein %>%
  distinct(Assay, .keep_all = TRUE) %>% select(OlinkID, UniProt, Assay) %>%
  arrange(OlinkID)
ZirFlu_proteinAnnot <- rawDat_ZirFluprotein_rename %>%
  distinct(Assay, .keep_all = TRUE) %>% select(OlinkID, UniProt, Assay) %>%
  arrange(OlinkID)

identical(iMED_proteinAnnot, ZirFlu_proteinAnnot) # TRUE, so the samme OlinkID and protein for both cohorts
proteinAnnot <- iMED_proteinAnnot

# convert Olink normalized data to table -------------------------
# # test if could directly convert to protein name
# dat_Olink <- protein_normOlink$normed_ZirFlu %>%
#   mutate(SampleID = stringr::str_remove(SampleID, "^0+")) %>% # to matched with sample ID in the metadata
#   filter(Assay_Warning == 'PASS') %>% 
#   filter(QC_Warning == 'PASS') %>% 
#   filter(MissingFreq < 0.30) %>% 
#   reshape2::dcast(data = ., SampleID ~ OlinkID, value.var = 'NPX') %>% 
#   slice(-c(grep("CONTROL", SampleID), grep("human", SampleID))) %>% # remove Control and brigde samples
#   tibble::column_to_rownames(var = 'SampleID')
# 
# dat_ProteinName <- protein_normOlink$normed_ZirFlu %>%
#   mutate(SampleID = stringr::str_remove(SampleID, "^0+")) %>% # to matched with sample ID in the metadata
#   filter(Assay_Warning == 'PASS') %>% 
#   filter(QC_Warning == 'PASS') %>% 
#   filter(MissingFreq < 0.30) %>% 
#   reshape2::dcast(data = ., SampleID ~ Assay, value.var = 'NPX') %>% 
#   slice(-c(grep("CONTROL", SampleID), grep("human", SampleID))) %>% # remove Control and brigde samples
#   tibble::column_to_rownames(var = 'SampleID') %>% 
#   # convert back to Olink ID to compare with dat_Olink
#   t() %>% as.data.frame %>%
#   rownames_to_column("Assay") %>% 
#   left_join(proteinAnnot %>% select(Assay, OlinkID)) %>%
#   select(-Assay) %>% column_to_rownames("OlinkID") %>% t() %>% as.data.frame %>%
#   select(names(dat_Olink))
# 
# identical(dat_Olink, dat_ProteinName) # TRUE, so could directly convert the Olink NPX data to protein name

protein_normDat <- list()
protein_normDat$iMED <- protein_normOlink$normed_iMED %>%
  filter(Assay_Warning == 'PASS') %>% 
  filter(QC_Warning == 'PASS') %>% 
  filter(MissingFreq < 0.30) %>% 
  reshape2::dcast(data = ., SampleID ~ Assay, value.var = 'NPX') %>% 
  slice(-grep("CONTROL", SampleID)) %>% # remove Control samples
  tibble::column_to_rownames(var = 'SampleID')

protein_normDat$ZirFlu <- protein_normOlink$normed_ZirFlu %>%
  # mutate(SampleID = stringr::str_remove(SampleID, "^0+")) %>% # to matched with sample ID in the metadata
  filter(Assay_Warning == 'PASS') %>% 
  filter(QC_Warning == 'PASS') %>% 
  filter(MissingFreq < 0.30) %>% 
  reshape2::dcast(data = ., SampleID ~ Assay, value.var = 'NPX') %>% 
  slice(-c(grep("CONTROL", SampleID), grep("human", SampleID))) %>% # remove Control and brigde samples
  tibble::column_to_rownames(var = 'SampleID')

# save data ----------------------------------------------------
save(protein_normDat, file = "protein_normDat.RData")

# # check data before and after normalization (need to clean the code)------------------------------------
# rm(iMED, ZirFlu) # remove the current object with the same name but different content
# load("../ZirFlu_NhanNguyen/ZirFlu.RData")
# load("../iMED_NhanNguyen/iMED.RData")
# 
# dat.preNorm <- ZirFlu$protein_dat %>% rownames_to_column("name") %>%
#   full_join(iMED$protein_dat %>% rownames_to_column("name")) %>%
#   column_to_rownames("name") %>% t()
# 
# boxplot(dat.preNorm)
# 
# dat.postNorm <- protein_normDat$ZirFlu %>% rownames_to_column("name") %>%
#   full_join(protein_normDat$iMED %>% rownames_to_column("name")) %>%
#   column_to_rownames("name") %>% t()
# 
# boxplot(dat.postNorm)
# 
# # check the normalization issue - unmatched normalization method sample ======================
# # iMED normalization - "Intensity", ZirFlu = "Plate control" and "Intensity"
# ZirFlu_plate.control_sample <- rawDat_ZirFluprotein_rename %>% 
#   dplyr::select(SampleID, Normalization) %>% distinct() %>% 
#   filter(Normalization == "Plate control")
# 
# ZirFlu_plate.control_sample2 <- ZirFlu_plate.control_sample  %>% 
#   filter(SampleID %in% ZirFlu_proteinPlate$probenID)
# 
# ZirFlu_intensity.control_sample  <- rawDat_ZirFluprotein_rename %>% 
#   dplyr::select(SampleID, Normalization) %>% distinct() %>% 
#   filter(Normalization == "Intensity")
# 
# oneSample <- rawDat_ZirFluprotein_rename %>% filter(SampleID == "0308800416")
# 
# oneSample2 <-	rawDat_ZirFluprotein_rename %>% filter(Assay == "PNLIPRP2") # plate control in the normalization column
# oneSample3 <-	rawDat_ZirFluprotein_rename %>% filter(Assay != "PNLIPRP2") # intensity in the normalization column
# unique(oneSample$Normalization)
# unique(oneSample2$Normalization)
# unique(oneSample3$Normalization) # only protein PNLIPRP2 use Plate control, the rest of proteins use intensity normalization
# # Answer from the Olink company: "The PNLIPRP2 is IPC normalized because of the bimodal distribution of this assay. You don???t need to worry about that when you are performing the bridge normalization"

