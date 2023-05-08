library(OlinkAnalyze)
library(openxlsx)
library(readxl)
library(tidyverse)

# iMED data --------------------------------------------------------------------
rawDat_iMEDprotein <- read_NPX('/vol/projects/CIIM/Influenza/iMED/proteomic/raw_data/20210222_Li_NPX_2021-05-18.csv') 

# ZirFlu data ------------------------------------------
rawDat_ZirFluprotein <- read_NPX(filename = "../proteomic/raw_data/20212645_Li_NPX_2022-02-02.csv")
ZirFlu_proteinPlate <- read_excel("../proteomic/raw_data/ZirrFlu plates final.xlsx",
                                  sheet = "Phenotypes")

ZirFlu_iMED_brideSamples <- ZirFlu_proteinPlate %>% 
  dplyr::select(patientID, probenID) %>% distinct() %>%
  dplyr::slice(which(grepl("human", patientID))) # 10 bridge samples with number and FR numbers

ZirFlu_iMED_brideSamples_v2 <- iMED$donorSample %>% filter(name %in% ZirFlu_proteinPlate$patientID) # 5 overlapped samples appear in the used iMED samples
intersect(ZirFlu_iMED_brideSamples_v2$name, ZirFlu_iMED_brideSamples$patientID) # 5 samples overlap

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
  pull(.) # 10 overlaped samples between 2 batches

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

# convert Olink normalized data to table
protein_normDat <- list()
protein_normDat$iMED <- protein_normOlink$normed_iMED %>% 
  filter(SampleID %in% iMED$donorSample$name) %>% 
  filter(Assay_Warning == 'PASS') %>% 
  filter(QC_Warning == 'PASS') %>% 
  filter(MissingFreq < 0.30) %>% 
  reshape2::dcast(data = ., SampleID ~ OlinkID, value.var = 'NPX') %>% 
  tibble::column_to_rownames(var = 'SampleID')

protein_normDat$ZirFlu <- protein_normOlink$normed_ZirFlu %>%
  mutate(SampleID = stringr::str_remove(SampleID, "^0+")) %>% 
  filter(SampleID %in% ZirFlu$donorSamples$probenID) %>% 
  filter(Assay_Warning == 'PASS') %>% 
  filter(QC_Warning == 'PASS') %>% 
  filter(MissingFreq < 0.30) %>% 
  mutate(SampleID = stringr::str_remove(SampleID, "^0+")) %>% 
  reshape2::dcast(data = ., SampleID ~ OlinkID, value.var = 'NPX') %>% 
  arrange(match(SampleID, ZirFlu$donorSamples$probenID)) %>% 
  tibble::column_to_rownames(var = 'SampleID')

# save data ----------------------------------------------------
# save(protein_normDat, file = "protein_normDat.RData")

# check data before and after normalization ------------------------------------
dat.preNorm <- ZirFlu$protein_dat %>% rownames_to_column("name") %>%
  full_join(iMED$protein_dat %>% rownames_to_column("name")) %>%
  column_to_rownames("name") %>% t()
  
boxplot(dat.preNorm)

dat.postNorm <- protein_normDat$ZirFlu %>% rownames_to_column("name") %>%
  full_join(protein_normDat$iMED %>% rownames_to_column("name")) %>%
  column_to_rownames("name") %>% t()

boxplot(dat.postNorm)

# check the normalization issue - unmatched normalization method sample ------------------------
# iMED normalization - "Intensity", ZirFlu = "Plate control" and "Intensity"
ZirFlu_plate.control_sample <- rawDat_ZirFluprotein_rename %>% 
  dplyr::select(SampleID, Normalization) %>% distinct() %>% 
  filter(Normalization == "Plate control")

ZirFlu_plate.control_sample2 <- ZirFlu_plate.control_sample  %>% 
  filter(SampleID %in% ZirFlu_proteinPlate$probenID)

ZirFlu_intensity.control_sample  <- rawDat_ZirFluprotein_rename %>% 
  dplyr::select(SampleID, Normalization) %>% distinct() %>% 
  filter(Normalization == "Intensity")

oneSample <- rawDat_ZirFluprotein_rename %>% filter(SampleID == "0308800416")

oneSample2 <-	rawDat_ZirFluprotein_rename %>% filter(Assay == "PNLIPRP2") # plate control in the normalization column
oneSample3 <-	rawDat_ZirFluprotein_rename %>% filter(Assay != "PNLIPRP2") # intensity in the normalization column
unique(oneSample3$Normalization) # only protein PNLIPRP2 use Plate control, the rest of proteins use intensity normalization


