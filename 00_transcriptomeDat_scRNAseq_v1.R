library(Seurat)
# one data file
Raw_data <- Read10X(data.dir = "data/matrix/I0-1.matrix")
dat_temp <- CreateSeuratObject(counts = Raw_data)

# all data files
sampleList <- list.files("data/matrix")
seauratObj <- list()
for (sampleName in sampleList) {
  rawDat_temp <- Read10X(data.dir = paste0("data/matrix/", sampleName))
  seauratObj[[sampleName]] <- CreateSeuratObject(counts = rawDat_temp)
}

dat <- merge(x = seauratObj$`I0-1.matrix`,
             y = seauratObj[2:6],
             add.cell.ids = names(seauratObj)[1:6])
dat <- RunTFIDF(dat)
dat <- FindTopFeatures(dat, min.cutoff = 20)
dat <- RunSVD(dat)
dat <- RunUMAP(dat, dims = 2:50, reduction = 'lsi')

#dat <- NormalizeData(object = dat)
dat <- FindVariableFeatures(object = dat)
dat <- ScaleData(object = dat)
dat <- RunPCA(object = dat)
dat <- FindNeighbors(object = dat)
dat <- FindClusters(object = dat)

dat <- RunTSNE(object = dat)
DimPlot(object = dat, reduction = "tsne")


dat <- RunUMAP(object = dat, dims = 2:50, reduction = 'lsi')