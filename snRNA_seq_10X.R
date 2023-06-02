############################################################################
#####                  Single nuclei 10X Aedes aegypti                 #####
############################################################################

# Cui, Y., Behura, S.K. & Franz, A.W.E. Cellular diversity and gene expression profiles in the male and female brain of Aedes aegypti. 
# BMC Genomics 23, 119 (2022). https://doi.org/10.1186/s12864-022-08327-9




### Libraries and working directory ----

# Working directory
setwd("~/Desktop/data/0_Current_work/Brasilian_samples/singlenuclei/")
# Libraries
library(Seurat)  # CRAN v4.3.0 
library(ggplot2) # CRAN v3.4.2 


### Import Male brain data ----

# GSM4878383
list.files('./male/')
malebrains = Read10X(data.dir = './male/')
malebrains_seurat_object = CreateSeuratObject(counts = malebrains)
# Filtered low quality nuclei (<500 count & >4000 counts)
pbmc_male <- subset(malebrains_seurat_object, subset = nFeature_RNA > 500 & nFeature_RNA < 4000)




### Import Female brain data ----

# GSM4878384
list.files('./female/')
femalebrains = Read10X(data.dir = './female/')
femalebrains_seurat_object = CreateSeuratObject(counts = femalebrains)
# Filtered low quality nuclei (<500 count & >4000 counts)
pbmc_female <- subset(femalebrains_seurat_object, subset = nFeature_RNA > 500 & nFeature_RNA < 4000)




### Merge the 2 sexes and normalise ----

pbmc.combined <- merge(pbmc_male, y = pbmc_female, add.cell.ids = c("Male", "Female"), project = "test")
pbmc.combined

# Normalisation 
pbmc.combined <- NormalizeData(pbmc.combined, normalization.method = "LogNormalize", scale.factor = 10000)

# Find genes of interest
pbmc.combined <- FindVariableFeatures(pbmc.combined, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(pbmc.combined)
# Scaling data
all.genes <- rownames(pbmc.combined)
pbmc.combined <- ScaleData(pbmc.combined, features = all.genes)

# PCA
pbmc.combined <- RunPCA(pbmc.combined, features = VariableFeatures(object = pbmc.combined))
DimPlot(pbmc.combined, reduction = "pca")


# Check for dimension dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc.combined <- JackStraw(pbmc.combined, num.replicate = 100)
pbmc.combined <- ScoreJackStraw(pbmc.combined, dims = 1:20)
JackStrawPlot(pbmc.combined, dims = 1:20)
ElbowPlot(pbmc.combined)
# seems to be ~9-10

pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:20)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 1)




#### Linear tSNE ----
pbmc.combined <- RunTSNE(pbmc.combined)

TSNEPlot(pbmc.combined)

FeaturePlot(pbmc.combined, features = c('AAEL010711', 'AAEL002879', 'AAEL005950', 'AAEL025907', "AAEL019804", "AAEL004731", "AAEL009895", "AAEL009895", "AAEL009895", "AAEL018348", "AAEL020244", "AAEL018306", "AAEL021054", "AAEL004098", "AAEL010483", "AAEL006283", "AAEL008507", "AAEL009211", "AAEL022180", "AAEL000181", "AAEL019711", "AAEL008603", "AAEL022806", "AAEL005478", "AAEL027400"))
FeaturePlot(pbmc.combined, features = c('AAEL010711', 'AAEL002879', 'AAEL005950', "AAEL019804", "AAEL004731", "AAEL009895", "AAEL009895", "AAEL009895", "AAEL018348", "AAEL018306", "AAEL021054", "AAEL004098", "AAEL010483", "AAEL006283", "AAEL008507", "AAEL009211", "AAEL022180", "AAEL000181", "AAEL019711", "AAEL008603", "AAEL022806", "AAEL005478", "AAEL027400"))





#### Non-linear tSNE ----

pbmc.combined=RunTSNE(pbmc.combined, dims=1:20)
TSNEPlot(pbmc.combined)


FeaturePlot(pbmc.combined, features = c("AAEL010711", "AAEL027189", "AAEL002879", "AAEL005950", "AAEL008603", "AAEL025907", "AAEL019804", "AAEL019898", "AAEL004731", "AAEL022806", "AAEL009895", "AAEL018348", "AAEL020244", "AAEL025228", "AAEL018306", "AAEL021054", "AAEL005478", "AAEL027790", "AAEL004098", "AAEL027400", "AAEL010483", "AAEL006283", "AAEL008507", "AAEL009211", "AAEL022180", "AAEL000181", "AAEL019711"))
FeaturePlot(pbmc.combined, features = c("AAEL010711", "AAEL027189", "AAEL002879", "AAEL005950", "AAEL008603", "AAEL019804", "AAEL019898", "AAEL004731", "AAEL022806", "AAEL009895", "AAEL018348", "AAEL025228", "AAEL018306", "AAEL021054", "AAEL005478", "AAEL027790", "AAEL004098", "AAEL027400", "AAEL010483", "AAEL006283"))
FeaturePlot(pbmc.combined, features = c("AAEL008507", "AAEL009211", "AAEL022180", "AAEL000181", "AAEL019711"))




###











