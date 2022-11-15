############################################################################
##### RNA-seq whole head swarming and non-swarming males Aedes aegypti #####
#####                   GENE LEVEL ANALYSIS                            #####
############################################################################



### Libraries and working directory ----

# Working directory
setwd("~/Desktop/Current_work/Brasilian_samples/")
# Libraries
library(DESeq2)          # Bioconductor v1.36.0
library(ggplot2)         # CRAN v3.3.6
library(ggrepel)         # CRAN v0.9.1
library(RColorBrewer)    # CRAN v1.1-3
library(pheatmap)        # CRAN v1.0.12
library(PoiClaClu)       # CRAN v1.0.2.1
library(apeglm)          # Bioconductor v1.18.0
library(tximportData)    # Bioconductor v1.24.0
library(tximport)        # Bioconductor v1.24.0
library(factoextra)      # CRAN v1.0.7
library(tidyverse)       # CRAN v1.3.1
library(dplyr)           # CRAN v1.0.9
library(readxl)          # CRAN v1.4.1
library(MODA)            # Bioconductor v1.22.0
library(mgsa)            # Bioconductor v1.44.0 
library(clusterProfiler) # Bioconductor v4.4.4



### Files and experimental design ----
samples=c("S1","S2","S3", "S4","N1","N2", "N3", "N4")
behaviour=c(rep("Swarming",4), rep("NonSwarming",4))
seq=c("s11", "s12", "s13", "s14", "s21", "s22", "s23", "s24")
path=list.files("./counts/", pattern="*.genes.results")
sample_info=data.frame(samples, behaviour, path)
rownames(sample_info) <- sample_info$sample
# Get all paths from wd
files_paths = file.path("./counts", paste0(sample_info$path))
# Put samples names as path name
names(files_paths) <- paste0(sample_info$sample)
# import the counts
txi.rsem <- tximport(files_paths, type = "rsem", txIn = FALSE, txOut = FALSE)
# remove null values
txi.rsem$length[txi.rsem$length == 0] <- 1
# Make deseq2 object
dds <- DESeqDataSetFromTximport(txi.rsem,
                                colData = sample_info,
                                design = ~ behaviour)
# remove features with low counts
dds = dds[ rowMeans(counts(dds)) > 10, ] 
nrow(dds) # 19804 => 12028
# rlog transform counts (by average Tx length and correcting for library size)
rld = rlog(dds, blind=FALSE)




#### Descriptive exploration ----

### PCA
x_pca = prcomp(t(assay(rld)), scale. = T)
# Sum pca
summary(x_pca)
# Exploring axes explanation
fviz_screeplot(x_pca, addlabels = T) + ggtitle('')
fviz_pca_ind(x_pca, asp=1, pointshape=20, habillage = behaviour) + scale_color_discrete() + ggtitle('PCA transcriptomes')

### Explore PCA results with dataframe
pca_df <- x_pca$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sample,
             group = behaviour)
pca_pivot <- pivot_longer(pca_df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)
ggplot(pca_pivot) +
  aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot") +
  theme_bw() +
  coord_flip()


#### Check normalisation parameters DESeq ----

# two 1st samples plotted against each other to check consistency (for rlog and log2) 
par( mfrow = c( 1, 2 ) )
dds = estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,c(1,2)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
plot(assay(rld)[,c(1,2)],
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
par( mfrow = c( 1, 1 ) )

# estimate size factors = normalize for dispersion
dds = DESeq2::estimateSizeFactors(dds)
dds = estimateDispersions(dds)
plotDispEsts(dds, xlab= "Mean of Normalised Counts",
             ylab= "Dispersion", cex=1.0, cex.lab=1.45, cex.axis=1.45)

# check sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         show_colnames = T,
         fontsize = 18)

# check sample distances using the poisson method
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep=" - " )
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         show_colnames = T,
         fontsize = 18)



#### Differential Gene Expression analysis ----
# differential expression (used Benjamini-Hochberg adjustment and a p-val of <0.1)
dds = DESeq2::DESeq(dds, parallel=TRUE)
resultsNames(dds)
# Swarming = +Log2FC and non-swarming = -Log2FC
res=results(dds, name="behaviour_Swarming_vs_NonSwarming")
summary(res)

#distribution of coefficents of the model
plotMA(res, ylim=c(-5,5),cex=1.0, cex.lab=1.45, 
       cex.axis=1.45, xlab='Mean of Normalised Counts',
       ylab='Log2 Fold Change')
# plot of p-vals excluding genes with very small counts
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", xlab= "P-value", ylab="Frequency",
     cex.axis=1.45, cex.lab=1.45, main="")

# Order the results by fold change to see the biggest changing genes
res_ordered=res[order(res$log2FoldChange),]
head(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_ordered$gene<-rownames(res_ordered)
rownames(res_ordered)<-c()
write.table(as.data.frame(res_ordered), file="diff_exp_output_log2FC_all_genes.txt",sep="\t", quote = F, col.names = T, row.names = F)

#out of 12420 with nonzero total read count
res_significant <- subset(res, log2FoldChange > 0.75 | log2FoldChange < -0.75)
res_significant <- subset(res_significant, padj < 0.05)
summary(res_significant) # 1 gene differentially express.... = Pdp (Pyruvate dehydrogenase phosphatase) overexpress in Swarming male's heads

### Corresponding boxplot
n=1
selGenes = head(rownames(res)[order(res$padj)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, plotCounts(dds, gene=gene, intgroup=c("behaviour"), returnData=TRUE))))
ggplot(data, aes(x=behaviour, y=count, fill=behaviour)) + 
  geom_boxplot(outlier.color="black", position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("behaviour") + ylab("Normalised read count") + 
  scale_y_log10() + 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))+
  scale_fill_manual("", breaks=c("Swarming","NonSwarming"), values = c("#DDCC77","#44AA99"))




#### Differential expression Module analysis -----

# Generate modules in the 2 conditions using a community detection method (Louvain algorithm (Blondel et al. 2008))
intModules1=WeightedModulePartitionLouvain(t(dds@assays@data$counts), '~/Desktop/Current_work/Brasilian_samples/MODA3/', 'S', rownames(dds@assays@data$counts))
intModules2=WeightedModulePartitionLouvain(t(dds@assays@data$counts), '~/Desktop/Current_work/Brasilian_samples/MODA3/', 'N', rownames(dds@assays@data$counts))
# Compare the 2 sets of modules
CompareAllNets('~/Desktop/Current_work/Brasilian_samples/MODA3/',intModules1,'S', intModules2, 'N', specificTheta=0.1,conservedTheta=0.1)
# Get the similarity matrix
JaccardMatrix <- comparemodulestwonets('~/Desktop/Current_work/Brasilian_samples/MODA3/',intModules1,intModules2, 'DenseModuleGene_S', 'DenseModuleGene_N')
# Select non conserved modules
nonconserved_modules = which(rowSums(JaccardMatrix)<(1/3))
# Select the genes from non conserved modules
nonconserved_genes = unlist(lapply(nonconserved_modules, function(x) read_csv(paste0("MODA3/DenseModuleGene_N_", x, ".txt"), col_names = F)))
# GO enrichment analysis on genes in the non conserved modules 
# Get the universe and syntax
aedes_GO=readGAF("~/Desktop/Current_work/Brasilian_samples/VectorBase-57_AaegyptiLVP_AGWG_GO.gaf")
aedesgo_annotation=as.data.frame(cbind(row.names(aedes_GO@setAnnotations), aedes_GO@setAnnotations$term))
colnames(aedesgo_annotation)=c('GOterm','function')
# Do GO enrichment analysis
GOen=enricher(nonconserved_genes, TERM2GENE=cbind(GO_file$X5, GO_file$X2), TERM2NAME=aedesgo_annotation)
# Extract significantly enrich GO terms (3 enrich GO terms without big insights)
res=as.data.frame(GOen@result)
res=res[res$p.adjust<0.05,]
revigo=res[,c(1,6)]

### Test module analysis with other module delimitation method
# Generate modules in the 2 conditions using spectral clustering (White and Smyth 2005)
intModules1=WeightedModulePartitionSpectral(t(dds@assays@data$counts), '~/Desktop/Current_work/Brasilian_samples/MODA4/', 'S', rownames(dds@assays@data$counts), k=10)
intModules2=WeightedModulePartitionSpectral(t(dds@assays@data$counts), '~/Desktop/Current_work/Brasilian_samples/MODA4/', 'N', rownames(dds@assays@data$counts), k=10)
# Compare the 2 sets of modules
CompareAllNets('~/Desktop/Current_work/Brasilian_samples/MODA4/',intModules1,'S', intModules2, 'N', specificTheta=0.1,conservedTheta=0.1)
# Get the similarity matrix
JaccardMatrix <- comparemodulestwonets('~/Desktop/Current_work/Brasilian_samples/MODA4/',intModules1,intModules2, 'DenseModuleGene_S', 'DenseModuleGene_N')
# Select non conserved modules
nonconserved_modules = which(rowSums(JaccardMatrix)<(1/3))
# Select the genes from non conserved modules
nonconserved_genes = unlist(lapply(nonconserved_modules, function(x) read_csv(paste0("MODA4/DenseModuleGene_N_", x, ".txt"), col_names = F)))








