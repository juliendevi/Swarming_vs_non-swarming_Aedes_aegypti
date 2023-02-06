############################################################################
##### RNA-seq whole head swarming and non-swarming males Aedes aegypti #####
#####                TRANSCRIPT LEVEL ANALYSIS                         #####
############################################################################




### Libraries and working directory ----

# Working directory
setwd("~/Desktop/Current_work/Brasilian_samples/kallisto/")
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
library(readr)           # CRAN v2.1.2 
library(mgsa)            # Bioconductor v1.44.0 
library(clusterProfiler) # Bioconductor v4.4.4 
library(enrichplot)      # Bioconductor v1.16.2 
library(MODA)            # Bioconductor v1.22.0




### Import experiment information and transcript counts from kallisto ----

# import experimental scheme
exp = read_excel("~/Desktop/Current_work/Brasilian_samples/experiment.xlsx", col_names = T)
# Make experiment design matrix
group = as.factor(exp$group)
design = model.matrix(~0 + group)
colnames(design)=levels(group)
# find file path
path = file.path(exp$sequence, "abundance.tsv")
# check path exist
all(file.exists(path))
exp$path=path

# import kallisto counts
kallisto_out = tximport(path, type = 'kallisto', txOut = T, geneIdCol = T)

# Make deseq2 object
dds <- DESeqDataSetFromTximport(kallisto_out,
                                colData = exp,
                                design = ~ group)
# remove features with low counts
dds = dds[ rowMeans(counts(dds)) > 10, ] 
nrow(dds) # 34477/25059

# rlog transform counts (by average Tx length and correcting for library size)
rld = rlog(dds, blind=FALSE)



#### Descriptive exploration ----

### PCA
x_pca = prcomp(t(assay(rld)), scale. = T)
# Sum pca
summary(x_pca)
# Exploring axes explanation
fviz_screeplot(x_pca, addlabels = T) + ggtitle('')
fviz_pca_ind(x_pca, asp=1, pointshape=20, habillage = group) + scale_color_discrete() + ggtitle('PCA transcriptomes')

### Explore PCA results with dataframe
# Arrange the df in tibble
pca_df <- x_pca$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = exp$sample,
             group = group)
# Pivot the dataset
pca_pivot <- pivot_longer(pca_df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)
# Show the plot
ggplot(pca_pivot) +
  aes(x=sample, y=loadings, fill=group) + 
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
res=results(dds, name="group_Swarming_vs_NonSwarming")
summary(res)

#distribution of coefficents of the model
plotMA(res, ylim=c(-5,5),cex=1.0, cex.lab=1.45, 
       cex.axis=1.45, xlab='Mean of Normalised Counts',
       ylab='Log2 Fold Change')

# plot of p-vals excluding genes with very small counts
hist(res$padj[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", xlab= "P-value", ylab="Frequency",
     cex.axis=1.45, cex.lab=1.45, main="")

# Order the results by fold change to see the biggest changing genes
res_ordered=res[order(res$log2FoldChange),]
head(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_ordered$gene<-rownames(res_ordered)
rownames(res_ordered)<-c()
write.table(as.data.frame(res_ordered), file="diff_exp_output_log2FC_all_genes.txt",
            sep="\t", quote = F, col.names = T, row.names = F)

#out of 25059 with nonzero total read count
res_significant<- subset(res, log2FoldChange > 2 | log2FoldChange < -2)
res_significant <- subset(res_significant, padj < 0.05)
nrow(res_significant) #17 total differentially expressed genes
nrow(res_significant[res_significant$log2FoldChange > 0.75,]) #9 upregulated in males
nrow(res_significant[res_significant$log2FoldChange < -0.75,]) #8 upregulated in females

#out of 25059 with nonzero total read count
res_significant<- subset(res, log2FoldChange > 0.75 | log2FoldChange < -0.75)
res_significant <- subset(res_significant, padj < 0.05)
nrow(res_significant) #27 total differentially expressed genes
nrow(res_significant[res_significant$log2FoldChange > 0.75,]) #12 upregulated in males
nrow(res_significant[res_significant$log2FoldChange < -0.75,]) #15 upregulated in females

#out of 25059 with nonzero total read count
res_significant<- subset(res, log2FoldChange > 0.75 | log2FoldChange < -0.75)
res_significant <- subset(res_significant, padj < 0.01)
nrow(res_significant) #14 total differentially expressed genes
nrow(res_significant[res_significant$log2FoldChange > 0.75,]) #6 upregulated in males
nrow(res_significant[res_significant$log2FoldChange < -0.75,]) #8 upregulated in females




#### Visualisation differential expression ----

#### Heatmap

# Set number of transcripts to look at
n=27
topdiff = head(c(1:nrow(res))[order(res$padj)],n)

# Set colours
library(MetBrewer)       # [github::BlakeRMills/MetBrewer] v0.2.0
a=met.brewer('Egypt', n=4, 'discrete')
my_colors = list(
  behaviour = c(Swarming = a[3] , NonSwarming =a[1]))

# Take the counts 
mat = assay(rld)[ topdiff, ]
mat = mat - rowMeans(mat)
df2 = as.data.frame(colData(rld)[,c("behaviour"),drop=FALSE])
colnames(df2)<-c("behaviour")

# Take the annotation
annotgenes = read_excel('~/Desktop/Brasilian_samples/graphs_test/significant_genes.xlsx', sheet = 'Feuil2', col_names = F)
rownames(mat)=annotgenes$...6

# Generate the heatmap
pheatmap(mat, annotation_col=df2,
         show_rownames = T,annotation_colors = my_colors,
         col=met.brewer("Manet",n=100, type="continuous"),
         legend_breaks = c(-0.2, 0, 0.2,0.4, max(mat)),
         legend_labels = c("-0.2", "0", "0.2","0.4", "z-score"),
         fontsize = 16, annotation_legend = F, angle_col = 0, annotation_names_col = F, border_color = "white")

#### Boxplot

# Set number of transcripts to look at
n=27
selGenes = head(rownames(res)[order(res$padj)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, plotCounts(dds, gene=gene, intgroup=c("behaviour"), returnData=TRUE))))

# Generate the boxplot
ggplot(data, aes(x=behaviour, y=count, fill=behaviour)) + 
  geom_boxplot(outlier.color="black", position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Behaviour") + ylab("Normalised read count") + 
  scale_y_log10() + 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22),
        legend.position = 'none')+
  scale_fill_manual("", breaks=c("Swarming","NonSwarming"),
                    values = c(a[3],a[1]))

#### Volcano plot
res_df <- as.data.frame(res)
res_df$gene <- row.names(res_df)
res_significant_df <- as.data.frame(res_significant)
res_significant_df$gene <- row.names(res_significant_df)
res_df$sig <- "no"
res_df$sig[res_df$gene %in% res_significant_df$gene] <- "yes"

# Generate the plot
v=ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig, text = paste("Symbol:", gene))) + 
  geom_point(size=1.5)+
  scale_colour_manual("", breaks=c("no","yes"),
                      values = c("black","red"))+
  xlim(-10,10)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0.75, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -0.75, linetype="longdash", colour="#2C467A", size=1) +
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = "none")
v
# For interactive plot: 
#ggplotly(v)



#### Export list significantly differentially expressed transcripts ----

significant_genes = which(res_df$sig=='yes')
df_significant_genes = as.data.frame(rbind(res_df[significant_genes,]))
write.table(as.data.frame(df_significant_genes), file="significant_genes_like_paper.txt", sep="\t", quote = F, col.names = T, row.names = F)


### GO enrichment analysis --------

# Put annotations back to gene level
a=res_significant@rownames
b=sapply(a, function(x){
  B=strsplit(x, "[[:punct:]]")[[1]]
  return(B[1])
}
)
res_significant@rownames=b

c=res@rownames
b=sapply(c, function(x){
  B=strsplit(x, "[[:punct:]]")[[1]]
  return(B[1])
}
)
res@rownames=b

# Import GO term annotation file
aedes_GO=readGAF("~/Desktop/Current_work/Brasilian_samples/VectorBase-57_AaegyptiLVP_AGWG_GO.gaf")
aedesgo_annotation=as.data.frame(cbind(row.names(aedes_GO@setAnnotations), aedes_GO@setAnnotations$term))
colnames(aedesgo_annotation)=c('GOterm','function')
GO_file = read_delim("~/Desktop/Current_work/Brasilian_samples/VectorBase-57_AaegyptiLVP_AGWG_GO.gaf", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
# Get significant genes
GeneList = as.data.frame(cbind(row.names(res_significant),res_significant$log2FoldChange))
# GO enrichment analysis
a=enricher(GeneList[,1], TERM2GENE=cbind(GO_file$X5, GO_file$X2), TERM2NAME=aedesgo_annotation)
res=as.data.frame(a@result)
res=res[res$p.adjust<0.05,] #=0





#### Module analysis -----

# Generate modules in the 2 conditions using a community detection method (Louvain algorithm (Blondel et al. 2008))
intModules1=WeightedModulePartitionLouvain(t(dds@assays@data$counts), '~/Desktop/Current_work/Brasilian_samples/MODA_transcript/', 'S', rownames(dds@assays@data$counts))
intModules2=WeightedModulePartitionLouvain(t(dds@assays@data$counts), '~/Desktop/Current_work/Brasilian_samples/MODA_transcript/', 'N', rownames(dds@assays@data$counts))
CompareAllNets('~/Desktop/Current_work/Brasilian_samples/MODA_transcript/',intModules1,'S', intModules2, 'N', specificTheta=0.1,conservedTheta=0.1)
# Get the similarity matrix
JaccardMatrix <- comparemodulestwonets('~/Desktop/Current_work/Brasilian_samples/MODA_transcript/',intModules1,intModules2, 'DenseModuleGene_S', 'DenseModuleGene_N')
# Select non conserved modules
nonconserved_modules = which(rowSums(JaccardMatrix)<(1/3))
# Select the genes from non conserved modules
nonconserved_transcripts = unlist(lapply(nonconserved_modules, function(x) read_csv(paste0("../MODA_transcript/DenseModuleGene_N_", x, ".txt"), col_names = F)))
nonconserved_genes=sapply(nonconserved_transcripts, function(x){
  B=strsplit(x, "[[:punct:]]")[[1]]
  return(B[1])})
# GO enrichment analysis on genes in the non conserved modules 
# Get the universe and syntax
aedes_GO=readGAF("~/Desktop/Current_work/Brasilian_samples/VectorBase-57_AaegyptiLVP_AGWG_GO.gaf")
aedesgo_annotation=as.data.frame(cbind(row.names(aedes_GO@setAnnotations), aedes_GO@setAnnotations$term))
colnames(aedesgo_annotation)=c('GOterm','function')
GO_file = read_delim("~/Desktop/Current_work/Brasilian_samples/VectorBase-57_AaegyptiLVP_AGWG_GO.gaf", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
# Do GO enrichment analysis
GOen=enricher(nonconserved_genes, TERM2GENE=cbind(GO_file$X5, GO_file$X2), TERM2NAME=aedesgo_annotation)
# Extract significantly enrich GO terms (3 enrich GO terms without big insights)
res=as.data.frame(GOen@result)
res=res[res$p.adjust<0.05,]
write.table(res, '~/Desktop/Current_work/Brasilian_samples/go_anal_module.csv', row.names = F,sep = "\t",quote = F)
res_latex=res[,c(1:4,6)]
print(xtable(res_latex, type='latex'), include.rownames=FALSE)
revigo=res[,c(1,6)]
write.table(revigo, '~/Desktop/revigo.csv', row.names = F,sep = "\t",quote = F)





