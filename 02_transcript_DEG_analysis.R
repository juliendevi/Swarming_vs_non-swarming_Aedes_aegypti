############################################################################
##### RNA-seq whole head swarming and non-swarming males Aedes aegypti #####
#####                TRANSCRIPT LEVEL ANALYSIS                         #####
############################################################################

#### Load all analysis directly ----

load("~/Desktop/data/Brasilian_samples/RNA-seq_analysis/R/kallisto/brasilian_samples2.RData")

### Libraries and working directory ----

# Working directory
setwd("~/Desktop/data/Brasilian_samples/RNA-seq_analysis/R/kallisto/")
# Libraries
library(DESeq2)          # Bioconductor v1.40.2          
library(ggplot2)         # CRAN v3.4.4         
library(ggrepel)         # CRAN v0.9.5         
library(RColorBrewer)    # CRAN v1.1-3    
library(pheatmap)        # CRAN v1.0.12        
library(PoiClaClu)       # CRAN v1.0.2.1       
library(apeglm)          # Bioconductor v1.22.1          
library(tximportData)    # Bioconductor v1.28.0    
library(tximport)        # Bioconductor v1.28.0        
library(factoextra)      # CRAN v1.0.7      
library(tidyverse)       # CRAN v2.0.0       
library(dplyr)           # CRAN v1.1.4           
library(readxl)          # CRAN v1.4.3          
library(readr)           # CRAN v2.1.5            
library(mgsa)            # Bioconductor v1.48.0             
library(clusterProfiler) # Bioconductor v4.8.3  
library(enrichplot)      # Bioconductor v1.20.3       

### Import experiment information and transcript counts from kallisto ----

# import experimental scheme
exp = read_excel("~/Desktop/data/Brasilian_samples/RNA-seq_analysis/R/experiment.xlsx", col_names = T)
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

# Make DESeq2 object
dds <- DESeqDataSetFromTximport(kallisto_out,
                                colData = exp,
                                design = ~ group)
# remove features with low counts
dds = dds[ rowMeans(counts(dds)) > 10, ] 
nrow(dds) # 34477/25059

# extract GO background
bg=rownames(dds)
write(bg, '~/Desktop/background.txt')

# rlog transform counts (by average Tx length and correcting for library size)
rld = rlog(dds, blind=FALSE)



#### Descriptive exploration ----

### PCA
x_pca = prcomp(t(counts(dds)), scale. = T)
# Sum pca
summary(x_pca)
# Exploring axes explanation
fviz_screeplot(x_pca, addlabels = T) + ggtitle('')
fviz_pca_ind(x_pca, asp=1, pointshape=20, habillage = group) + scale_color_discrete()+ggtitle('')

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



#### Check normalisation parameters DESeq2 ----

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
res_significant<- subset(res, log2FoldChange > 0.75 | log2FoldChange < -0.75)
res_significant <- subset(res_significant, padj < 0.05)
nrow(res_significant) #27 total differentially expressed genes
nrow(res_significant[res_significant$log2FoldChange > 0.75,]) #12 upregulated in swarming
nrow(res_significant[res_significant$log2FoldChange < -0.75,]) #15 upregulated in non-swarming

# Save the table
write.table(as.data.frame(res_significant), file="~/Desktop/DEG_FC0.75_p0.05.txt", sep="\t", quote = F, col.names = T, row.names = T)


#### Visualisation differential expression ----

#### Heatmap

# Set number of transcripts to look at
n=27

# Set colours
library(MetBrewer)       # [github::BlakeRMills/MetBrewer] v0.2.0
a=met.brewer('Egypt', n=4, 'discrete')
my_colors = list(
  group = c(Swarming = a[3] , NonSwarming =a[1]))

# Take the counts 
mat = assay(dds)[which(rownames(dds) %in% rownames(res_significant)),]
mat = mat - rowMeans(mat)
df2 = as.data.frame(colData(dds)[,c("group"),drop=FALSE])
colnames(df2)<-c("group")
rownames(df2)=c('S1', 'S2', 'S3', 'S4', 'NS1', 'NS2', 'NS3', 'NS4')
# Reorder rows in the table
mat=mat[match(logorder, rownames(mat)),]

# Get transcript annotation
annotgenes = read_excel('~/Desktop/data/Brasilian_samples/RNA-seq_analysis/R/new_list_genes.xlsx', sheet = 'Feuil1', col_names = F)
for (i in 1:nrow(mat)){
  a=which(annotgenes$...1==rownames(mat)[i])
  rownames(mat)[i]=annotgenes$...6[a]
}
# Change sample names
colnames(mat)=c('S1', 'S2', 'S3', 'S4', 'NS1', 'NS2', 'NS3', 'NS4')


# Generate the heatmap
pheatmap(mat, annotation_col=df2,
         show_rownames = T,annotation_colors = my_colors,
         col=met.brewer("Manet",n=100, type="continuous"),
         legend_breaks = c(-400, 0, 400,800, max(mat)),
         legend_labels = c("-0.2", "0", "0.2","0.4", "z-score"),
         cluster_cols = T,
         cluster_rows = F,
         fontsize = 16, annotation_legend = F, angle_col = 0, annotation_names_col = F, border_color = "white")


#### Volcano plot ----

# This is the code for figure 1B # 

res_df <- as.data.frame(res)
res_df$gene <- row.names(res_df)
res_significant_df <- as.data.frame(res_significant)
res_significant_df$gene <- row.names(res_significant_df)
res_df$sig <- "no"
res_df$sig[res_df$gene %in% res_significant_df$gene] <- "yes"

# New names for graph
deg_names=read_excel("~/Desktop/data/Brasilian_samples/RNA-seq_analysis/R/new_list_genes.xlsx", sheet="Feuil1", col_names = F)
to_change=which(res_df$gene %in% deg_names$...1)

for (i in 1:length(to_change)){
  res_df$gene[to_change[i]]=deg_names$...3[which(res_df$gene[to_change[i]]==deg_names$...1)]
}

# Get 'per' and 'tim' transcript IDs for plot 
per=c("AAEL008141-RA", "AAEL008141-RB", "AAEL008141-RC", "AAEL008141-RD")
tim=c("AAEL019461-RA", "AAEL019461-RB", "AAEL019461-RC", "AAEL019461-RD", "AAEL019461-RE", "AAEL019461-RF", "AAEL019461-RG", "AAEL019461-RH", "AAEL019461-RI", "AAEL019461-RJ", "AAEL019461-RK", "AAEL019461-RL", "AAEL019461-RM")

orto=data.frame(
  stringsAsFactors = FALSE,
  VectorBase.ID = c("AAEL025907-RA","AAEL010711-RD","AAEL008603-RH",
                    "AAEL009895-RI","AAEL009211-RC","AAEL000181-RF",
                    "AAEL027977-RB","AAEL019804-RF",
                    "AAEL018306-RF","AAEL022180-RB","AAEL027189-RA",
                    "AAEL003228-RB","AAEL020244-RA","AAEL025228-RD",
                    "AAEL027400-RG","AAEL020264-RD",
                    "AAEL005950-RS","AAEL002879-RS","AAEL025228-RC",
                    "AAEL021054-RA","AAEL004731-RE","AAEL027790-RF",
                    "AAEL018348-RD","AAEL019711-RC","AAEL006283-RA",
                    "AAEL022806-RC","AAEL019898-RF"),
  LogFC = c(8.81120219008974,8.05491679463046,7.6970423099866,
            6.98034185062882,6.89770142690902,
            6.85649475915555,6.07788099164431,5.26486228341179,
            2.01910788720389,1.28374810701443,
            1.23244554902235,0.964354629261942,-9.90706779481855,
            -9.02913486903167,-8.98939278971224,
            -8.77196159862634,-7.3053249062435,-7.04851442679623,
            -6.39078895206306,-5.54499930568334,
            -1.51420598905933,-1.399664586898,-1.23873378750299,
            -1.17663412042304,-1.14112366907515,
            -1.12954115169289,-0.86508645496451),
  VectorBase.description = c("unspecified product","eph receptor tyrosine kinase",
                             "unspecified product","neprilysin",
                             "unspecified product","polybromo-1",
                             "unspecified product","unspecified product",
                             "unspecified product","unspecified product",
                             "unspecified product",
                             "mitotic protein phosphatase 1 regulator",
                             "Serine/threonine-protein phosphatase 2A activator","unspecified product",
                             "unspecified product","unspecified product",
                             "chloride channel protein 2",
                             "heterogeneous nuclear ribonucleoprotein r","unspecified product",
                             "unspecified product","unspecified product",
                             "unspecified product","unspecified product",
                             "unspecified product","GPCR Myosuppressin",
                             "unspecified product","unspecified product"),
  Fly.gene.name = c("Pyruvate dehydrogenase phosphatase",
                    "eph receptor tyrosine kinase",
                    "FBtr0346461/FBtr0079092/FBtr0342941/FBtr0342940","neprilysin",
                    "Bardet-Biedl syndrome 1","polybromo","lethal (3) neo43",
                    "5-hydroxytryptamine (serotonin) receptor 2A","slowpoke","Cdc42-interacting protein 4",
                    "no orthologs","flyers-cup",
                    "Phosphotyrosyl phosphatase activator","no orthologs",
                    "FBtr0331339/FBtr0331340/FBtr0113327/FBtr0070181",
                    "sidekick","Chloride channel-a","Syncrip",
                    "no orthologs",
                    "Inositol polyphosphate-5-phosphatase type I","Zizimin","no orthologs",
                    "Autophagy-related 17",
                    "Phosphoinositide phospholipase C","Myosuppressin receptor",
                    "FBtr0340116/FBtr0346189/FBtr0074104","no orthologs"),
  Fly.gene.symbol = c("Pdp",
                      "Eph","CG7742","Nep","BBS1","polybromo",
                      "l(3)neo43","5-HT2A","slo","Cip4","n.a.",
                      "f-cup","Ptpa","n.a.","CG3655","sdk","ClC-a",
                      "Syp","n.a.","5Ptasel","Ziz-Zir","n.a.",
                      "Atg17","Plc21C","MsR","CG8578","n.a.")
)
to_change=which(orto$Fly.gene.symbol=='n.a.')
orto$Fly.gene.symbol[to_change]=orto$VectorBase.ID[to_change]
to_change=which(res_df$gene %in% orto$VectorBase.ID)

for (i in 1:length(to_change)){
  res_df$gene[to_change[i]]=orto$Fly.gene.symbol[which(res_df$gene[to_change[i]]==orto$VectorBase.ID)]
}
# Tag 'tim' and 'per' expression data
to_change=which(res_df$gene %in% per)
res_df$gene[to_change]='per'
res_df$sig[to_change] = "per"
to_change=which(res_df$gene %in% tim)
res_df$gene[to_change]='tim'
res_df$sig[to_change] <- "tim"
res_df$order <- ifelse(res_df$sig=="tim", 1, 2)
# Generate the plot
v=ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig), text = paste("Symbol:", gene)) + 
  geom_point(size=1.5)+
  geom_point(data = subset(res_df, sig == 'per'|sig=='tim'),
             aes(size = 2))+
  scale_colour_manual("", breaks=c("no","yes", "tim", "per", 'all'),
                      values = c("black","red", "green", "blue", 'green'),
                      labels=c('non significant', 'significant', 'tim', 'per', 'all'))+
  xlim(-10,10)+
  geom_text_repel(aes(label=ifelse(sig=='no', NA, gene)))+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0.75, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = -0.75, linetype="longdash", colour="grey", size=1) +
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20))#,
#legend.position = "none")
v+guides(colour="legend", size = "none")+theme(legend.position = c(0.1, 0.9), legend.text = element_text(size=12), legend.box.background = element_rect(colour = "black"), legend.title = element_blank() )




