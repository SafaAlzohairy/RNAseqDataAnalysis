# get the table of read counts
read.counts <- read.table("/Users/safa/Documents/PhD-Hausbeck lab/Results/summer2017/RNAseq-project/SeqAnalysis/htseqCount_results.txt", header = TRUE)
#the gene IDs should be stored as row.names
row.names(read.counts) <- read.counts$Geneid
#give meaningful sample names
names(read.counts) <- c("C7-1", "C7-2", "C7-3", "C10-1", "C10-2", "C10-3", "C14-1", "C14-2", "C14-3", "C21-1", "C21-2", "C21-3", "D7-1", "D7-2", "D7-3", "D10-1", "D10-2", "D10-3", "D14-1", "D14-2", "D14-3", "D21-1", "D21-2", "D21-3")
str(read.counts)
head(read.counts)
# make a data.frame with meta-data where row.names should match the individual sample names
sample.info <- data.frame(condition = c( rep("C7", 3), rep("C10", 3), rep("C14", 3), rep("C21", 3), rep("D7", 3), rep("D10", 3), rep("D14", 3), rep("D21", 3)), row.names=names(read.counts))
Error in data.frame(condition = c(rep("C7", 3), rep("C10", 3), rep("C14",  : 
  row names supplied are of the wrong length





source ("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
sampleFiles <- list.files(path="/Users/safa/Documents/PhD-Hausbeck lab/Results/summer2017/RNAseq-project/SeqAnalysis/htseq-output",pattern="*.txt")
sampleCondition<- read.table(file="/Users/safa/Documents/PhD-Hausbeck lab/Results/summer2017/RNAseq-project/SeqAnalysis/squash_phenodata.txt",head=TRUE)
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)
directory <- c("/Users/safa/Documents/PhD-Hausbeck lab/Results/summer2017/RNAseq-project/SeqAnalysis/htseq-output/")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition)
##1.Normalization for sequencing depth differences
#get the table of read counts
read.counts <- read.table("/Users/safa/Documents/PhD-Hausbeck lab/Results/summer2017/RNAseq-project/SeqAnalysis/htseqCount_results.txt", header = TRUE)
dim(read.counts) #dimension of the file
str(read.counts)
head(read.counts)
#generate the DESeqDataSet
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition)
# remove genes without any counts
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 0, ]
# investigate different library sizes
colSums(counts(ddsHTSeq)) # should be the same as colSums(read.counts)
#calculate the size factors; normalized read counts can be retrieved via counts(..., normalized = TRUE)
ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
counts.sf_normalized <- counts(ddsHTSeq, normalized = TRUE)
##2.Transformation of sequencing-depth-normalized read counts
##2.1 Log2 transformation of read counts
# transform size-factor normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts <- log2(counts.sf_normalized + 1)
par(mfrow=c(2,1)) # to plot the following two images underneath each other
par(mar = rep(2, 4))
# first, boxplots of non-transformed read counts (one per sample)
boxplot(counts.sf_normalized, notch = TRUE,main = "untransformed read counts", ylab = "read counts")
# box plots of log2-transformed read counts
boxplot(log.norm.counts, notch = TRUE,main = "log2-transformed read counts",ylab = "log2(read counts)")
#To get an impression of how similar read counts are between replicates
plot(log.norm.counts[,1:2], cex=.1, main = "Normalized log2(read counts)")
#Many statistical tests and analyses assume that data is homoskedastic, i.e. that all variables have similar variance. However, data with large differences among the sizes of the individual observations often shows heteroskedastic behavior. One way to visually check for heteroskedasticity is to plot the mean vs. the standard deviation
#install the vsn package
source("http://bioconductor.org/biocLite.R")
biocLite("vsn")
library(vsn)
meanSdPlot(log.norm.counts, ranks = FALSE, ylim = c(0,3), main="sequencing depth normalized log2(read counts)")
#Transformation of read counts including variance shrinkage
#To reduce the amount of heteroskedasticity, DESeq2 offers means to shrink the variance of low read counts by using the dispersion-mean trend seen for the entire data set as a reference. That means that genes with low, but highly variable read counts across the different samples and replicates will be assigned more homogeneous read counts so that their variance resembles the variance of genes with higher read counts. DESeq2’s rlog() function returns values that are both normalized for sequencing depth and transformed to the log2 scale where the values are adjusted to fit the experiment-wide trend of the variance-mean relationship.
# obtain regularized log-transformed values
rlog.DESeq.sumExp <- rlog(ddsHTSeq, blind = TRUE)
rlog.norm.counts <- assay(rlog.DESeq.sumExp)
# mean-sd plot for rlog-transformed data
meanSdPlot(rlog.norm.counts, ranks = FALSE, ylim = c(0,3), main = "rlog-transformed read counts")
#Read count correlations, the expression values of individual genes should be fairly similar for biological replicates
#Hierarchical clustering To determine whether the different sample types can be separated in an unsu- pervised fashion
# cor() calculates the correlation between columns of a matrix
distance.m_rlog <- as.dist(1 - cor(rlog.norm.counts, method = "pearson" ))
# plot() can directly interpret the output of hclust()
plot( hclust(distance.m_rlog), labels = colnames(rlog.norm.counts), main = "rlog transformed read counts\ndistance: Pearson correlation")
#Principle Component analysis PCA, to determine whether samples display greater variability between experimental conditions than between replicates of the same treatment
library(ggplot2)
P <- plotPCA(rlog.DESeq.sumExp)
# plot cosmetics
P <- P + theme_bw() + ggtitle("Rlog transformed counts")
p1<- P + theme_bw() + ggtitle("seq.depth normalized")
print(P)
print(p1)
