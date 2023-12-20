##=============
## PCA script
##=============
library('DESeq2')
library('RColorBrewer')
library("BiocParallel")
library('ggplot2')
library('ggrepel')

##===============
## Parameters
##===============
## file containing raw counts
## Assume that the file has header
CountFile <- "Inp_Count.txt"    ## custom count file name - user can edit with the appropriate file path
AllCountData <- read.table(CountFile, header=T, stringsAsFactors=F, sep="\t", check.names=F)

## File containing sample-wise information
## Assume that the file has header
## In this sample file, the first column should contains the sample names
## Note: the sample names (first column) should match with the column names of the "CountFile"
SampleFile <- "Sample_Info.txt" ## custom sample file name - user can edit with the appropriate file path
AllSampleData <- read.table(SampleFile, header=T, row.names=1, stringsAsFactors=F, check.names=F)

## In the sample file, this is the column where condition information (i.e. sample specific variability) is provided
## initialized as "Condition" means that the column "Condition" contains the sample specific variability information
## user can edit this column name according to the input file
## however, then he/she needs to also replace the strings "Condition_1" with the appropriate column name from the code below
condition = "Condition_1"

## output plot file - user can edit this file name
outplotfile <- "Out_PCA_Plot.pdf"

##===============
## DESeq2
##===============

## sanity check
length(which( names(AllCountData) == rownames(AllSampleData) ))

## Select condition-specific data
sampleData <- AllSampleData
countData <- AllCountData[, row.names(sampleData)]
countData <- subset(countData, rowSums(countData) > 0)
length(which( names( countData ) == rownames( sampleData ) ) )

## convert design columns into factors
sampleData$Condition_1 <- as.factor(sampleData$Condition_1)

## DESeq2 design condition
dds <- DESeqDataSetFromMatrix(countData=countData, colData= sampleData, design = ~ Condition_1 )

## perform DESeq2 
dds <- estimateSizeFactors(dds)
counts.Norm <- counts(dds, normalized=TRUE)
counts.Norm <- data.frame(counts.Norm, check.names=F)
#register(MulticoreParam(25))
dds <- DESeq(dds, parallel =T)

##===========
## plot PCA
##===========

se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1), colData=colData(dds))
p <- plotPCA( DESeqTransform( se ), intgroup = condition, ntop = 1000) 

## if sample names to be displayed
# p <- p + geom_text_repel(aes_string(x="PC1", y="PC2", label="name"), color="black")

## increasing the point size
p <- p + geom_point(size=6)

ggsave(outplotfile, width = 12, height = 9)

