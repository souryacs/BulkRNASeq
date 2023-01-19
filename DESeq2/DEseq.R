library('DESeq2')
library('RColorBrewer')
library("BiocParallel")

#========================
# DESEQ constants
# FDR threshold
FDR_Th_DESeq <- 0.05
# log2 fold change threshold
FOLD_Change_Thr <- 1
FOLD_Change_Thr2 <- 2
#========================

#========================
# FILENAME and PARAMETER SETTINGS - TO EDIT  !!!
#========================

# input gene expression file
# first two fields are geneID and geneName
## we assume that column names are provided
InpMatrixFile <- 'Input_Matrix.bed'

# this variable indicates the number of columns which contain gene related information
# here it is 2 since the first two fields of the input count data are geneID and geneName
HEADERFIELDS <- 2

# this file contains categorical information for all samples (contains header)
# 1st column: sample name
# when reading this file, set row.names=1
# so that the sample names become the row names of input data
InpCategoryFile <- 'SampleCategoryInfo.txt'

# Conditions 
## the "condition" value denotes the column name in the "InpCategoryFile" 
## storing the category information
condition = "Condition"

## two input categories
## check the input category file to put the corresponding categories
conditionA = "Category1"
conditionB = "Category2"

comparison_name = paste0(conditionA, "_vs_", conditionB)

# output directory corresponding to the current DeSeq comparison
DESEQOutDir <- paste0('DESeq2_Output_', comparison_name)
system(paste("mkdir -p", DESEQOutDir))

#========================
# END FILENAME and PARAMETER SETTINGS
#========================

# read input count information
AllCountData <- read.table(InpMatrixFile, header=T, sep="\t", stringsAsFactors=F, check.names=F)
Colname_Mat <- colnames(AllCountData)
cat(sprintf("\n\n **** AllCountData dimension --- nrow : %s ncol : %s ", nrow(AllCountData), ncol(AllCountData)))

# read input category information
AllSampleData <- read.table(InpCategoryFile, header=T, sep="\t", row.names=1, stringsAsFactors=F, check.names=F)

# Select condition-specific data ## Edit
sampleData <- AllSampleData[( AllSampleData$Condition == conditionA ) | ( AllSampleData$Condition == conditionB ), ]
countData <- AllCountData[,row.names(sampleData)]

# select only those instances where not all count values are zero
Selected_Idx <- which(rowSums(countData) > 0)
countData <- countData[Selected_Idx, ]
cat(sprintf("\n\n ===>> reduced number of samples in countData : %s ==== \n", nrow(countData)))

# also truncate the gene information accordingly
Selected_GeneInfo <- AllCountData[Selected_Idx, 1:HEADERFIELDS]

# sanity check
length(which( names( countData ) == rownames( sampleData ) ) )

# convert design columns into factors
sampleData$Condition <- as.factor(sampleData$Condition)

dds <- DESeqDataSetFromMatrix(countData=countData, colData= sampleData, design = ~ Condition )

dds <- estimateSizeFactors(dds)
counts.Norm <- counts(dds, normalized=TRUE)
counts.Norm <- data.frame(counts.Norm)
write.csv(counts.Norm, file=paste0(DESEQOutDir, "/countDataDEseqNorm.csv"))

#register(MulticoreParam(25))
dds <- DESeq(dds, parallel =T)

res <- results(dds, parallel =T, contrast= c( condition, conditionA, conditionB), cooksCutoff=FALSE, addMLE=TRUE )

# append DESeq results with the original count data structure
# and write it in a separate file
DESEQData <- cbind.data.frame(Selected_GeneInfo, countData, res)
colnames(DESEQData) <- c(colnames(Selected_GeneInfo), colnames(countData), colnames(res))
write.table(DESEQData, paste0(DESEQOutDir, '/DESeq_Results_ALL.bed'), row.names = F, col.names = T, sep = "\t", quote=F, append=F)

# obtain significant results in terms of FDR <= FDR_Th_DESeq 
SigIdx <- which(DESEQData$padj <= FDR_Th_DESeq)
write.table(DESEQData[SigIdx, ], paste0(DESEQOutDir, '/DESeq_Results_SIG_FDR_', FDR_Th_DESeq, '.bed'), row.names = F, col.names = T, sep = "\t", quote=F, append=F)

# obtain significant results in terms of FDR <= FDR_Th_DESeq and log2 fold change >= FOLD_Change_Thr
SigIdx <- which((DESEQData$padj <= FDR_Th_DESeq) & (abs(DESEQData$log2FoldChange) >= FOLD_Change_Thr))
write.table(DESEQData[SigIdx, ], paste0(DESEQOutDir, '/DESeq_Results_SIG_FDR_', FDR_Th_DESeq, '_LogFC_', FOLD_Change_Thr, '.bed'), row.names = F, col.names = T, sep = "\t", quote=F, append=F)

# obtain significant results in terms of FDR <= FDR_Th_DESeq and log2 fold change >= FOLD_Change_Thr2
SigIdx <- which((DESEQData$padj <= FDR_Th_DESeq) & (abs(DESEQData$log2FoldChange) >= FOLD_Change_Thr2))
write.table(DESEQData[SigIdx, ], paste0(DESEQOutDir, '/DESeq_Results_SIG_FDR_', FDR_Th_DESeq, '_LogFC_', FOLD_Change_Thr2, '.bed'), row.names = F, col.names = T, sep = "\t", quote=F, append=F)

# plot MA
pdf(file = paste0(DESEQOutDir, '/plotMA.pdf'), width = 12, height = 9);
plotMA(res, main=paste0(comparison_name), ylim=c(-5,5))
dev.off()

# plot PCA
pdf(file = paste0(DESEQOutDir, '/plotPCA.pdf'), width = 12, height = 9);
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1), colData=colData(dds))
plotPCA( DESeqTransform( se ), intgroup = 'Condition', ntop = 5000 )
dev.off()

# plot Dispersions
dds <- estimateDispersions(dds)
pdf(file = paste0(DESEQOutDir, '/plotDispersions.pdf'), width = 12, height = 9);
plotDispEsts(dds, genecol = "black", fitcol = "red", finalcol = "dodgerblue", legend = TRUE, log = "xy", cex = 0.45)
dev.off()

# plot Sparsity
pdf(file = paste0(DESEQOutDir, '/plotSparsity.pdf'), width = 12, height = 9);
plotSparsity(dds)
dev.off()

# Plot counts for a single gene
pdf(file = paste0(DESEQOutDir, '/singleGene.pdf'), width = 12, height = 9);
plotCounts(dds, gene=which.min(res$padj), intgroup = 'Condition', pch = 19)
dev.off()

# sorting the DESeq results
resOrd <- res[order(res$padj),]

# plot Heatmap
resOrdered <- resOrd[!is.na(resOrd$padj) & resOrd$padj < FDR_Th_DESeq, ]
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0.585, ], resOrdered[ resOrdered[,'log2FoldChange'] < -0.585, ] )
hmcol <- brewer.pal(11,'RdBu')
nCounts <- counts(dds, normalized=TRUE)
pdf(file = paste0(DESEQOutDir, '/topHeatmap.pdf'), width = 12, height = 9);
heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(10,2))
dev.off()


