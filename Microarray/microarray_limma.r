##=======================
## sample script to process microarray data
##=======================

library(limma) 
library(affycoretools)
library(oligo)
## Need to install the Platform Design Info for 
## The Manufacturer's Name HG-U133_Plus_2 for this experiment
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)

##========================
## Parameters
##========================
## file containing the metadata (sample-wise) information
InputFile <- "Affy_targets_Sjogren_GSE84844.txt"

## directory containing the downloaded .CEL.gz files
DataDir <- "../Data/"

## boolean indicator if the original gene expression data is already log transformed
LogExpr <- 0

## Output file (tab delimited text file) to contain the affy expression values. 
OutFile_GeneExpr <- "affy_all_Sjogren_GSE84844.txt"

## output file from limma
Limma_outfile <- "limma_complete.txt"

##========================
## main code
##========================

targets <- readTargets(InputFile) 

## Reads CEL files (specified in 'targets') into AffyBatch object.
data <- read.celfiles(filenames=paste0(DataDir, targets$FileName))

## Normalizes data with 'rma' function and assigns them to ExpressionSet object (see above).
eset <- rma(data)

# if the original expression data does not contain fractional values / not normalized
if (LogExpr == 0) {
	exprs(eset) <- log2(exprs(eset)) 
}

eset <- annotateEset(eset, hgu133plus2.db)

## Exports all affy expression values to tab delimited text file. 
write.exprs(eset, file=OutFile_GeneExpr) 

## Creates appropriate design matrix. 
## Alternatively, such a design matrix can be created in any spreadsheet program
## "1" for Disease and "2" for Healthy. 
design <- model.matrix(~ -1+factor(c(rep(1,nrow(targets[targets$Genotype=="Disease",])),rep(2,nrow(targets[targets$Genotype=="Healthy",])))))
colnames(design) <- c("Disease", "Healthy")

## Fits a linear model for each gene based on the given series of arrays.
fit <- lmFit(eset, design)

## Calculate the Disease vs Healthy difference
## contrast matrix
contrast.matrix <- makeContrasts(Disease-Healthy, levels=design) 

## Computes estimated coefficients and standard errors for a given set of contrasts.
fit2 <- contrasts.fit(fit, contrast.matrix) 

# Computes moderated t-statistics and log-odds of differential expression 
## by empirical Bayes shrinkage of the standard errors towards a common value.
fit2 <- eBayes(fit2) 

## Generates list of top 10 ('number=10') differentially expressed genes sorted by B-values ('sort.by=B') for each of the three comparison groups ('coef=1') in this sample set. 
## The summary table contains the following information: 
## 1) logFC: log2-fold change, 
## 2) AveExpr: average expression value accross all arrays and channels
## 3) moderated t-statistic (t): logFC to its standard error, 
## 4) P.Value is the associated p-value,  
## 5) adj.P.Value is the p-value adjusted for multiple testing 
## 6) B-value (B) is the log-odds that a gene is differentially expressed (the-higher-the-better). 
## Usually one wants to base gene selection on the adjusted P-value rather than the t- or B-values. 
## More details on this can be found in the limma PDF manual (type 'limmaUsersGuide()') or on this FAQ page.
topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=10) 

# Exports complete limma statistics table for first comparison group ('coef=1') to tab delimited text file.
write.table(topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=nrow(exprs(eset))), file=Limma_outfile, row.names=T, sep="\t") 

# Filters out candidates that have P-values < 0.05 in each group ('coef=1') 
# and provides the number of candidates for each list.
x <- topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=nrow(exprs(eset))); 
y <- x[x$adj.P.Val < 0.05,];
write.table(y,file="Differential_candidates_fdr0.05.txt",row.names=F,quote=F,sep="\t")

# Same as above but with complex filter: P-value < 0.05 AND at least 2-fold change 
x <- topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=nrow(exprs(eset))); 
y <- x[x$adj.P.Val < 0.05 & (x$logFC > 1 | x$logFC < -1),];
write.table(y,file="Differential_candidates_fdr0.05_and_logFC1.txt",row.names=F,quote=F,sep="\t")

# Same as above but with complex filter: P-value < 0.01
x <- topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=nrow(exprs(eset))); 
y <- x[x$adj.P.Val < 0.01,];
write.table(y,file="Differential_candidates_fdr0.01.txt",row.names=F,quote=F,sep="\t")

# Same as above but with complex filter: P-value < 0.01 AND at least 2-fold change
x <- topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=nrow(exprs(eset))); 
y <- x[x$adj.P.Val < 0.01 & (x$logFC > 1 | x$logFC < -1),];
write.table(y,file="Differential_candidates_fdr0.01_and_logFC1.txt",row.names=F,quote=F,sep="\t")


