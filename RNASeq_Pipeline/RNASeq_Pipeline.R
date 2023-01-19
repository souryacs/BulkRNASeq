library(Rsubread)

countToTpm <- function(counts, effLen)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    val <- exp(rate - denom + log(1e6))
    return(val)
}
 
countToFpkm <- function(counts, effLen)
{
    N <- sum(counts)
    val <- exp( log(counts) + log(1e9) - log(effLen) - log(N) )
    return(val)
}
 
# currently not used
fpkmToTpm <- function(fpkm)
{
    val <- exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
    return(val)
}

args <- commandArgs(TRUE)

alignfile <- args[1]
outdir <- args[2]
GTFAnnotFile <- args[3]
PEinfo <- as.integer(args[4])
GeneLenFile <- args[5]

if (PEinfo == 0) {
	FeatCntOut <- Rsubread::featureCounts(alignfile, annot.ext=GTFAnnotFile, isGTFAnnotationFile=TRUE)
} else {
	FeatCntOut <- Rsubread::featureCounts(alignfile, annot.ext=GTFAnnotFile, isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, requireBothEndsMapped=TRUE)
}

OutAnnot <- FeatCntOut$annotation
CntDF <- cbind.data.frame(OutAnnot[, 1], OutAnnot[, ncol(OutAnnot)], FeatCntOut$counts)
colnames(CntDF) <- c('gene_id', 'gene_len_annot', 'count')

GeneLenData <- read.table(GeneLenFile, header=F, sep="\t", stringsAsFactors=F)
colnames(GeneLenData) <- c('gene_id', 'gene_name', 'gene_len_GTF')

OutCntDF <- merge(GeneLenData, CntDF)
write.table(OutCntDF, paste0(outdir, '/FeatureCounts_FINAL.bed'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

write.table(FeatCntOut$counts_junction, paste0(outdir, '/Counts_Junction.bed'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
write.table(FeatCntOut$annotation, paste0(outdir, '/annotation.bed'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
write.table(FeatCntOut$stat, paste0(outdir, '/stat.bed'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

tpm_val_GTF <- countToTpm(OutCntDF$count, OutCntDF$gene_len_GTF)
fpkm_val_GTF <- countToFpkm(OutCntDF$count, OutCntDF$gene_len_GTF)

tpm_val_annot <- countToTpm(OutCntDF$count, OutCntDF$gene_len_annot)
fpkm_val_annot <- countToFpkm(OutCntDF$count, OutCntDF$gene_len_annot)

OutCntDF_Expr <- cbind.data.frame(OutCntDF, tpm_val_GTF, fpkm_val_GTF, tpm_val_annot, fpkm_val_annot)
colnames(OutCntDF_Expr) <- c(colnames(OutCntDF), 'TPM_GTF', 'FPKM_GTF', 'TPM_Annot_Length', 'FPKM_Annot_Length')
write.table(OutCntDF_Expr, paste0(outdir, '/FeatureCounts_Expr_FINAL.bed'), row.names=F, col.names=T, sep="\t", quote=F, append=F)




