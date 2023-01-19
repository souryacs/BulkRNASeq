library(clusterProfiler)
library(enrichplot)
library(mygene)
library(ggplot2)

args <- commandArgs(TRUE)
GeneListFile <- args[1]
OutDir <- args[2]
SpeciesType <- as.integer(args[3])
if (length(args) > 3) {
	pval_thr <- as.numeric(args[4])
} else {
	pval_thr <- 0.05
}
cat(sprintf("\n *** pval_thr : %s ", pval_thr))

system(paste("mkdir -p", OutDir))

genelistDF <- read.table(GeneListFile, header=F, sep="\t", stringsAsFactors=F)
colnames(genelistDF) <- c('GeneName', 'FC')
cat(sprintf("\n dimension genelistDF : %s X %s  ", nrow(genelistDF), ncol(genelistDF)))

## convert to the Entrez ID 
if (SpeciesType == 1) {
	p <- queryMany(genelistDF[,1], scopes="symbol", fields="entrezgene", species="human")
} else {
	p <- queryMany(genelistDF[,1], scopes="symbol", fields="entrezgene", species=10090)	#"mouse")
}
outdf <- p[, c("query", "entrezgene")]
colnames(outdf) <- c('GeneName', 'EntrezID')
outdf <- outdf[which(outdf[,2] != "NA"), ]
outdf <- outdf[!duplicated(outdf$GeneName), ]
cat(sprintf("\n dimension outdf : %s X %s  ", nrow(outdf), ncol(outdf)))

mergeDF <- merge(genelistDF, outdf)
cat(sprintf("\n dimension mergeDF : %s X %s  ", nrow(mergeDF), ncol(mergeDF)))

## remove any entries with NA in 'entrezID'
mergeDF <- mergeDF[which(mergeDF[,3] != "NA"), ]

write.table(mergeDF, paste0(OutDir, '/mergeDF.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

genelistDF <- unique(mergeDF[, c('EntrezID', 'FC')])
cat(sprintf("\n dimension genelistDF (after unique) : %s X %s  ", nrow(genelistDF), ncol(genelistDF)))

## sort the genes according to decreasing order of the fold change
geneList = as.vector(genelistDF[, 2])
names(geneList) = as.character(genelistDF[,1])
geneList = sort(geneList, decreasing = TRUE)

## GO Gene set enrichment analysis
cat(sprintf("\n\n ==>> start GO GSEA analysis "))
if (SpeciesType == 1) {
	ego2 <- gseGO(gene=geneList, OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=pval_thr, pAdjustMethod="none", verbose=FALSE)
} else {
	ego2 <- gseGO(gene=geneList, OrgDb="org.Mm.eg.db", ont="BP", pvalueCutoff=pval_thr, pAdjustMethod="none", verbose=FALSE)
}
write.table(ego2, paste0(OutDir, '/ego2.txt'), row.names=T, col.names=T, sep="\t", quote=F, append=F)

## dotplot
plotfile2 <- paste0(OutDir, '/dotplot_GSEA_GO.pdf')
dotplot(ego2, showCategory=30) + ggtitle("dotplot for GO GSEA")
ggsave(plotfile2, width=8, height=10)

## running score and preranked list of GSEA result
plotfile3A <- paste0(OutDir, '/running_score_GSEA_GO_1.pdf')
gseaplot(ego2, geneSetID = 1, by = "runningScore", title = ego2$Description[1])
ggsave(plotfile3A, width=8, height=6)

plotfile3B <- paste0(OutDir, '/running_score_GSEA_GO_2.pdf')
gseaplot(ego2, geneSetID = 1, by = "preranked", title = ego2$Description[1])
ggsave(plotfile3B, width=8, height=6)

plotfile3C <- paste0(OutDir, '/running_score_GSEA_GO_3.pdf')
gseaplot(ego2, geneSetID = 1, title = ego2$Description[1])
ggsave(plotfile3C, width=8, height=6)

cat(sprintf("\n\n ==>> end GO GSEA analysis "))

## KEGG pathway enrichment analysis
cat(sprintf("\n\n ==>> start KEGG pathway enrichment analysis "))
if (SpeciesType == 1) {
	kk2 <- gseKEGG(geneList=geneList, organism="hsa", pvalueCutoff=pval_thr, pAdjustMethod = "none", verbose=FALSE)
} else {
	kk2 <- gseKEGG(geneList=geneList, organism="mmu", pvalueCutoff=pval_thr, pAdjustMethod = "none", verbose=FALSE)
}
plotfile3 <- paste0(OutDir, '/dotplot_GSEA_KEGG.pdf')
dotplot(kk2, showCategory=30) + ggtitle("dotplot for KEGG GSEA")
ggsave(plotfile3, width=8, height=10)
cat(sprintf("\n\n ==>> end KEGG pathway enrichment analysis "))

