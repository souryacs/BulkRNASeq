library(WGCNA)
library(parallel)

#==========================
# parameters
#==========================

# file containing genes and sample features used for WGCNA 
# preferable in tab delimited bed format
WGCNA_InpFile <- 'WGCNA_Input.bed'

# file containing clinical - trait information
# should be in .csv format
Clinical_Trait_FileName <- "ClinicalTraits_RNASeq.csv"

# output directory for WGCNA
WGCNAOutDir <- paste0(getwd(), '/WGCNA_Out')
system(paste("mkdir -p", WGCNAOutDir))


#==========================
# main code
#==========================

# number of columns (starting from the first column)
# which are non-numbers (i.e. do not represent data but rather gene names etc)
Number_NonData_Columns <- 2

# vector of column indices for the clinical trait data which will be subsequently processed
# here we require first 2 columns for subsequent processing
Clinical_Trait_Data_Columns <- seq(1,2)

# in the second module (network construction)
# this is the minimum module size (recommended by the tutorial)
NETWORK_MIN_MODULE_SIZE <- 30

# in the second module (network construction)
# objective is to use one block for all the input data
# so the maximum size of module equals to the number of lines of the input file
NETWORK_MAX_BLOCK_SIZE <- as.integer(system(paste("cat", WGCNA_InpFile, "| wc -l"), intern = TRUE))

#==========================
# output text file
OutTextFile <- paste0(WGCNAOutDir, '/Out_Summary.txt')
fp_out <- file(OutTextFile, "w")

# this file stores the R session and related variables
# used in the first stage (data input)
Data_Input_stage_1_R_Variables_File <- paste0(WGCNAOutDir, '/WGCNA_01_Input.RData')

# this file stores the R session and related variables
# used in the second stage (automatic network construction)
Auto_Network_Construction_stage_2_R_Variables_File <- paste0(WGCNAOutDir, '/WGCNA_02_network.RData')

# enabling multi-thread
allowWGCNAThreads()

# number of cores in the system
ncore <- detectCores()
cat(sprintf("\n\n\n *** Number of cores in the system: %s \n\n", ncore))

#==========================
# 1. Data input and cleaning
#==========================
if (file.exists(Data_Input_stage_1_R_Variables_File) == FALSE) {

	femData <- read.table(WGCNA_InpFile, header=T, sep="\t", stringsAsFactors=F)
	dim(femData);
	names(femData);

	# leave the non data (number) columns
	# and store the expression (data) in transpose 
	# so that now rows indicate samples, columns indicate transcript ID
	datExpr0 = as.data.frame(t(femData[, -c(1:Number_NonData_Columns)]));

	# sanity check the input data after transpose
	# store the columns corresponding to gene related information (gene ID, gene name)
	nondata_cols <- femData[, 1:Number_NonData_Columns]
	# assign datExpr0
	names(datExpr0) = femData[, 1]
	rownames(datExpr0) = names(femData)[-c(1:Number_NonData_Columns)];

	# filtering of genes - remove low weight value transcripts
	gsg = goodSamplesGenes(datExpr0, verbose = 3);
	outtext <- paste("\n\n **** checking goodSamplesGenes -- gsg$allOK : ", gsg$allOK, "**** \n\n")
	writeLines(outtext, con=fp_out, sep="\n")	

	if (!gsg$allOK) { 
		# Optionally, print the gene and sample names that were removed:
		if ( sum(!gsg$goodGenes) > 0 ) {
			outtext <- paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", "))
			writeLines(outtext, con=fp_out, sep="\n")			
		}
		if (sum(!gsg$goodSamples)>0) {
			outtext <- paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", "))
			writeLines(outtext, con=fp_out, sep="\n")			
		}
		# Remove the offending genes and samples from the data:
		datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes] 
		# also operate the nondata_cols and filter out the bad genes
		nondata_cols <- nondata_cols[gsg$goodGenes, ]
	}

	# clustering function - hierarchical
	sampleTree = hclust(dist(datExpr0), method = "average");

	pdf(file = paste0(WGCNAOutDir, '/sampleClustering1.pdf'), width = 12, height = 9);
	par(cex = 0.6);
	par(mar = c(4,6,2,2))
	plot(sampleTree, main = "Sample clustering to detect outliers", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
	dev.off()

	# hierarchical clustering is applied
	datExpr = datExpr0;
	nGenes = ncol(datExpr)
	nSamples = nrow(datExpr)
	outtext <- paste("\n\n **** after default hierarchical clustering --- nGenes : ", nGenes, "  nSamples : ", nSamples, " **** \n\n")
	writeLines(outtext, con=fp_out, sep="\n")

	# read the disease / control classification + weight matrix
	traitData = read.csv(Clinical_Trait_FileName, check.names=FALSE);
	outtext <- paste("\n\n **** dim(traitData) : ", paste(as.vector(dim(traitData)), sep=" "), " *** \n\n")
	writeLines(outtext, con=fp_out, sep="\n")

	outtext <- paste("\n\n **** names(traitData) : ", paste(as.vector(names(traitData)), sep=" "), " *** \n\n")
	writeLines(outtext, con=fp_out, sep="\n")

	# we include the columns specified in Clinical_Trait_Data_Columns
	# for subsequent processing of traitData
	allTraits = traitData[, Clinical_Trait_Data_Columns];
	# allTraits
	outtext <- paste("\n\n **** dim(allTraits) : ", paste(as.vector(dim(allTraits)), sep=" "), " *** \n\n")
	writeLines(outtext, con=fp_out, sep="\n")

	outtext <- paste("\n\n **** names(allTraits) : ", paste(as.vector(names(allTraits)), sep=" "), " *** \n\n")
	writeLines(outtext, con=fp_out, sep="\n")

	# check sanity of the classification matrix
	# and the complete transcript data
	# row names of datExpr indicate sample names
	# which should also be present in the first column of allTraits
	# corresponding rows in allTraits (with matching sample names) would be processed
	femaleSamples = rownames(datExpr);
	traitRows = match(femaleSamples, allTraits[,1]);
	outtext <- paste("\n\n **** sanity check - number of traitRows with matching sample names : ", length(traitRows), " ** \n\n")
	writeLines(outtext, con=fp_out, sep="\n")

	# 1st column of allTraits are the sample names
	# create a new data frame "datTraits" whose entries are only the weights
	# and whose row names are the sample names (1st column of allTraits)
	# as.data.frame() typecast is essential when the weight data is of only one column
	# also assigning the colnames to datTraits is essential
	datTraits =  as.data.frame(allTraits[traitRows, -1])
	rownames(datTraits) = allTraits[traitRows, 1];
	colnames(datTraits) <- colnames(allTraits)[-1]
	# datTraits

	collectGarbage();

	# Convert traits to a color representation: 
	# white means low, red means high, grey means missing entry
	traitColors = numbers2colors(datTraits, signed = FALSE);

	# Plot the sample dendrogram and the colors underneath.
	# we are using the previous cluster sampleTree 
	# since we did not discard / filter out any sample
	pdf(file = paste0(WGCNAOutDir, '/sampleClustering2.pdf'), width = 16, height = 12);
	par(cex = 0.6);
	par(mar = c(4,6,2,2))
	plotDendroAndColors(sampleTree, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
	dev.off()

	save(datExpr, datTraits, nondata_cols, file = Data_Input_stage_1_R_Variables_File)

}	# end if

#==========================
# 2. Network construction and module detection
# automatic network construction (tutorial section 2.1)
#==========================
if (file.exists(Auto_Network_Construction_stage_2_R_Variables_File) == FALSE) {

# Allow multi-threading within WGCNA.
enableWGCNAThreads()

lnames = load(file = Data_Input_stage_1_R_Variables_File);
# lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
# pick the soft Threshold parameter out of the possible values in "powers" vector
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#*****************
x_soft_power_vec <- sft$fitIndices[,1]
y_topology_model_vec <- (-sign(sft$fitIndices[,3]) * sft$fitIndices[,2])
maxval_y_topology_model_vec <- max(y_topology_model_vec)

# option 1 - 
# picking the soft threshold power manually.Currently not running.
if (0) {
	target_power <- 7
}

# option 2
# picking the soft threshold power automatically
# target power value should be the minimum to achieve 90% of the maximum value
# source: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
if (1) {	
	target_power <- x_soft_power_vec[min(which(y_topology_model_vec >= (0.9 * maxval_y_topology_model_vec)))]
}
#*****************

outtext <- paste("\n\n ***** selected scale free topology power : ", target_power, " **** \n\n")
writeLines(outtext, con=fp_out, sep="\n")

# scale free topology plot - topology vs soft Threshold 
pdf(file = paste0(WGCNAOutDir, '/Soft_Threshold.pdf'), width = 9, height = 6);
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
if (0) {
	# old code - used 0.9 as static threshold
	abline(h=0.90,col="red")
}
if (1) {
	# new code - used 90% of the maximum y axis value
	abline(h=(0.9 * maxval_y_topology_model_vec), col="red")	
}

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# One-step network construction and module detection
# previously default scale free power 7 was used
# also we manually specify the nThreads parameter
# according to the number of cores
net = blockwiseModules(datExpr, checkMissingData = TRUE, power = target_power, TOMType = "unsigned", minModuleSize = NETWORK_MIN_MODULE_SIZE, maxBlockSize = NETWORK_MAX_BLOCK_SIZE, reassignThreshold = 0, mergeCutHeight = 0.10, numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE, saveTOMFileBase = "TOM", nThreads=ncore, verbose = 3)

table(net$colors)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
pdf(file = paste0(WGCNAOutDir, '/cluster_dendrogram_colors.pdf'), width = 15, height = 15)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# save the data
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = Auto_Network_Construction_stage_2_R_Variables_File)
table(moduleColors)

}	# end dummy comment

#=======================
# 3a. Relating modules to external clinical traits and identifying important genes
# note: it uses automatic network construction related data
# Tutorials section 3
#=======================

if (1) {

# Allow multi-threading within WGCNA. 
enableWGCNAThreads()	

# read the disease / control classification + weight matrix
traitData = read.csv(Clinical_Trait_FileName, check.names=FALSE);
dim(traitData)
names(traitData)

# columns in traitData is specified in "Clinical_Trait_Data_Columns" structure
allTraits = traitData[, Clinical_Trait_Data_Columns];
allTraits
dim(allTraits)
names(allTraits)

lnames = load(file = Data_Input_stage_1_R_Variables_File);
# lnames

lnames = load(file = Auto_Network_Construction_stage_2_R_Variables_File);
lnames

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
# datExpr: gene expression data for all the given genes
# moduleColors: per gene assignment of their modules (in terms of color information)
# output structure has a field "eigengenes" which contains the eigen gene per module
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs0
datTraits
# computing correlation per module with different disease categories
MEs = orderMEs(MEs0)

#*****************
# use = "p" arguemnt specifies the pearson correlation
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitCor
# significance of the correlation is tested via student's p-value testing
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix
dim(textMatrix) = dim(moduleTraitCor)

pdf(file = paste0(WGCNAOutDir, '/Module_trait_relationships_Pearson.pdf'), width = 20, height = 20);
par(mar = c(8, 12, 1, 4));
# plot the heatmap of module specific correlation with different disease categories
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = redWhiteGreen(50), textMatrix = textMatrix,setStdMargins = FALSE, cex.text = 1.5, zlim = c(-1,1), main = paste("Module-trait relationships (pearson)"))
dev.off()

# dump the number of genes per module in a .csv file
table(moduleColors)
write.csv(MEs, file=paste0(WGCNAOutDir, '/MEs.csv'), quote=F)

# dump per gene the module assignment
# geneInfo0 = data.frame(gene_id = probes, moduleColor = moduleColors)
probes = names(datExpr)
if (Number_NonData_Columns > 1) {
	# nondata_cols is obtained by loading the Data_Input_stage_1_R_Variables_File
	geneInfo0 <- cbind.data.frame(probes, nondata_cols[, 2:ncol(nondata_cols)], moduleColors)
	colnames(geneInfo0) <- c('gene_id', colnames(nondata_cols)[2:ncol(nondata_cols)], 'moduleColor')
} else {
	geneInfo0 = data.frame(gene_id = probes, moduleColor = moduleColors)	
}
write.csv(geneInfo0, file=paste0(WGCNAOutDir, '/moduleColors.csv'), quote=F)
head(moduleColors)

# dump the peason correlation per module to different disease categories
write.csv(moduleTraitCor, file = paste0(WGCNAOutDir, '/moduleTraits_Pearson_Correlation.csv'), row.names = TRUE)

# dump the significance (p-value) of the above mentioned correlation values
write.csv(moduleTraitPvalue, file = paste0(WGCNAOutDir, '/moduleTrait_Pearson_Pvalue.csv'), row.names = TRUE)

#*****************
# similarly compute the correlation with respect to spearman correlation
# and compute the statistical significance (p-value)
# of the correlation statistics for different modules
moduleTraitCor = cor(MEs, datTraits, use = "p", method="spearman");

# significance of the correlation is tested via student's p-value testing
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

textMatrix =  paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix
dim(textMatrix) = dim(moduleTraitCor)

# sizeGrWindow(20,20)
pdf(file = paste0(WGCNAOutDir, '/Module_trait_relationships_Spearman.pdf'), width = 20, height = 20);
par(mar = c(8, 12, 1, 4));
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(datTraits),yLabels = names(MEs),ySymbols = names(MEs),colorLabels = FALSE,colors = redWhiteGreen(50),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 1.5, zlim = c(-1,1), main = paste("Module-trait relationships (spearman)"))
dev.off()

# dump the spearman correlation per module to different disease categories
write.csv(moduleTraitCor, file = paste0(WGCNAOutDir, '/moduleTraits_Spearman_Correlation.csv'), row.names = TRUE)

# dump the significance (p-value) of the above mentioned correlation values
write.csv(moduleTraitPvalue, file = paste0(WGCNAOutDir, '/moduleTrait_Spearman_Pvalue.csv'), row.names = TRUE)

#************************
# 3.b - Gene relationship to trait and important modules: 
# Gene Significance and Module Membership
# should be simultaneously executed as 3.a
#************************

# names (colors) of the modules
names(allTraits)
modNames = substring(names(MEs), 3)
modNames

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
geneModuleMembership
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# Intramodular analysis: identifying genes with high GS and MM
probes = names(datExpr)
# names(allTraits)

# initialize a data structure which combines all the module membership and gene significance values of individual traits
# this structure will be updated by the below computed
# trait specific correlation statistics
if (Number_NonData_Columns > 1) {
	# nondata_cols is obtained by loading the Data_Input_stage_1_R_Variables_File
	geneInfo0_original <- cbind.data.frame(probes, nondata_cols[, 2:ncol(nondata_cols)], moduleColors)
	colnames(geneInfo0_original) <- c('gene_id', colnames(nondata_cols)[2:ncol(nondata_cols)], 'moduleColor')
} else {
	geneInfo0_original = data.frame(gene_id = probes, moduleColor = moduleColors)	
}

# now check the structure datTraits
# it contains weights for individual traits / categories
# column names denote individual traits
# scroll through individual columns / traits
# extract those columns as a separate data frame
# and perform the correlation and significance (p-value) analysis 
# with respect to the complete set of samples
COLNAMES_datTraits <- colnames(datTraits)

for (colnameIdx in (1:length(COLNAMES_datTraits))) {
	CurrTraitName <- COLNAMES_datTraits[colnameIdx]
	cat(sprintf("\n\n Extracting traits from the structure datTraits -- corresponding to the column : %s \n\n", CurrTraitName))

	# extract the specific column
	# Note: check the syntax
	# if column name is a variable (say v) then the syntax DF$v does not work
	# rather we use DF[[v]]
	currTraitDF <- as.data.frame(datTraits[[CurrTraitName]])
	names(currTraitDF) = CurrTraitName

	# find the correlation between the current trait and the gene expression data
	# store as a data frame
	geneTraitSignificance_CurrTraitName = as.data.frame(cor(datExpr, currTraitDF, use = "p"));
	names(geneTraitSignificance_CurrTraitName) = paste("GS.", names(currTraitDF), sep="");

	# also get the corrlation and p-value with respect to the current trait
	# store as a data frame
	GSPvalue_CurrTraitName = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_CurrTraitName), nSamples))
	names(GSPvalue_CurrTraitName) = paste("p.GS.", names(currTraitDF), sep="");

	# add these two data frame information
	# in the geneInfo0_original data frame
	# again we use [[]] operator to insert one field (column)
	# Note: assigning the last columns of individual data frames
	fieldname <- names(geneTraitSignificance_CurrTraitName)	
	geneInfo0_original[[fieldname]] <- geneTraitSignificance_CurrTraitName[, ncol(geneTraitSignificance_CurrTraitName)]

	fieldname <- names(GSPvalue_CurrTraitName)	
	geneInfo0_original[[fieldname]] <- GSPvalue_CurrTraitName[, ncol(GSPvalue_CurrTraitName)]
}	# end column name processing loop

# write the updated geneInfo0_original, 
# for each gene its modules
# and for each trait, its correlation and p-values
if (1) {
	geneInfo0_original_Dump_FileName <- paste0(WGCNAOutDir, '/geneInfo0_original_Dump.bed')
	write.table(geneInfo0_original, geneInfo0_original_Dump_FileName, row.names=F, col.names=T, sep="\t", quote=F, append=F)
}

# now apply sorting on the module membership data
# based on individual traits of the datTraits structure
# for each trait, first keep a separate copy of the geneInfo0_original structure 
# and then apply sorting
for (colnameIdx in (1:length(COLNAMES_datTraits))) {
	CurrTraitName <- COLNAMES_datTraits[colnameIdx]
	cat(sprintf("\n\n Sorting of module membership for different traits in datTraits -- corresponding to the column : %s \n\n", CurrTraitName))

	# extract the specific column
	# Note: check the syntax
	# if column name is a variable (say v) then the syntax DF$v does not work
	# rather we use DF[[v]]
	currTraitDF <- as.data.frame(datTraits[[CurrTraitName]])
	names(currTraitDF) = CurrTraitName

	# first copy the geneInfo0_original structure
	geneInfo0 = data.frame(geneInfo0_original)

	# now find the correlation between the module eigengenes (for each module)
	# and the current trait (stored in currTraitDF)
	modOrder = order(-abs(cor(MEs, currTraitDF, use = "p")));

	for (mod in 1:ncol(geneModuleMembership)){
		# geneInfo0 structure is dumped in file "geneInfo0_original_Dump_FileName"
		# according to modOrder, module membership of the current trait is 
		# appended at the last two columns		
		oldNames = names(geneInfo0)
		geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], MMPvalue[, modOrder[mod]]);
		names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),paste("p.MM.", modNames[modOrder[mod]], sep=""))
	}

	# current trait significance column name in the geneInfo0 structure
	currColName <- paste("GS.", names(currTraitDF), sep="")

	geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0[[currColName]]))
	geneOrder
	geneInfo = geneInfo0[geneOrder, ]

	# write per module membership information for the current trait
	write.csv(geneInfo, file = paste0(WGCNAOutDir, '/MASTER_GS_MM_trait_', CurrTraitName, '.csv'))

	# also write in a tabular structure
	# write.table(geneInfo, paste0(WGCNAOutDir, '/MASTER_GS_MM_trait_', CurrTraitName, '_auto_network_construction.bed'), row.names=F, col.names=T, sep="\t", append=F, quote=F)

}	# end column name processing loop

# now plot for individual modules (represented by different colors)
# module membership related plots for individual traits

for (colnameIdx in (1:length(COLNAMES_datTraits))) {
	CurrTraitName <- COLNAMES_datTraits[colnameIdx]
	cat(sprintf("\n\n Processing trait : %s for module specific plots \n\n", CurrTraitName))	

	# column name in the structure "geneInfo0_original"
	# corresponding to the current trait
	colname_geneinfo <- paste("GS.", CurrTraitName, sep="")

	# extract the gene trait significance for the current trait
	# using the column name field
	geneTraitSignificance_CurrTraitName <- geneInfo0_original[[colname_geneinfo]]
	
	# loop through all available modules
	# when all modules need to be checked
	for (modIdx in (1:length(modNames))) {

		# current module name
		module <- modNames[modIdx]
		# hard code module name
		# module <- "brown"

		cat(sprintf("\n plotting - checking module : %s ", module))

		# column containing this particular module information
		column = match(module, modNames);
		
		# row index for the genes belonging to the current module
		# moduleGenes = moduleColors==module;
		geneIdx <- which(moduleColors == module)

		# plot the current module membership vs current trait gene significance
		plotfile <- paste0(WGCNAOutDir, '/MM_vs_GS_', CurrTraitName, '_', module, '_module.pdf')
		pdf(file = plotfile, width = 10, height = 8);
		par(mfrow = c(1,1));
		verboseScatterplot(abs(geneModuleMembership[geneIdx, column]), abs(geneTraitSignificance_CurrTraitName[geneIdx]), xlab = paste("Module Membership in", module, "module"), ylab = paste("Gene significance for ", CurrTraitName), main = paste("Module membership vs. gene significance\n"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
		dev.off()

	}	# end module loop

}	# end trait loop


# write the module eigen genes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
write.csv(MEs, file=paste0(WGCNAOutDir, '/MEs.csv'), quote=F)
MET = orderMEs(MEs)

# plot the eigen networks
plotfile <- paste0(WGCNAOutDir, '/EigengeneNetworks.pdf')
pdf(file = plotfile, width = 10, height = 8);
# sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()


}	# end dummy if 0

# close the output text file
close(fp_out)


#===============================================
# step - 4 - visualizing networks 
#===============================================

if (1) {

# Allow multi-threading within WGCNA.
enableWGCNAThreads()	

# start of step 4

lnames = load(file = Data_Input_stage_1_R_Variables_File);

lnames = load(file = Auto_Network_Construction_stage_2_R_Variables_File);

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
# datExpr: gene expression data for all the given genes
# moduleColors: per gene assignment of their modules (in terms of color information)
# output structure has a field "eigengenes" which contains the eigen gene per module
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs0
datTraits
# computing correlation per module with different disease categories
MEs = orderMEs(MEs0)

#**********
# again repeat some steps of stage 2
# to compute TOM and similarity matrix

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
# pick the soft Threshold parameter out of the possible values in "powers" vector
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 1)

x_soft_power_vec <- sft$fitIndices[,1]
y_topology_model_vec <- (-sign(sft$fitIndices[,3]) * sft$fitIndices[,2])
maxval_y_topology_model_vec <- max(y_topology_model_vec)
target_power <- x_soft_power_vec[min(which(y_topology_model_vec >= (0.9 * maxval_y_topology_model_vec)))]

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
TOM <- TOMsimilarityFromExpr(datExpr, power = target_power);
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
plotfile <- paste0(WGCNAOutDir, '/Network_heatmap_all_genes.pdf')
pdf(file = plotfile, width = 9, height = 9);
# sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
plotfile <- paste0(WGCNAOutDir, '/Network_heatmap_selected_genes.pdf')
pdf(file = plotfile, width = 9, height = 9);
# sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

} 	# end dummy if 0



