Sample code to process microarray data
=======================================

Devloped by : Sourya Bhattacharyya

Supervisors: Dr. Ferhat Ay and Dr. Pandurangan Vijayanand

La Jolla Institute for Immunology

La Jolla, San Diego, CA 92037, USA

Code to process the microarray data.



Sample dataset
==============

We use the GEO: GSE84844 data in this sample example
(Homo Sapiens)
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84844


Installation requirements
===========================

As the input data has homo sapiens as the reference species, we need to execute the following commands in a R terminal:

	source("https://bioconductor.org/biocLite.R")
	biocLite("limma")
	biocLite("affycoretools")
	biocLite("oligo")
	biocLite("pd.hg.u133.plus.2")
	biocLite("hgu133plus2.db")

Data
=======

First download the ".CEL.gz" files from the above mentioned GEO repository (using custom wget commands) 
and save them in a specified data directory.

Then create a metadata file (for example, we have shared a file "Affy_targets_Sjogren_GSE84844.txt") containing the sample-wise information. Here, the columns are sampleID, filename, and "Genotype" field (containing disease / healthy information).

Specifically, in the "Genotype" field, we marked two conditions as "Disease" and "Healthy". We request the users to keep these settings.


Execution
=============

Check the script "microarray_limma.r"

Edit the "Parameters" section. 

Currently, sample input files and parameters are provided along with this repository.

Specifically,  

	1) "InputFile" points to the sample-specific metadata file, 
	
	2) "DataDir" points to the folder containing the downloaded ".CEL.gz" files.
	
	3) "LogExpr": sometimes the shared data already contains log transformed normalized gene expression. Then specify this value as 1. Otherwise, 0 (in such a case, log transformation is applied on the input gene expression data, to normalize it).

	4) "OutFile_GeneExpr": Output file (tab delimited text file) to contain the gene expression values. 

	5) "Limma_outfile": Limma generated output file


Then execute the script by typing the command in R terminal:

	Rscript microarray_limma.r


Output
========

	1) The output file (variable "Limma_outfile") contains the complete output.

	2) In addition, the files "Differential_candidates_fdr*.txt" within the output directory contains differential genes subject to different FDR and logFC thresholds.


Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org

