
DESeq2 script
-----------------


Installation requirements:

Install the following R packages:

	1) DESeq2, 
	2) RColorBrewer,
	3) BiocParallel

Sample Gene expression File:

	Input_Matrix.bed

Sample category information File:

	SampleCategoryInfo.txt

	Specifically, 
		first column: sample names, 
		third column: categories.
	
	User should maintain the format of the sample category information (including the column names).


The script
----------

DEseq.R

First check and edit the parameters and file names within the R code, according to the user's input files:
(check lines 25 - 52) 


Then run the command:
Rscript DEseq.R

or run the following script:
DEseq_script.sh


Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org

