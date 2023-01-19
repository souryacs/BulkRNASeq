Cluster profiler (GSEA)
-------------------------

Gene set enrichment analysis (GSEA) using the tool "Clusterprofiler"

https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html

https://guangchuangyu.github.io/software/clusterProfiler/


Installation requirements:

Install the following R packages:

	1) clusterProfiler, 
	2) enrichplot,
	3) mygene,
	4) ggplot2

The script
----------

Run the following script:
GSEA_ClusterProfiler.sh

Check the parameters.

	*SpeciesType*: 1: human, 2: mouse

	*genelistfile*: sample gene list file. 

	*OutDir*: output directory.

	The final parameter is the p-value significance threshold.

Output
----------

	Within the specified output directory *OutDir*, following files are created:

		*mergeDF.txt* : Input genes and their entrez IDs.

		*dotplot_GSEA_GO.pdf* : Dotplot for GSEA.

		*running_score_GSEA_*.pdf* : Running score for GSEA.

		*dotplot_GSEA_KEGG.pdf* : Dotplot for GSEA using KEGG.








Sample Gene expression File:

	Input_Matrix.bed

Sample category information File:

	SampleCategoryInfo.txt

	Specifically, 
		first column: sample names, 
		third column: categories.
	
	User should maintain the format of the sample category information (including the column names).
