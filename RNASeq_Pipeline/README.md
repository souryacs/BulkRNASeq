RNA-seq pipeline
===================

Devloped by : Sourya Bhattacharyya

Supervisors: Dr. Ferhat Ay and Dr. Pandurangan Vijayanand

La Jolla Institute for Immunology

La Jolla, San Diego, CA 92037, USA


Pipeline to analyze RNA-seq data, starting from input Fastq files.


Installation requirements
===========================

	1) HISAT2 (http://daehwankimlab.github.io/hisat2/)

	2) Genome indices from HISAT2 (corresponding to the reference genomes like hg19 or hg38)

	3) GTF file for the reference genome. Can be obtained from the UCSC table browser.

	4) Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic) and the adapters 

	5) fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

	6) samtools (http://www.htslib.org/) with version > 1.6

	7) R (> 3.4 version) and the R package "Rsubread"


Execution
============

	Edit the configuration file "configfile" to indicate the paths of various installed packages and files.

	Run the script "Run_RNASeq_Pipeline_script.sh" 

		- first, edit the parameters in the first few lines.

			- *BaseDir*: Directory containing the fastq files 

			- *ReadName1* and *ReadName2* : .fastq.gz file name formats for paired end reads.

			- *ResultDir* : Output directory.

			- *configfilename* : Configuration file path.


Output
========

	The outputs are stored within the directory *ResultDir*

		- Check the folder *5_HiSAT2*

			- The file *accepted_hits_sorted.bam* contains the sorted alignment file.

		- Check the folder *6_FeatureCounts*

			- The file *FeatureCounts_Expr_FINAL.bed* contains the output feature counts.
	

Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org




