
Various scripts for processing bulk RNA-seq data
-----------------------------------------------------

Devloped by : Sourya Bhattacharyya

Supervisors: Dr. Ferhat Ay and Dr. Pandurangan Vijayanand

La Jolla Institute for Immunology

La Jolla, San Diego, CA 92037, USA


Scripts:
==========

Check the following folders for respective scripts and README:

1. RNASeq_Pipeline

    - Pipeline to process bulk RNA-seq fastq files, using HISAT2 as the aligner. 
    - Output: gene expression
    - Details are provided in the respective documentation

2. DESeq2

    - Script to run DESeq2 to obtain differentially expressed genes.
    - Details are provided in the respective documentation

3. WGCNA:

    - Script to run WGCNA to identify networks of interconnected gene regulatory modules.
    - Details are provided in the respective documentation

4. PathwayAnalysis

    - Sample script to perform the pathway anaysis on a set of genes using the package *clusterprofiler* (https://guangchuangyu.github.io/software/clusterProfiler/)

5. PCA

    - Sample script to perform PCA on RNA-seq counts (and corresponding samples) according to the given groups of input samples.
    - User needs to check the Parameters section of the R script, edit the parameters if needed, and run PCA.
    - This PCA is based on DESeq2 routines, so user needs to install DESeq2 and other dependencies mentioned in the R script.


Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org
