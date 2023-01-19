#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=10GB
#PBS -l walltime=04:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

## Species type
## 1: human, 2: mouse
SpeciesType=1

genelistfile='sample_genelist.txt'
OutDir='out_sample_genelist'
mkdir -p $OutDir

## the final parameter is the p-value significance threshold
Rscript GSEA_ClusterProfiler.r $genelistfile $OutDir $SpeciesType 0.05

