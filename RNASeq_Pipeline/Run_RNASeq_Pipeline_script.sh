#!/bin/bash

##======================
## parameters
##======================

# directory containing the fastq files 
BaseDir='/home/sourya/Data'

# executable of RNA Seq SLE piepeline code
codeexec='RNASeq_Pipeline.sh'

# file name - read 1
ReadName1="_1.fastq.gz"

# file name - read 2
ReadName2="_2.fastq.gz"

# directory to store the results
ResultDir='/home/sourya/Results'

## configuration file storing the parameters
configfilename=`pwd`'/configfile'

##======================
## loop
##======================
 
# script counter
c=0

# search for the fastq files 
for ff in `find $BaseDir -type f -name "*"$ReadName1`; do

	echo 'processing input file: '$ff
	
	inpdir=$(dirname "$ff")

	fastq1=`find $inpdir -type f -name "*"$ReadName1`
	fastq2=`find $inpdir -type f -name "*"$ReadName2`

	echo "fastq1: "$fastq1
	echo "fastq2: "$fastq2

	sample=$(basename "$inpdir")

	outdir=$ResultDir'/'$sample
	mkdir -p $outdir

	# output directory
	echo 'output dir: '$outdir

	# use this fastq file as an input to the RNA seq pipeline
	c=`expr $c + 1`
	scriptfile='temp_script_RNASeq_fastq_'${c}'.sh'

echo '' > ${scriptfile}
cat <<EOT >> ${scriptfile}
#!/bin/bash -ex
#PBS -m ae
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -l walltime=48:00:00
#PBS -j eo

${codeexec} -f ${fastq1} -r ${fastq2} -d ${outdir} -C ${configfilename}

EOT

     chmod +x ${scriptfile}
     qsub ${scriptfile}

done

