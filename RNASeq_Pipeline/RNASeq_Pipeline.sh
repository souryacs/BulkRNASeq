#!/bin/bash

# usage info
usage(){
cat << EOF

Options:    

	-C  ConfigFile		 Name of the configuration file for executing RNA-seq pipeline
	-f  FASTQ1           R1 of pair-end sequencing data  [.fq|.gz|.bz2]. 
	-r  FASTQ2           R2 of pair-end sequencing data [.fq|.gz|.bz2]. 
						 If not provided, the input is assumed to be single ended.	
	-d  OutDir 			 Set the output directory which will contain all the results

EOF
}

while getopts "C:f:r:g:d:" opt;
do
	case "$opt" in
		C) ConfigFile=$OPTARG;;
		f) FASTQ1=$OPTARG;;
		r) FASTQ2=$OPTARG;;
		d) OutDir=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $ConfigFile ]]; then
	echo 'Configuration file is not provided - check the option -C - quit !! '
	exit 1
fi

if [[ ${OutDir: -1} == "/" ]]; then
	OutDir=${OutDir%?}
fi
mkdir -p $OutDir
echo '**** Output directory: '$OutDir

#----------------------------------
current_dir=$(pwd)
script_dir=$(dirname $0)
cd $script_dir

#=================================
# parse the configuration file
#=================================
echo -e "\n ================ Parsing input configuration file ================="

# separator used in the config file
IFS="="
while read -r name value
do
	param=$name
	paramval=${value//\"/}
	if [[ -n $param ]]; then
		if [[ $param != \#* ]]; then
			echo -e "Content of $param is $paramval"
			if [ $param == "GENOME" ]; then
				HISAT_GENOME=$paramval
			fi
			if [ $param == "GTFFile" ]; then
				RefGTFFile=$paramval
			fi
			if [ $param == "TrimmomaticExec" ]; then
				Trimmomatic_exec=$paramval
			fi
			if [ $param == "TrimmomaticAdapter" ]; then
				Trimmomatic_adapter_path=$paramval
			fi		
		fi
	fi
done < $ConfigFile


echo "*** HISAT_GENOME: "$HISAT_GENOME
echo "*** RefGTFFile: "$RefGTFFile
#=========================

#=============================
# reference packages
#=============================
fastqc_exec=`which fastqc`

# fastq_illumina_filter_exec=`which fastq_illumina_filter`

# HiSAT2 package usage
HiSAT2Exec=`which hisat2`

# samtools package usage 
samtools_Exec=`which samtools`

#=============================
# determining single end or paired end input
#=============================
if [[ -z "$FASTQ2" ]]; then
	echo "Single end read fastq file is provided as input"
	PE=0
else
	echo "Paired end read fastq file is provided as input"
	PE=1
fi

#=============================
# applying fastQC 1 for quality check
#=============================
fastqc_1_out_dir=$OutDir'/1_fastQC_1'
mkdir -p $fastqc_1_out_dir

outzipfilecount=`find $fastqc_1_out_dir -type f -name "*.zip" | wc -l`
outhtmlfilecount=`find $fastqc_1_out_dir -type f -name "*.html" | wc -l`

if [[ $outzipfilecount == 0 || $outhtmlfilecount == 0 ]]; then
	if [[ $PE == 0 ]]; then
		$fastqc_exec -o $fastqc_1_out_dir -f fastq $FASTQ1
	else
		$fastqc_exec -o $fastqc_1_out_dir -f fastq $FASTQ1 $FASTQ2
	fi
fi

#=============================
# applying Trimmomatic - adapter trimming
#=============================
Trim_out_dir=$OutDir'/3_Trimmomatic'
mkdir -p $Trim_out_dir

if [[ $PE == 0 ]]; then

	# The following command will execute for Single End Read.
	Trim_out_file=$Trim_out_dir'/Trimmomatic_out.fastq.gz'
	if [[ ! -f $Trim_out_file ]]; then	
		java -jar $Trimmomatic_exec SE -phred33 -trimlog $Trim_out_dir'/Trim_logfile.log' $FASTQ1 $Trim_out_file ILLUMINACLIP:$Trimmomatic_adapter_path'/TruSeq3-SE.fa':2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75
	fi
else 

	# The following command will execute for paired End Read.
	Trim_R1=$Trim_out_dir'/Trimmomatic_out_R1.fastq.gz'
	Trim_R1_Unpair=$Trim_out_dir'/Trimmomatic_out_R1_unpaired.fastq.gz'
	Trim_R2=$Trim_out_dir'/Trimmomatic_out_R2.fastq.gz'
	Trim_R2_Unpair=$Trim_out_dir'/Trimmomatic_out_R2_unpaired.fastq.gz'
	if [[ ! -f $Trim_R1 || ! -f $Trim_R2 ]]; then	
		java -jar $Trimmomatic_exec PE -phred33 -trimlog $Trim_out_dir'/Trim_logfile.log' $FASTQ1 $FASTQ2 $Trim_R1 $Trim_R1_Unpair $Trim_R2 $Trim_R2_Unpair ILLUMINACLIP:$Trimmomatic_adapter_path'/TruSeq3-PE-2.fa':2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75
	fi
fi

#=============================
# applying fastQC-2 for quality check
#=============================
fastqc_out_dir=$OutDir'/4_fastQC_2'
mkdir -p $fastqc_out_dir

if [[ $PE == 0 ]]; then
	if [[ ! -f $fastqc_out_dir'/Trimmomatic_out_fastqc.zip' ]]; then
		$fastqc_exec -o $fastqc_out_dir -f fastq $Trim_out_file
	fi
else
	if [[ ! -f $fastqc_out_dir'/Trimmomatic_out_R1_fastqc.zip' || ! -f $fastqc_out_dir'/Trimmomatic_out_R2_fastqc.zip' ]]; then
		$fastqc_exec -o $fastqc_out_dir -f fastq $Trim_R1 $Trim_R2
	fi
fi

#=============================
# applying HISAT2
#=============================
HiSAT2_out_dir=$OutDir'/5_HiSAT2'
mkdir -p $HiSAT2_out_dir

HISAT2_OutSAMFile=$HiSAT2_out_dir'/accepted_hits.sam'
HISAT2_SummaryFile=$HiSAT2_out_dir'/alignment_summary.txt'

if [[ $PE == 0 ]]; then
	if [[ ! -f $HISAT2_OutSAMFile ]]; then
		$HiSAT2Exec --threads 8 --no-mixed --no-discordant -x ${HISAT_GENOME} -q -U $Trim_out_file -S $HISAT2_OutSAMFile --summary-file $HISAT2_SummaryFile --add-chrname
	fi
else
	if [[ ! -f $HISAT2_OutSAMFile ]]; then
		$HiSAT2Exec --threads 8 --no-mixed --no-discordant -x ${HISAT_GENOME} -q -1 $Trim_R1 -2 $Trim_R2 -S $HISAT2_OutSAMFile	--summary-file $HISAT2_SummaryFile --add-chrname
	fi
fi

# sort the HISAT2 generated output SAM file into BAM format
if [ ! -f $HiSAT2_out_dir'/accepted_hits.bam' ]; then
	samtools view -bhS -o $HiSAT2_out_dir'/accepted_hits.bam' $HISAT2_OutSAMFile
fi

# sort the HISAT2 package generated alignment file
if [[ ! -f $HiSAT2_out_dir'/accepted_hits_sorted.bam' ]]; then
	# works for samtools version >= 1.6
	$samtools_Exec sort -o $HiSAT2_out_dir'/accepted_hits_sorted.bam' $HiSAT2_out_dir'/accepted_hits.bam'
	# works for older samtools version
	# samtools sort $HiSAT2_out_dir/accepted_hits.bam $HiSAT2_out_dir/accepted_hits_sorted
fi

#===================
# FeatureCounts package use
#===================
FeatureCountOutDir=$OutDir'/6_FeatureCounts'
mkdir -p $FeatureCountOutDir

GeneLengthFile=$FeatureCountOutDir'/Gene_Length.bed'
awk -F[' \t'] '{if ((substr($1,1,1) != "#") && ($3=="gene")) {printf "%d",($5-$4); for (i=9; i<=NF; i++) {if ($i=="gene_name") {printf "\t%s\n",$(i+1)}; if ($i=="gene_id") {printf "\t%s",$(i+1)}}}}' $RefGTFFile | awk '!seen[$2]++' - | awk '{print $2"\t"$3"\t"$1}' - > $GeneLengthFile
sed -i 's/"//g' $GeneLengthFile
sed -i 's/;//g' $GeneLengthFile

Rscript RNASeq_Pipeline_Using_Hisat2_featureCount.R $HiSAT2_out_dir'/accepted_hits_sorted.bam' $FeatureCountOutDir $RefGTFFile $PE $GeneLengthFile

cd $current_dir

