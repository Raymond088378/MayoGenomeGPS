#!/bin/bash

####
# samse.sh (bwa alignment)
# 4/24/2013
# script to run sampe command 
# baheti.saurabh@mayo.edu
####

### function
function check_variable()	{
	message=$1
	if [[ "$2" == "" ]] 
	then 
		echo "$message is not set correctly."
		exit 1
	fi		
}	
if [ $# =! 7 ]
then
	echo -e "\nscript to run samse on a input file \
		\nUsage: ./fastqc.sh <sample name> <sai file for read1><sai file for read2><full path to fastq file read1> <full path to fastq read2><full path to output bam file> \
			<tool info file>"
	exit 1;
fi
START=$(date +%s)

sample=$1
sai_R1=$2
sai_R2=$3
fastq_R1=$4
fastq_R2=$5
outputbam=$6
tool_info=$7


bwa=$( cat $tool_info | grep -w '^BWA' | cut -d '=' -f2)
check_variable "$tool_info:BWA" $bwa
center=$( cat $tool_info | grep -w '^CENTER' | cut -d '=' -f2 )
check_variable "$tool_info:CENTER" $center
platform=$( cat $tool_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
check_variable "$tool_info:PLATFORM" $platform
genome_bwa=$( cat $tool_info | grep -w '^BWA_REF' | cut -d '=' -f2)	
check_variable "$tool_info:BWA_REF" $genome_bwa
samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
check_variable "$tool_info:SAMTOOLS" $samtools
script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
check_variable "$tool_info:WORKFLOW_PATH" $script_path


$bwa/bwa sampe -r "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:$platform\tCN:$center" \
	$genome_bwa $sai_R1 $SAI_r2 $fastq_R1 $fastq_R2 | $samtools/samtools view -bS - > $outputbam

$samtools/samtools view -H $outputbam 1>$outputbam.header 2> $outputbam.log
if [[ `cat $outputbam.log | wc -l` -gt 0 || `cat $outputbam.header | wc -l` -le 0 ]]	
then
	$script_path/errorlog.sh $outputbam sampe.sh ERROR "truncated or corrupt"
	exit 1;
else
	rm $fastq_R1 $fastq_R2
	rm $sai_R1 $sai_R2
fi    	
rm $outputbam.header	
		
				
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "sampe for $sample took $DIFF seconds"		