#!/bin/bash

####
# BaseCoverageDistribution.sh
# 5/2/2013
# script to run BaseCoverageDistribution on a bam file
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
if [ $# -le 4 ]
then
	echo -e "\nscript to run  on a BaseCoverageDistribution input file(bam )\
		\nUsage: ./BaseCoverageDistribution.sh </path/to/output folder> <fullpath/to/input BAM> <output file name>\
			<tool info file><memory info file><region bed file(optional)>"
	exit 1;
fi
START=$(date +%s)
output=$1
bam=$2
output_file=$3
tool_info=$4
memory_info=$5

if [ $6 ]
then
	region=$6
	interval="-L $region"
else
	interval=""	
fi	

if [ ! -d $output/temp ]
then
	mkdir -p $output/temp
	temp=$output/temp
	sleep 10s
fi	

mem=$( cat $memory_info | grep -w '^BaseCoverageDistribution_JVM' | cut -d '=' -f2)
check_variable "$memory_info:BaseCoverageDistribution_JVM" $mem
command_line_params=$( cat $tool_info | grep -w '^BaseCoverageDistribution_params' | cut -d '=' -f2 )

java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
check_variable "$tool_info:JAVA" $java
gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
check_variable "$tool_info:GATK" $gatk
ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
check_variable "$tool_info:REF_GENOME" $ref

gatk_param="-R $ref -et NO_ET -K $gatk/Hossain.Asif_mayo.edu.key "

$java/java $mem -Djava.io.tmpdir=$temp -jar $gatk/GenomeAnalysisTK.jar \
-T BaseCoverageDistribution \
--out $output_file \
-I $bam \
$command_line_params $gatk_param $interval
		
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "BaseCoverageDistribution for $bam took $DIFF seconds"
																												