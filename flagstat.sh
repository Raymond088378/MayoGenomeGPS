#!/bin/bash

####
# flagstat.sh
# 4/24/2013
# script to run flagstat on a bam file
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
if [ $# != 5 ]
then
	echo -e "\nscript to run flagstat on a input BAM \
		\nUsage: ./flagstat.sh <sample name> <fullpath/to/input bam file> \
			<tool info file> </path/to/output folder> <which tool(samtools/picard)>"
	exit 1;
fi

START=$(date +%s)
sample=$1
inputbam=$2
output=$3
tool_info=$4
tool=`echo $5 | tr "[A-Z]" "[a-z]"`

##local variables
samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
check_variable "$tool_info:SAMTOOLS" $samtools
picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2)
check_variable "$tool_info:PICARD" $picard
java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
check_variable "$tool_info:JAVA" $java

if [ $tool == "samtools" ]
then
	$samtools/samtools flagstat $inputbam > $output/$sample.flagstat
elif [ $tool == "picard" ]
then
	$java/java -Xmx2g -Xms512m \
    -jar $picard/CollectAlignmentSummaryMetrics.jar \
    INPUT=$inputbam \
    OUTPUT=$output/$sample.flagstat \
    TMP_DIR=$output VALIDATION_STRINGENCY=SILENT
fi

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "flagstat for $sample took $DIFF seconds" 

   