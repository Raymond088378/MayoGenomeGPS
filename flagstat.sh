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
if [ $# != 4 ]
then
	echo -e "\nscript to run flagstat on a input BAM \
		\nUsage: ./flagstat.sh <fullpath/to/input bam file> </full path/outputfile name>\
			<tool info file> <which tool(samtools/picard)>"
	exit 1;
fi

START=$(date +%s)
inputbam=$1
outputfile=$2
tool_info=$3
tool=`echo $4 | tr "[A-Z]" "[a-z]"`

##local variables
samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
check_variable "$tool_info:SAMTOOLS" $samtools
picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2)
check_variable "$tool_info:PICARD" $picard
java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
check_variable "$tool_info:JAVA" $java

if [ $tool == "samtools" ]
then
	$samtools/samtools flagstat $inputbam > $outputfile
elif [ $tool == "picard" ]
then
	$java/java -Xmx2g -Xms512m \
    -jar $picard/CollectAlignmentSummaryMetrics.jar \
    INPUT=$inputbam \
    OUTPUT=$outputfile \
    TMP_DIR=$output VALIDATION_STRINGENCY=SILENT
fi

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "flagstat for $sample took $DIFF seconds" 

   