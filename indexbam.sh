#!/bin/bash

####
# indexbam.sh
# 4/24/2013
# script to index a bam file
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
if [ $# -le 2 ]
then
	echo -e "\nscript to generate bam index \
		\nUsage: ./indexbam.sh <fullpath/to/input file bam> \
			<tool info file>"
	exit 1;
fi
START=$(date +%s)
bam=$1
tool_info=$2


samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
check_variable "$tool_info:SAMTOOLS" $samtools

$samtools/samtools index $bam

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "indexing $bam took $DIFF seconds"
