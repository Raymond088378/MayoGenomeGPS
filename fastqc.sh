#!/bin/sh

####
# fastqc.sh
# 4/24/2013
# script to run fastqc on a bam or fastq (zipped or unzipped) file
# baheti.saurabh@mayo.edu
####

### function
function check_variable()	{
	message=$1
	if [[ "$2" == "" ]] 
	then 
		echo "$message is not set."
		exit 1
	fi		
}	
if [ $# -le 2 ]
then
	echo -e "\nscript to run fastqc on a input file(bam or fastq)\nUsage: ./fastqc.sh </path/to/output folder> <fullpath/to/input file(bam or fastq)> <tool info file> </path/to/folder with previously computed results(optional)>"
	exit 1;
fi
START=$(date +%s)
output=$1
input_file=$2
tool_info=$3
if [ $4 ]
then
	FOLDER_FASTQC=$4
else
	FOLDER_FASTQC="NA"
fi	

## local variables
fastqc_path=$( cat $tool_info | grep -w '^FASTQC' | cut -d '=' -f2)
check_variable "$tool_info:FASTQC" $fastqc_path
script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
check_variable "$tool_info:WORKFLOW_PATH" $script_path


### check for extension of file 
ext=$(echo $input_file | sed 's/.*\.//')
if [ $ext != "gz" ]
then
	file1=$(echo $input_file | sed 's/\.[^\.]*$//')
else
	file1=$(echo $input_file | sed 's/\.[^\.]*$//'| sed 's/\.[^\.]*$//')
fi	    
            
### run fastqc depending on the flags
if [[ ${#FOLDER_FASTQC} -ne 0 && $FOLDER_FASTQC == "NA" ]]
then
	$fastqc_path/fastqc -o $output $input_file
	if [ -s $file.zip ]
	then
		rm $file.zip
	fi
elif [ $FOLDER_FASTQC != "NA" ]
then
	ln -s $FOLDER_FASTQC $file
else
	$script_path/errorlog.sh $FOLDER_FASTQC fastq.sh ERROR "doesn't have fastqc results"		
fi
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "fastqc for $input_file took $DIFF seconds"
