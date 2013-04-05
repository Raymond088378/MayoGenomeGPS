#!/bin/bash

if [ $# != 3 ]
then
	echo -e "script to run fastqc on a input file(bam or fastq)\nUsage: </path/to/output folder> <fullpath/to/input file(bam or fastq)> <tool infofile>"
	exit 1;
fi

set -x
echo `date`
output=$1
input_file=$2
tool_info=$3

fastqc_path=$( cat $tool_info | grep -w '^FASTQC' | cut -d '=' -f2)
FOLDER_FASTQC=$( cat $tool_info | grep -w '^FOLDER_FASTQC' | cut -d '=' -f2 )
script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )

ext=$(echo $input_file | sed 's/.*\.//')
if [ $ext != "gz" ]
then
	file1=$(echo $input_file | sed 's/\.[^\.]*$//')
else
	file1=$(echo $input_file | sed 's/\.[^\.]*$//'| sed 's/\.[^\.]*$//')
fi	    
            
if [[ ${#FOLDER_FASTQC} -ne 0 && $FOLDER_FASTQC == "NA" ]]
then
	$fastqc_path/fastqc -o $output $input_file
	rm $file.zip
elif [ $FOLDER_FASTQC != "NA" ]
then
	ln -s $FOLDER_FASTQC $file
else
	$script_path/errorlog.sh $FOLDER_FASTQC fastq.sh ERROR "doesn't have fastqc results"		
fi
echo `date`	