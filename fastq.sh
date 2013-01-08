#!/bin/bash

if [ $# != 5 ]
then
    echo -e "SCRIPT  to run fastqc on zipped or unzipped fastq\nUsage: ./fastq.sh <input fastq file> </path/to/input dir> </path/to/output dir> </path/to/run info> </path/to/FASTQC directory>"
else
    set -x
    echo `date`
    read=$1
    input=$2
    output=$3
    run_info=$4
    fastqc_dir=$5

    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    fastqc_path=$( cat $tool_info | grep -w '^FASTQC' | cut -d '=' -f2)
    FASTQC=$( cat $run_info | grep -w '^FASTQC' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    FOLDER_FASTQC=$( cat $run_info | grep -w '^FOLDER_FASTQC' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    
    if [ ! -s $input/$read ]
    then
        $script_path/errorlog.sh $input/$read fastq.sh ERROR "not found"
        exit 1;
    fi
    
    ext=$(echo $read | sed 's/.*\.//')
    if [ $ext != "gz" ]
    then
		file1=$(echo $read | sed 's/\.[^\.]*$//')
    else
		file1=$(echo $read | sed 's/\.[^\.]*$//'| sed 's/\.[^\.]*$//')
    fi	
    ## make a soft link fastq's
    ln -s $input/$read $output/$read	
    
    ##FASTQC
    if [ $FASTQC == "YES" ]
    then
        $fastqc_path/fastqc -o $fastqc_dir/ $output/$read
		rm $fastqc_dir/${file1}_fastqc.zip
    else
        if [ -f $FOLDER_FASTQC/${file1}_fastqc ]
		then
			ln -s $FOLDER_FASTQC/${file1}_fastqc $fastqc_dir/
		else
			 $script_path/errorlog.sh $FOLDER_FASTQC fastq.sh ERROR "doesn't have fastqc results"
			exit 1;
		fi	
    fi		
    echo `date`
fi	