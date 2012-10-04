#!/bin/bash

if [ $# != 5 ]
then
	echo -e "Usage: In buit QC script \n <filename> <error message><job name><job id><runinfo>"
else	
	echo `date`
	file=$1
	message=$2
	job_name=$3
	job_id=$4
	run_info=$5
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	workflow=$( cat $run_info | grep '^TOOL=' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)

	TO=`id |awk -F '(' '{print $2}' | cut -f1 -d ')'`
	SUB="$tool workflow In error State for RunID: ${run_num}"
	date=`date`
	SGE=$SGE_TASK_ID
	if [ ! $SGE ]
	then
		SGE="-"
	fi	
	MESG="Date: $date\nFilename: $file\nError: $message from the predecessor job\njobname: $job_name\njobid: $job_id\nArrayjobid: $SGE\nPlease fix the error and delete the $file.fix.log file so that job can resume properly\n\nCourtesy: $workflow $version"
	## send the completion email
	echo -e "$MESG" | mailx -v -s "$SUB" "$TO" 
	echo `date`
fi	
	