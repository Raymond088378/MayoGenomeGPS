#!/bin/bash

if [ $# != 4 ]
then
	echo -e "In buit QC script to send an email if something fails\nUsage: ./email.sh <filename> <error message><previous script><runinfo>"
else	
	echo `date`
	file=$1
	message=$2
	previous=$3
	run_info=$4
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	workflow=$( cat $run_info | grep '^TOOL=' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
	job_name=$JOB_NAME 
	job_id=$JOB_ID
	TO=$USER
	email=`finger $USER | awk -F ';' '{print $2}'`
	SUB="$tool workflow In error State for RunID: ${run_num}"
	date=`date`
	SGE=$SGE_TASK_ID
	if [ ! $SGE ]
	then
		SGE="-"
	fi	
	MESG="Date: $date\n============================\n\nFilename: $file\nError: $message from the predecessor job\nScriptToCheck: $previous\nJobName: $job_name\nJobId: $job_id\nArrayJobId: $SGE\n\nPlease fix the error and delete the $file.fix.log file so that job can resume properly\n\nCourtesy: $workflow $version"
	## send the completion email
	echo -e "$MESG" | mailx -s "$SUB" "$email" 
	echo `date`
fi	
	
