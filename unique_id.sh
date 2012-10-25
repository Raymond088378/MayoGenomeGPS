#!/bin/bash

if [ $# != 1 ]
then
	echo -e "script to get the unique id for execution of the workflow\nUsage: ./unique_id.sh </path/to/run info file>"
else
	echo `date`	
	run_info=$1
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	type=$( cat $run_info | grep -w '^TOOL' | cut -d '=' -f2|tr "[a-z]" "[A-Z]")
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	TO=$USER
	cat $run_info | grep -w -v -E '^IDENTIFICATION_NUMBER' > $run_info.tmp
	mv $run_info.tmp $run_info
	
	unique_id=`$java/java -Xmx1g -jar $script_path/AddGPSMetadata.jar -p $script_path/AddGPSMetadata.properties -t $type -a begin -I -r $run_num -s $run_num -u $TO`
	if echo $unique_id | egrep -q '^[0-9]+$'
	then
		echo -e "IDENTIFICATION_NUMBER=$unique_id" >> $run_info
	else
		echo "ERROR : unique identification for the workflow was not generated"
		exit 1;
	fi	
	echo `date`
fi	
