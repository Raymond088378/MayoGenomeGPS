#!/bin/bash

if [ $# != 4 ]
then
	echo -e "script to upload the status to the dash board when megng diffrent flowcell\nUsage: ./manual_dashboard.sh <','seperated flowcells> <type of data(Exome/WholeGenome)><status <Beginning/Complete/Delivered> <full/path/to/tool_info file>"
else
	flowcells=`echo $1 | tr "," " "`
	type=$2
	status=$3
	tool_info=$4
	
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	for flow in $flowcells
	do
		id=`echo $flow |awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
		$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -c -f $id -r $flow -s $status -a $type
	done
fi
	