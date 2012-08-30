#!/bin/bash

if [ $# != 6 ]
then
	echo "usage:<analysis type><sample name><filename><job id ><size of the file><run info >"
else
	echo `date`
	analysis=$1
	sample=$2
	filename=$3
	job=$4
	size=$5
	run_info=$6
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	type=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	identify=$( cat $run_info | grep -w '^IDENTIFICATION_NUMBER' | cut -d '=' -f2)
	TO=`id |awk -F '(' '{print $2}' | cut -f1 -d ')'`
	
	$java/java -jar $script_path/AddGPSMetadata.jar -p $script_path/AddGPSMetadata.properties -S $identify -t $type -a $analysis -b $size -j $job -r $run_num -s $sample -n $filename -u $TO
	echo `date`
fi