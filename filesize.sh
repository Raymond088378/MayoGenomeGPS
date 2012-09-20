#!/bin/bash

if [ $# != 6 ]
then
	echo "usage:<analysis type><sample name><filename><job id ><size of the file><run info >"
else
	echo `date`
	analysis=$1
	sample=$2
	dirname=$3
	filename=$4
	job=$5
	run_info=$6
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	type=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	identify=$( cat $run_info | grep -w '^IDENTIFICATION_NUMBER' | cut -d '=' -f2)
	output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	out=$output/$PI/$tool/$run_num
	TO=`id |awk -F '(' '{print $2}' | cut -f1 -d ')'`
	size=`du -b $dirname/$filename | sed 's/\([0-9]*\).*/\1/'`
	id=$SGE_TASK_ID
	if [ $id ]
	then
		file=$out/size/filesize.$job.$id.csv
	else
		file=$out/size/filesize.$job.csv
	fi	
	$java/java -Xmx32M -jar $script_path/AddGPSMetadata.jar -p $script_path/AddGPSMetadata.properties -S $identify -t $type -a $analysis -b $size -j $job -r $run_num -s $sample -n $filename -u $TO -F $file
	echo `date`
fi