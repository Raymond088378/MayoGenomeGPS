#!/bin/bash

if [ $# != 2 ]
then
	echo -e "script to upload the runs on table browser\nUsage: ./upload_TTB.sh </path/to/output folder></aothto run info file>"
else
	echo `date`
	output_dir=$1
	run_info=$2
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	upload_tb=$( cat $tool_info | grep -w '^UPLOAD_TABLEBROWSER' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	
	if [[ $upload_tb == "YES"  && $analysis != "alignment" ]]
	then
		PI_LANID=$( echo $PI | cut -d '_' -f 3 )
		mem=$( cat $memory_info | grep -w '^TREATUploader_JVM' | cut -d '=' -f2)
		$java/java $mem -jar $script_path/TREATUploader.jar -n $PI_LANID -u $run_num -i $output_dir/Reports/INDEL.xls -s $output_dir/Reports/SNV.xls -r $run_num
		echo -e "Variants uploaded to TableBrowser" >> $output_dir/log.txt
	else
		echo -e "Variants Not uploaded to TableBrowser" >> $output_dir/log.txt
	fi		
	echo `date`
fi	