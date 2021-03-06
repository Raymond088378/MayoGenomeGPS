#!/bin/bash
###
# copy the config files and other stable files
###

if [ $# != 2 ]
then
	echo "Usage: wrapper to copy the config files\n <output dir ><run_info>"
	exit 1;
fi
    echo `date`	
    output_dir=$1
    run_info=$2
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
    
	config=$output_dir/config
	cp $tool_info $config/tool_info.txt
	cp $sample_info $config/sample_info.txt
	cp $memory_info $config/memory_info.txt
    
    cp $script_path/${tool}_workflow.png $output_dir/${tool}_workflow.png
    cp $script_path/IGV_Setup.doc $output_dir/IGV_Setup.doc
    echo `date`
	