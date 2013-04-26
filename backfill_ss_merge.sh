#!/bin/bash

####
# backfill_ss_merge.sh
# creates a merged variant file containing all samples
# 4/24/2013
# usage:
# This should be submitted as single task to create a merged vcf for backfilling
#
####

function print_usage()
{
	
}

function check_vars()
{
	if [[ $1 != $2 ]] 
	then
		echo -e "\nExpected $2 entries on the command line, got $1.\n"
		print_usage
		exit 1
	fi		
}

## exist_or_fail "variable name" variable
## 		check if variable is set; if not, exit script with error message
function exist_or_fail ()
{
	local message=$1;
	local found="";
	shift
	if [[ "$1" == "" ]] 
	then 
		echo "$message is not set."
		exit 1
	fi
}

function load_configuration_or_fail ()
{
	local run_info=$1
		
	if [[ ! -f $run_info ]]
	then
		echo "Failed to find run_info. Exiting."
		exit 1
	fi
	
	## Load directory names from the run_info and associated configuration files
	local tool_info=$( grep -w '^TOOL_INFO' $run_info | cut -d '=' -f2 )
	local script_path=$( grep -w '^WORKFLOW_PATH' $tool_info | cut -d '=' -f2)
	local memory_info=$( grep -w '^MEMORY_INFO' $run_info | cut -d '=' -f2)
    local base_output=$( grep -w '^BASE_OUTPUT_DIR' $run_info | cut -d '=' -f2)
	local run_num=$( grep -w '^OUTPUT_FOLDER' $run_info | cut -d '=' -f2)
	local javapath=$( grep -w '^JAVA' $tool_info | cut -d '=' -f2)
    local samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2)
    
    exist_or_fail "$run_info:TOOL_INFO" "$tool_info" 
    exist_or_fail "$tool_info:WORKFLOW_PATH" "$script_path" 
    exist_or_fail "$run_info:MEMORY_INFO" "$memory_info"
    exist_or_fail "$run_info:BASE_OUTPUT_DIR" "$base_output"
    exist_or_fail "$tool_info:OUTPUT_FOLDER" "$run_num" 
    exist_or_fail "$tool_info:JAVA" "$javapath"
	exist_or_fail "$run_info:SAMPLENAMES" "$samples";
	
	local numsamples=$(grep -w '^SAMPLENAMES' $run_info | cut -d '=' -f2 | tr ":" "\n" | wc -l)
	
	for sample in `echo $samples | tr ":" "\n"`
	do
	
			
	done      
}

load_configuration_or_fail





