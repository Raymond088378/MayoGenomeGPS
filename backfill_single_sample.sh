#!/bin/bash

### This should be submitted as an array job across number of samples, SGE_TASK_ID will determine the sample to backfill
### From the merged VCF

### Case 1: GATK EXOME ALL_SITES
### Case 2: GATK EXOME SOME_SITES
### Case 3: GATK WHOLE_GENOME
### Case 4: SNV_MIX, EXOME, ALL_SITES
### Case 5: SNV_MIX, EXOME, SOME_SITES
### Case 6:	SNV_MIX, WHOLE_GENOME 
### Case 6: BEAUTY EXOME
### Case 7: BEAUTY WHOLE_GENOME
			
### if [[ "$all_sites" == "YES" ]] $output/$sample.chr${chr}.all.vcf
###		else $output/$sample.chr${chr}.vcf

#!/bin/bash

####
# backfill_single_sample.sh
# accepts a merged variant file and backfills a single sample 
# 4/24/2013
# usage:
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
    
    exist_or_fail "$run_info:TOOL_INFO" "$tool_info" 
    exist_or_fail "$tool_info:WORKFLOW_PATH" "$script_path" 
    exist_or_fail "$run_info:MEMORY_INFO" "$memory_info"
    exist_or_fail "$run_info:BASE_OUTPUT_DIR" "$base_output"
    exist_or_fail "$tool_info:OUTPUT_FOLDER" "$run_num" 
    exist_or_fail "$tool_info:JAVA" "$javapath"
     
}

option=$1
shift
load_configuration_or_fail

if [[ "$option" == "GRID" ]]
then
	check_vars 0
	exist_or_fail "SGE_TASK_ID" $SGE_TASK_ID
 	samplenumber=$SGE_TASK_ID	 	 
elif [[ "$option" == "STANDALONE" ]]
	check_vars 1
	samplenumber=$1			
else 
	echo -e "Unrecognized or missing option.\n"
	print_usage
	exit 1
fi


