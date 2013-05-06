#!/bin/bash

####
# backfill_ss_merge.sh
# creates a merged variant file containing all samples
# 4/24/2013
# usage:
# This should be submitted as single task to create a merged vcf for backfilling
####

function print_usage()
{
	echo "Usage: backfill_ss_merge.sh <path/to/run_info.txt>"
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
	run_info=$1
		
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
    
    local PI=$( grep -w '^PI' $run_info | cut -d '=' -f2)
    local tool=$( grep -w '^TYPE' $run_info | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
    local run_num=$( grep -w '^OUTPUT_FOLDER' $run_info | cut -d '=' -f2)
    
    
    exist_or_fail "$run_info:TOOL_INFO" "$tool_info" 
    exist_or_fail "$tool_info:WORKFLOW_PATH" "$script_path" 
    exist_or_fail "$run_info:MEMORY_INFO" "$memory_info"
    exist_or_fail "$run_info:BASE_OUTPUT_DIR" "$base_output"
    exist_or_fail "$tool_info:OUTPUT_FOLDER" "$run_num" 
    exist_or_fail "$tool_info:JAVA" "$javapath"
	exist_or_fail "$run_info:SAMPLENAMES" "$samples";
	exist_or_fail "$run_info:PI" "$PI";
	exist_or_fail "$run_info:TYPE" "$tool";
	exist_or_fail "$run_info:OUTPUT_FOLDER" "$run_num";
	
	local numsamples=$(grep -w '^SAMPLENAMES' $run_info | cut -d '=' -f2 | tr ":" "\n" | wc -l)

	local i=1	
	for sample in `echo $samples | tr ":" "\n"`
	do
		sampleArray[$i]=$sample
		let i=i+1
	done
	
	output_variant=$output/$PI/$tool/$run_num/variants/
	output_IGV=$output/$PI/$tool/$run_num/IGV_BAM/
		
}

check_vars $# 1
load_configuration_or_fail $1

### Merge
input_var=""
for i in $(seq 1 ${#sampleArray[@]})
do
	sample=${sampleArray[$i]}
	vcf=$sample.variants.final.vcf
	input_var="${input_var} -V $output_variant/$vcf"
done

$script_path/combinevcf.sh "$input_var" $output_variant/backfill.final.vcf $run_info no

if [ ! -s $output_variant/backfill.final.vcf ]
then
	$script_path/errorlog.sh $output_variant/backfill.final.vcf backfill_ss_merge.sh ERROR "Error: failed to create merged all sample variant file."
    exit 1;
fi




