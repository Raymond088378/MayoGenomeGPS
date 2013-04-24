#!/bin/bash

####
# checkpoint.sh
# Generic Checkpointing Script for Workflows
# 4/23/2013
# usage:
#	
#	checkpoint.sh set /path/to/run_info <ProcessName>
#		creates ProcessName.chkpt file in the checkpoint directory
#
# 	checkpoint.sh finish /path/to/run_info <groupID> <Stage> <ProcessName>:<ProcessName>...
#		checks for presence of all ProcessName.chkpt files in the checkpoint directory and
#		calls AddSecondaryAnalysis <GroupID> <Stage> complete
#	
####

function print_usage()
{
	echo -e "checkpoint.sh -set /path/to/run_info <ProcessName>"
	echo -e "\tcreates ProcessName.check file in the checkpoint directory. \n"
	echo -e "checkpoint.sh -finish /path/to/run_info <groupID> <Stage> <ProcessName>:<ProcessName>..."
	echo -e "\tchecks for presence of all ProcessName.check files in the checkpoint directory and"
	echo -e "\tcalls AddSecondaryAnalysis <GroupID> <Stage> complete. \n"
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

## makedir_or_fail <directory>
## 		check if <directory> exists, and if not, create it
## 		if directory can not be created, exit script with an error
function makedir_or_fail () 
{
	if [ ! -d $1 ]
	then
		mkdir $1
		if [ $? != 0 ]
		then
			echo Failed to create directory $1.
			exit 1
		fi
	fi
}

## exist_or_fail "error message" variable1 variable2 ... variableN
## 		check if variables1..N are set; if not, exit script with error message
function exist_or_fail ()
{
	local message=$1;
	local found="";
	shift
	while (( "$#" )); do
		if [[ "$1" == "" ]] 
		then 
			echo "$message is not set."
			exit 1
		fi
		shift
	done
}

function load_run_info_or_fail ()
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
    local javamem=$( grep -w '^AddSecondaryAnalysis_JVM' $memory_info | cut -d '=' -f2)
    
    exist_or_fail "$run_info:TOOL_INFO" "$tool_info" 
    exist_or_fail "$tool_info:WORKFLOW_PATH" "$script_path" 
    exist_or_fail "$run_info:MEMORY_INFO" "$memory_info"
    exist_or_fail "$run_info:BASE_OUTPUT_DIR" "$base_output"
    exist_or_fail "$tool_info:OUTPUT_FOLDER" "$run_num" 
    exist_or_fail "$tool_info:JAVA" "$javapath"
	exist_or_fail "$memory_info:AddSecondaryAnalysis_JVM" "$javamem"
     
    addsecondary="$javapath/java $javamem $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties"    
	chkpt_path=$base_output/$run_num/.runstatus/
}

option=$1

if [[ $option == "-set" ]]
then
    ####
	# set a checkpoint 
    ####
	check_vars $# 3
	
	load_run_info_or_fail $2
	raw_name=$3
	makedir_or_fail $chkpt_path
	if [ ! -f $chkpt_path/$raw_name.chkpt ]
	then
		touch $chkpt_path/$raw_name.chkpt
		exit 0
	else
		echo "Warning - Checkpoint File for $raw_name already present."
		exit 2
	fi
elif [[ $option == "-finish" ]]
then
	check_vars $# 5
	load_run_info_or_fail $2
	group=$3
	stage=$4
	joblist=$5
	
	missing=""
	count=0
	
	for job in `echo $joblist | tr ":" "\n"`
	do
	    files=$(ls $chkpt_path/$job.chkpt 2> /dev/null | wc -l)
	    if [[ $files == 0 ]]
	    then 
	   		## missing a checkpoint
			(( count++ ))
	   		if [[ "$missing" == "" ]] 
	    	then
				missing=$job
			else
				missing=$missing:$job
			fi
		else
			## found this checkpoint, so remove the file
			rm -f $chkpt_path/$job.chkpt
		fi
	done
	
	if [[ count != 0 ]]
	then
		echo Missing $missing
		exit 1
	else
		## Call AddSecondaryAnalysis and set Stage Complete
    	exit 0
	fi
	
else
	echo -e "Command not recognized.\n"
    print_usage
    exit 1
fi
