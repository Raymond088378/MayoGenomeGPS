#!/bin/bash

####
# checkpoint.sh
# Generic Checkpointing Script for Workflows
# 4/23/2013
# usage:
#	
#	checkpoint.sh set /path/to/checkpoint/directory <ProcessName>
#		creates ProcessName.check file in the checkpoint directory
#
# 	checkpoint.sh finish <groupID> <Stage> /path/to/checkpoint/directory <ProcessName>:<ProcessName>...
#		checks for presence of all ProcessName.check files in the checkpoint directory and
#		calls AddSecondaryAnalysis <GroupID> <Stage> complete
#	
####

function print_usage()
{
	echo -e "\n"
	echo -e "checkpoint.sh set /path/to/checkpoint/directory <ProcessName>"
	echo -e "\tcreates ProcessName.check file in the checkpoint directory"
	echo -e "\n"
	echo -e "checkpoint.sh finish <groupID> <Stage> /path/to/checkpoint/directory <ProcessName>:<ProcessName>..."
	echo -e "\tchecks for presence of all ProcessName.check files in the checkpoint directory and"
	echo -e "\tcalls AddSecondaryAnalysis <GroupID> <Stage> complete."
	echo -e "\n"
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

function build_unique_name() 
{
	local name=$1
	
	if [[ $JOB_ID ]]
	then
		name=$name.$JOB_ID
	fi
	
	if [[ $SGE_TASK_ID ]]
	then
		name=$name.$SGE_TASK_ID
	fi
	
	_uniquename_=$name
}



typeset -u option=$1
echo $option

if [[ $option == "SET" ]]
then
    ####
	# set a checkpoint 
    ####
	check_vars $# 3
	chkpt_path=$2
	raw_name=$3
	makedir_or_fail $output_path
	build_unique_name $raw_name
	name=$_uniquename_
	touch $chkpt_path/$name.chkpt
    exit 0
elif [[ $option == "FINISH" ]]
then
	check_vars $# 5
	group=$2
	stage=$3
	chkpt_path=$4
	joblist=$5
	
	for job in `echo $joblist | tr ":" "\n"`
	do
	    files=$(ls $chkpt_path/$job.*.chkpt 2> /dev/null | wc -l)
	    echo $files found
	    if [[ $files != 0 ]]
	    then
		echo Found $job
	    else
		echo Missing $job 
	    fi
	done
    exit 0
else
	echo -e "Command not recognized.\n"
    print_usage
    exit 1
fi
