#!/bin/bash

####
# checkpoint.sh
# Generic Checkpointing Script for Workflows
# 4/23/2013
# usage:
#	
#	checkpoint.sh set /path/to/checkpoint/directory <ProcessName>
#		creates ProcessName.chkpt file in the checkpoint directory
#
# 	checkpoint.sh finish <groupID> <Stage> /path/to/checkpoint/directory <ProcessName>:<ProcessName>...
#		checks for presence of all ProcessName.chkpt files in the checkpoint directory and
#		calls AddSecondaryAnalysis <GroupID> <Stage> complete
#	
####

function print_usage()
{
	echo -e "checkpoint.sh -set /path/to/checkpoint/directory <ProcessName>"
	echo -e "\tcreates ProcessName.check file in the checkpoint directory\n"
	echo -e "checkpoint.sh -finish <groupID> <Stage> /path/to/checkpoint/directory <ProcessName>:<ProcessName>..."
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

if [[ $option == "-set" ]]
then
    ####
	# set a checkpoint 
    ####
	check_vars $# 3
	chkpt_path=$2
	raw_name=$3
	makedir_or_fail $output_path
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
	group=$2
	stage=$3
	chkpt_path=$4
	joblist=$5
	
	for job in `echo $joblist | tr ":" "\n"`
	do
	    files=$(ls $chkpt_path/$job.chkpt 2> /dev/null | wc -l)
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
