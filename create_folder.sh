#!/bin/bash

if [ $# != 1 ]
then
	echo -e "Usage: wrapper to create folder structure\nUsage: ./create_folder.sh <run info file>"
	exit 1;
fi
	echo `date`
	run_info=$1
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )	
	output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	
	if [ -d $output/$PI ]
	then
		echo "WARNING : $PI folder exists"
	else
		mkdir -p $output/$PI
	fi

	if [ -d $output/$PI/$tool ]
	then 	
		echo "WARNING : $tool analysis folder exists"
	else
		mkdir -p $output/$PI/$tool
	fi

	if [ -d $output/$PI/$tool/$run_num ]
	then
		echo "ERROR : $run_num folder exists"
		touch $output/$PI/$tool/$run_num/folder_exist.log
		exit 1;
	else 
		mkdir -p $output/$PI/$tool/$run_num
		output_dir=$output/$PI/$tool/$run_num
	fi
	
	mkdir -p $output_dir/logs
    mkdir -p $output_dir/config

	if [ $analysis == "external" -o $analysis == "mayo" -o $analysis == "alignment" ]
	then
		mkdir -p $output_dir/qc
		mkdir -p $output_dir/qc/fastqc
	fi
	if [[ $analysis != "variant" ]]
	then
		mkdir -p $output_dir/alignment
	fi
	mkdir -p $output_dir/IGV_BAM
	if [ $analysis != "alignment" ]
	then
    	mkdir -p $output_dir/realign      
    	mkdir -p $output_dir/variants
    	if [ $tool == "whole_genome" ]
    	then
        	mkdir -p $output_dir/variants/cnv
        	mkdir -p $output_dir/variants/struct
        	mkdir -p $output_dir/variants/circos
    	fi
	fi
	
	if [ $analysis != "alignment" ]
    then
        mkdir -p $output_dir/Reports_per_Sample
        mkdir -p $output_dir/Reports
    fi
    mkdir -p $output_dir/numbers
    echo `date`
fi	
