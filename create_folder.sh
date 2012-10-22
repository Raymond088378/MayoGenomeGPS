#!/bin/bash

if [ $# != 1 ]
then
	echo "Usage: wrapper to create folder structure\n <run info file>"
else
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
	mkdir -p $output_dir/size
        mkdir -p $output_dir/config

	if [[ $analysis != "annotation"  && $analysis != "ontarget" ]]
	then
		if [ $analysis == "external" -o $analysis == "mayo" -o $analysis == "alignment" ]
		then
			mkdir -p $output_dir/fastq
			mkdir -p $output_dir/fastqc
		fi
		if [[ $analysis != "variant" ]]
		then
			mkdir -p $output_dir/alignment
		fi
		mkdir -p $output_dir/IGV_BAM
		if [ $analysis != "alignment" ]
        then
            if [[ $analysis != "annotation" && $analysis != "ontarget"  ]]
            then
                mkdir -p $output_dir/realign      
            fi
            if [ $analysis != "annotation" ]
            then
                mkdir -p $output_dir/variants
            fi
            if [ $tool == "whole_genome" ]
            then
                mkdir -p $output_dir/cnv
                mkdir -p $output_dir/struct
                mkdir -p $output_dir/circos
            fi
        fi
    fi
	
	if [ $analysis != "alignment" ]
    then
        mkdir -p $output_dir/OnTarget
        mkdir -p $output_dir/annotation
        output_annot=$output_dir/annotation
		mkdir -p $output_annot/SIFT
        mkdir -p $output_annot/SNPEFF
        mkdir -p $output_annot/POLYPHEN
        mkdir -p $output_dir/TempReports
        mkdir -p $output_dir/Reports_per_Sample
        mkdir -p $output_dir/Reports
    fi
    mkdir -p $output_dir/numbers
    echo `date`
fi	
