#!/bin/sh

if [ $# != 2 ]
then
	echo "Usage: wrapper to create folder structure\n </path/to/output folder> <run info file>"
else
	set -x
	echo `date`
	output_dir=$1
	run_info=$2
	
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )	
	mkdir $output_dir/job_ids

	if [ $analysis != "annotation" ]
	then
		if [ $analysis == "external" -o $analysis == "mayo" ]
		then
			mkdir $output_dir/fastq
			mkdir $output_dir/fastqc
		fi
		mkdir $output_dir/alignment
		mkdir $output_dir/realign
		mkdir $output_dir/IGV_BAM
		mkdir $output_dir/variants
		if [ $tool == "whole_genome" ]
		then
			mkdir $output_dir/cnv
			mkdir $output_dir/struct
			mkdir $output_dir/circos
		fi
	fi
	
	mkdir $output_dir/OnTarget
	mkdir $output_dir/numbers
	mkdir $output_dir/annotation
	output_annot=$output_dir/annotation
	mkdir $output_annot/SIFT
	mkdir $output_annot/SSEQ
	mkdir $output_dir/TempReports
	mkdir $output_dir/Reports_per_Sample
	mkdir $output_dir/Reports
	echo `date`
fi	
