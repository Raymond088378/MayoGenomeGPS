#!/bin/sh

if [ $# != 2 ]
then
	echo "Usage: </path/to/output dir > <tool (exome/whole_genome)>"
else
	set -x
	echo `date`
	output_dir=$1
	tool=$2
	
	rm -R $output_dir/alignment
	rm -R $output_dir/fastq
	rm -R $output_dir/fastqc
	rm -R $output_dir/job_ids
	rm -R $output_dir/TempReports
	if [ $tool == "whole_genome" ]
	then
		rm -R $output_dir/Reports_per_Sample/plot
		rm -R $output_dir/Reports_per_Sample/temp
		rm -R $output_dir/cnv
		rm -R $output_dir/struct
	fi
	
	rm -R $output_dir/numbers
	rm -R $output_dir/OnTarget
	rm -R $output_dir/realign
	#if 
	#rm -R $output_dir/IGV_BAM
	echo `date`
fi	
		