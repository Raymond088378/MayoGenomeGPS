#!/bin/sh

if [ $# != 3 ]
then
	echo "Usage:";
	echo "1. input directory containing AllSamples_GeneCount.noreads";
	echo "2. output directory for plot";
	echo "3. script path";
else
	set -x
	echo `date`
	input_dir=$1
	output_dir=$2
	script_path=$3
	
	cd $output_dir
	cat $input_dir/AllSamples_GeneCount.noreads |  awk 'gsub(/_GeneCount.noreads/, "")1' | cut -f4 -d " " > $output_dir/samples
	echo "Count" >> $output_dir/numbers
	cat $input_dir/AllSamples_GeneCount.noreads |  awk 'gsub(/_GeneCount.noreads/, "")1' | cut -f3 -d " " >> $output_dir/numbers
	Rscript $script_path/QC_mRNASeq_noreadsmapped_Rplot.r $output_dir/numbers $output_dir/samples
	rm $output_dir/numbers $output_dir/samples
fi