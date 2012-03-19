#!/bin/sh

if [ $# != 3 ]
then
	echo -e "Usage: to plot coverage plot \n <input directory><output dir><run info >"
else
	set -x
	echo `date`
	input=$1
	output=$2
	run_info=$3
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2)
	tool=`echo "$tool" | tr "[A-Z]" "[a-z]"`
	CaptureKit=$( cat $tool_info | grep -w '^CAPTUREKIT' | cut -d '=' -f2 )
	master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
	samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" " " )
	if [ $tool == "whole_genome" ]
    then
        kit=$master_gene_file
    else
        kit=$CaptureKit
    fi    
	
	cd $input
	region=`awk '{sum+=$3-$2+1; print sum}' $kit | tail -1`
	
	Rscript $script_path/coverage_plot.r $region $samples
	mv $input/coverage.jpeg $output/Coverage.JPG
	echo `date`
	

fi	
