#!/bin/sh

if [ $# != 2 ]
then
    echo "Usage : <output_dir> <run_info>";
else
    set -x
    echo `date`
    output_dir=$1
    run_info=$2
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
	genome_version=$(cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
	variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
    samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" " ")
	## jar script to add IGV, pathway, TIssue specificity and Gene Card link
   
	
	### merge per sample files to make merged report to be uploaded to TBB
	##Merge the unfiltered file
	cd $output_dir/Reports_per_Sample/
	mkdir -p $output_dir/Reports/
	## SNV
	for sample in $samples
	do
		ls $sample.SNV.xls >> list
	done	
	perl $script_path/union.snv.pl list $output_dir/Reports/SNV.xls
	rm list
	
	for sample in $samples
	do
	ls $sample.SNV.filtered.xls >> list
	done
	perl $script_path/union.snv.pl list $output_dir/Reports/SNV.filtered.xls
	rm list
	## INDEL
	for sample in $samples
	do
		ls $sample.INDEL.xls >> list
	done
	perl $script_path/union.indel.pl list $output_dir/Reports/INDEL.xls
	rm list
	
	for sample in $samples
	do
		ls $sample.INDEL.filtered.xls >> list
	perl $script_path/union.indel.pl list $output_dir/Reports/INDEL.filtered.xls
	rm list
	echo `date`
fi	
