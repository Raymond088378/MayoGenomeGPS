#!/bin/bash

if [ $# != 2 ]
then
    echo -e "script to merge the per sample report\nUsage : <output_dir> <run_info>";
else
#    set -x
    echo `date`
    output_dir=$1
    run_info=$2
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
    genome_version=$(cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
    samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" " ")
    multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2 | tr ":" " ")
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2 |tr "[a-z]" "[A-Z]")
    snv_caller=$( cat $run_info | grep -w '^SNV_CALLER' | cut -d '=' -f2)
    ### merge per sample files to make merged report to be uploaded to TBB
    ##Merge the unfiltered file
    cd $output_dir/Reports_per_Sample/
    mkdir -p $output_dir/Reports/
    ## SNV
    if [ $multi_sample == "NO" ]
	then
		for sample in $samples
		do
			if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
			then
				if [ ! -s $sample.SNV.xls ]
				then
					$script_path/email.sh $sample.SNV.xls "doesn't exist" "sample_report.sh" $run_info
					touch $sample.SNV.xls.fix.log
					$script_path/wait.sh $sample.SNV.xls.fix.log 
				fi
				ls $sample.SNV.xls >> list.snv
				if [ ! -s $sample.SNV.filtered.xls ]
				then
					$script_path/email.sh $sample.SNV.filtered.xls "doesn't exist" "sample_report.sh" $run_info
					touch $sample.SNV.filtered.xls.fix.log
					$script_path/wait.sh $sample.SNV.filtered.xls.fix.log 
				fi
				ls $sample.SNV.filtered.xls >> list.filter.snv
			fi
			if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
			then
				if [ ! -s $sample.INDEL.xls ]
				then
					$script_path/email.sh $sample.INDEL.xls "doesn't exist" "sample_report.sh" $run_info
					touch $sample.INDEL.xls.fix.log
					$script_path/wait.sh $sample.INDEL.xls.fix.log 
				fi
				ls $sample.INDEL.xls >> list.indel
				if [ ! -s $sample.INDEL.filtered.xls ]
				then
					$script_path/email.sh $sample.INDEL.filtered.xls "doesn't exist" "sample_report.sh" $run_info
					touch $sample.INDEL.filtered.xls.fix.log
					$script_path/wait.sh $sample.INDEL.filtered.xls.fix.log 
				fi
				ls $sample.INDEL.filtered.xls >> list.filter.indel
			fi
		done	
		if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
		then
			perl $script_path/union.snv.pl list.snv single $output_dir/Reports/SNV.xls
			perl $script_path/union.snv.pl list.filter.snv single $output_dir/Reports/SNV.filtered.xls
			rm list.snv list.filter.snv
		fi
		if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
		then
			perl $script_path/union.indel.pl list.indel single $output_dir/Reports/INDEL.xls
			perl $script_path/union.indel.pl list.filter.indel single $output_dir/Reports/INDEL.filtered.xls
			rm list.indel list.filter.indel
		fi
	else
		for group in $groups
		do
			if [ ! -s $group.SNV.xls ]
			then
				$script_path/email.sh $group.SNV.xls "doesn't exist" "sample_report.sh" $run_info
				touch $group.SNV.xls.fix.log
				$script_path/wait.sh $group.SNV.xls.fix.log 
			fi
			ls $group.SNV.xls >> list.snv
			if [ ! -s $group.SNV.filtered.xls ]
			then
				$script_path/email.sh $group.SNV.filtered.xls "doesn't exist" "sample_report.sh" $run_info
				touch $group.SNV.filtered.xls.fix.log
				$script_path/wait.sh $group.SNV.filtered.xls.fix.log 
			fi
			ls $group.SNV.filtered.xls >> list.filter.snv
			if [ ! -s $group.INDEL.xls ]
			then
				$script_path/email.sh $group.INDEL.xls "doesn't exist" "sample_report.sh" $run_info
				touch $group.INDEL.xls.fix.log
				$script_path/wait.sh $group.INDEL.xls.fix.log 
			fi
			ls $group.INDEL.xls >> list.indel
			if [ ! -s $group.SNV.filtered.xls ]
			then
				$script_path/email.sh $group.INDEL.filtered.xls "doesn't exist" "sample_report.sh" $run_info
				touch $group.INDEL.filtered.xls.fix.log
				$script_path/wait.sh $group.INDEL.filtered.xls.fix.log 
			fi
			ls $group.INDEL.filtered.xls >> list.filter.indel
		done
		perl $script_path/union.snv.pl list.snv multi $output_dir/Reports/SNV.xls
		perl $script_path/union.snv.pl list.filter.snv multi $output_dir/Reports/SNV.filtered.xls
		perl $script_path/union.indel.pl list.indel multi $output_dir/Reports/INDEL.xls
		perl $script_path/union.indel.pl list.filter.indel multi $output_dir/Reports/INDEL.filtered.xls	
		rm list.snv list.filter.snv list.indel list.filter.indel	
    	### Merge the TUMOR files
    	for group in $groups
		do
			if [ ! -s TUMOR.$group.SNV.xls ]
			then
				$script_path/email.sh TUMOR.$group.SNV.xls "doesn't exist" "sample_report.sh" $run_info
				touch TUMOR.$group.SNV.xls.fix.log
				$script_path/wait.sh TUMOR.$group.SNV.xls.fix.log 
			fi
			ls TUMOR.$group.SNV.xls >> list.snv
			if [ ! -s TUMOR.$group.SNV.filtered.xls ]
			then
				$script_path/email.sh TUMOR.$group.SNV.filtered.xls "doesn't exist" "sample_report.sh" $run_info
				touch TUMOR.$group.SNV.filtered.xls.fix.log
				$script_path/wait.sh TUMOR.$group.SNV.filtered.xls.fix.log 
			fi
			ls TUMOR.$group.SNV.filtered.xls >> list.filter.snv
			if [ ! -s TUMOR.$group.INDEL.xls ]
			then
				$script_path/email.sh TUMOR.$group.INDEL.xls "doesn't exist" "sample_report.sh" $run_info
				touch TUMOR.$group.INDEL.xls.fix.log
				$script_path/wait.sh TUMOR.$group.INDEL.xls.fix.log 
			fi
			ls TUMOR.$group.INDEL.xls >> list.indel
			if [ ! -s $group.SNV.filtered.xls ]
			then
				$script_path/email.sh TUMOR.$group.INDEL.filtered.xls "doesn't exist" "sample_report.sh" $run_info
				touch TUMOR.$group.INDEL.filtered.xls.fix.log
				$script_path/wait.sh TUMOR.$group.INDEL.filtered.xls.fix.log 
			fi
			ls TUMOR.$group.INDEL.filtered.xls >> list.filter.indel
		done
		perl $script_path/union.snv.pl list.snv multi $output_dir/Reports/TUMOR.SNV.xls
		perl $script_path/union.snv.pl list.filter.snv multi $output_dir/Reports/TUMOR.SNV.filtered.xls
		perl $script_path/union.indel.pl list.indel multi $output_dir/Reports/TUMOR.INDEL.xls
		perl $script_path/union.indel.pl list.filter.indel multi $output_dir/Reports/TUMOR.INDEL.filtered.xls	
		rm list.snv list.filter.snv list.indel list.filter.indel	
	fi

	### Adding Beauty Anntation Module
	if [ $snv_caller == "BEAUTY_EXOME" ]
	then
		## Script to create additional annotation file, with pharmicogentic emphasis
		$script_path/PharmacoAnnotModule.sh $output_dir/Reports/SNV.filtered.xls $run_info
	fi

    echo `date`
fi	
