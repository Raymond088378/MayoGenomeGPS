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
    multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2 | tr ":" " ")
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2 |tr "[a-z]" "[A-Z]")
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
				ls $sample.SNV.xls >> list.snv
				ls $sample.SNV.filtered.xls >> list.filter.snv
			fi
			if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
			then
				ls $sample.INDEL.xls >> list.indel
				ls $sample.INDEL.filtered.xls >> list.filter.indel
			fi
		done	
		if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
		then
			perl $script_path/union.snv.pl list.snv $output_dir/Reports/SNV.xls
			perl $script_path/union.snv.pl list.filter.snv $output_dir/Reports/SNV.filtered.xls
			rm list.snv list.filter.snv
		fi
		if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
		then
			perl $script_path/union.indel.pl list.indel $output_dir/Reports/INDEL.xls
			perl $script_path/union.indel.pl list.filter.indel $output_dir/Reports/INDEL.filtered.xls
			rm list.indel list.filter.indel
		fi
	else
		for group in $groups
		do
			for sam in `cat $sample_info | grep "^$group" | cut -f2 -d '=' | tr "\t" " "`
			do
				ls $group.$sam.SNV.xls >> list.snv
				ls $group.$sam.SNV.filtered.xls >> list.filter.snv
				ls $group.$sam.INDEL.xls >> list.indel
				ls $group.$sam.INDEL.filtered.xls >> list.filter.indel
			done
			perl $script_path/union.snv.pl list.snv $output_dir/Reports/$group.SNV.xls
			perl $script_path/union.snv.pl list.filter.snv $output_dir/Reports/$group.SNV.filtered.xls
			perl $script_path/union.indel.pl list.indel $output_dir/Reports/$group.INDEL.xls
			perl $script_path/union.indel.pl list.filter.indel $output_dir/Reports/$group.INDEL.filtered.xls	
			rm list.snv list.filter.snv list.indel list.filter.indel
		done	
    fi
	
    if [ $multi_sample == "YES" ]
    then
		for i in $groups
		do
			sampleNames=$( cat $sample_info| grep -w "^$i" | cut -d '=' -f2 )
			k=1
			for sample in $sampleNames
			do
				sampleArray[$k]=$sample
				let k=k+1
			done
			for j in $(seq 2 ${#sampleArray[@]})
			do  
				tumor=${sampleArray[$j]}
				if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
				then
					ls TUMOR.$i.$tumor.SNV.xls >> list.snv
					ls TUMOR.$i.$tumor.SNV.filtered.xls >> list.filter.snv
				fi
				if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
				then
					ls TUMOR.$i.$tumor.INDEL.xls >> list.indel
					ls TUMOR.$i.$tumor.INDEL.filtered.xls >> list.filter.indel
				fi
			done
		done		
		if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
		then
			perl $script_path/union.snv.pl list.snv $output_dir/Reports/Paired.SNV.xls
			perl $script_path/union.snv.pl list.filter.snv $output_dir/Reports/Paired.SNV.filtered.xls
			rm list.snv list.filter.snv 
		fi
		if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
		then
			perl $script_path/union.indel.pl list.indel $output_dir/Reports/Paired.INDEL.xls
			perl $script_path/union.indel.pl list.filter.indel $output_dir/Reports/Paired.INDEL.filtered.xls
			rm list.indel list.filter.indel
		fi
	fi	
    echo `date`
fi	
