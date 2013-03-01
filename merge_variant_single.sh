#!/bin/bash

########################################################
###### 	Merges variants from vcf files by chromosome
######		Program:			merge_variant_group.sh
######		Date:				12/13/2011
######		Summary:			Using PICARD to sort and mark duplicates in bam 
######		Input files:		$1	=	/path/to/input directory
######					$2	=	sample name
######					$3	=	/path/to/run_info.txt
########################################################


### inputs
###	$input/$sample.variants.chr$chr.raw.vcf
###
### outputs
### $out/$sample.variants.raw.vcf
### $out/$sample.variants.filter.vcf


if [ $# != 4 ];
then
    echo -e "script to merge the varaiants and then apply filters to the vcf file\nUsage: ./merge_variant_single.sh </path/to/input directory> <sample name> </path/to/output folder> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    input=$1
    sample=$2
    out=$3
    run_info=$4
    
	### reading run_info.txt and assign variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" )
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    filter_variants=$( cat $tool_info | grep -w '^VARIANT_FILTER' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    depth=$( cat $tool_info | grep -w '^T_DEPTH_FILTER' | cut -d '=' -f2 )
    perllib=$( cat $tool_info | grep -w '^PERLLIB' | cut -d '=' -f2)
    export PERL5LIB=$perllib:$PERL5LIB
    export PATH=$PERL5LIB:$PATH
    export PATH=$java:$PATH
	
    
    ### loop over all chromosomes to create list of files to merge
    inputargs=""
    for i in $chrs
    do
		inputfile=$input/$sample/$sample.variants.chr$i.raw.vcf 
		if [ ! -s $inputfile ]
		then	
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" variants.sh $run_info
			$script_path/wait.sh $inputfile.fix.log
		fi
		inputargs=$inputargs"$inputfile "
    done
	
	### concatenate all vcf files in file list together
    $script_path/concatvcf.sh "$inputargs" $out/$sample.variants.raw.vcf $run_info no
	
    ### filter the variant calls
    if [ $filter_variants == "YES" ]
    then
        $script_path/filter_variant_vqsr.sh $out/$sample.variants.raw.vcf $out/$sample.variants.filter.vcf BOTH $run_info
    else
        cp $out/$sample.variants.raw.vcf $out/$sample.variants.filter.vcf
    fi
    
    ### check whether vqsr produced output / file was copied correctly
	if [ ! -s $out/$sample.variants.filter.vcf ]
    then
    	touch $out/$sample.variants.filter.vcf.fix.log
    	$script_path/email.sh $out/$sample.variants.filter.vcf "vqsr failed" filter_variant_vqsr.sh $run_info
		$script_path/wait.sh $out/$sample.variants.filter.vcf.fix.log
	fi	
	
	### Filter the variants using total depth 
	### use GATK variant filter to filter using DP
	if [[ $tool == "exome" && $filter_variants == "YES" ]]
	then
		$script_path/filtervcf.sh $out/$sample.variants.filter.vcf $run_info 
	fi
	
    if [ ! -s $out/$sample.variants.filter.vcf ]
    then
        $script_path/errorlog.sh $out/$sample.variants.filter.vcf merge_variant_single.sh ERROR "failed to create"
        exit 1;
    else
        for i in $chrs
        do
            cat $out/$sample.variants.filter.vcf | awk -v num=chr${i} '$0 ~ /^#/ || $1 == num' > $input/$sample/$sample.variants.chr$i.filter.vcf 
        done
    fi  
	echo `date`	
fi  
### at this point, we should have $sample.variants.filter.vcf in whatever output folder given
