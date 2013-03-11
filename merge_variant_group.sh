#!/bin/bash

########################################################
###### 	Merges variants from vcf files by chromosome
######		Program:			merge_variant_group.sh
######		Date:				12/13/2011
######		Summary:			Using PICARD to sort and mark duplicates in bam 
######		Input files:		$1	=	/path/to/input directory
######					$2	=	group name
######					$3	=	/path/to/run_info.txt
########################################################

### SOMATIC INPUTS
### inputs
###	$input/MergeAllSamples.chr$chr.vcf
###
### outputs
### $out/$group.somatic.variants.vcf
### $out/$group.somatic.variants.final.vcf
###
### INDIVIDUAL SAMPLES
###
### inputs
### $input/$group/variants.chr$chr.vcf
### outputs
### $out/$group.variants.vcf
### $out/$group.variants.final.vcf

if [ $# != 4 ];
then
    echo "Usage: </path/to/input directory> <group name> </path/to/output directory> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    input=$1
    group=$2
    out=$3
    run_info=$4
	
	
	###	read run_info.txt and assign variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" " " )
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
    output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    filter_variants=$( cat $tool_info | grep -w '^VARIANT_FILTER' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
	somatic_filter_variants=$( cat $tool_info | grep -w '^SOMATIC_VARIANT_FILTER' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	perllib=$( cat $tool_info | grep -w '^PERLLIB' | cut -d '=' -f2)
	depth=$( cat $tool_info | grep -w '^T_DEPTH_FILTER' | cut -d '=' -f2 )
	export PERL5LIB=$perllib:$PERL5LIB
	export PATH=$PERL5LIB:$PATH
	export PATH=$java:$PATH
    somatic_calling=$( cat $tool_info | grep -w '^SOMATIC_CALLING' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )    
	########################################################	

    if [ $somatic_calling == "YES" ]
	then
		inputargs=""
		inputargs_multi=""
		for i in $chrs
		do
			inputfile=$input/$group/MergeAllSamples.chr$i.vcf 
			if [ ! -s $inputfile ]
			then		
				touch $inputfile.fix.log
				$script_path/email.sh $inputfile "not exist" variants.sh $run_info
				$script_path/wait.sh $inputfile.fix.log
			fi	
			inputargs=$inputargs"$inputfile "
		done
		
		$script_path/concatvcf.sh "$inputargs" $out/$group.somatic.variants.vcf $run_info no
		
		if [ $somatic_filter_variants == "YES" ]
		then
			$script_path/filter_variant_vqsr.sh $out/$group.somatic.variants.vcf $out/$group.somatic.variants.final.vcf BOTH $run_info somatic
		else
			cp $out/$group.somatic.variants.vcf $out/$group.somatic.variants.final.vcf
		fi
		if [ ! -s $out/$group.somatic.variants.final.vcf ]
		then
			touch $out/$group.somatic.variants.final.vcf.fix.log
			$script_path/email.sh $out/$group.somatic.variants.final.vcf "vqsr failed" filter_variant_vqsr.sh $run_info
			$script_path/wait.sh $out/$group.somatic.variants.final.vcf.fix.log
		fi	
			cat $out/$group.somatic.variants.final.vcf | sed -e 's/-1/\./g' > $out/$group.somatic.variants.final.vcf.tmp.vcf
			mv $out/$group.somatic.variants.final.vcf.tmp.vcf $out/$group.somatic.variants.final.vcf
		if [ ! -s $out/$group.somatic.variants.final.vcf ]
		then
			$script_path/errorlog.sh $out/$group.somatic.variants.final.vcf merge_variant_greoup.sh ERROR "does not exist"
			exit 1;
		else
		for chr in $chrs
		do
			cat $out/$group.somatic.variants.final.vcf | awk -v num=chr${chr} '$0 ~ /^#/ || $1 == num' > $input/$group/$group.somatic.variants.chr$chr.final.vcf 
		done
		fi
	fi
	
	### multi sample calling	
    inputargs=""
    for i in $chrs
    do
        inputfile=$input/$group/variants.chr$i.vcf 
        if [ ! -s $inputfile ]
        then		
            touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" variants.sh $run_info
			$script_path/wait.sh $inputfile.fix.log
		fi
        inputargs=$inputargs"$inputfile "
    done

	$script_path/concatvcf.sh "$inputargs" $out/$group.variants.vcf $run_info no
	
    if [ $filter_variants == "YES" ]
    then
        $script_path/filter_variant_vqsr.sh $out/$group.variants.vcf $out/$group.variants.final.vcf BOTH $run_info
    else
        cp $out/$group.variants.vcf $out/$group.variants.final.vcf
    fi    
	
    if [ ! -s $out/$group.variants.final.vcf ]
    then
    	touch $out/$group.variants.final.vcf.fix.log
    	$script_path/email.sh $out/$group.variants.final.vcf "vqsr failed" filter_variant_vqsr.sh $run_info
		$script_path/wait.sh $out/$group.variants.final.vcf.fix.log
	fi	
		
	### Filter the variants using total depth 
	### use GATK variant filter to filter using DP
	if [ $tool == "exome" ]
	then
		if [ $depth -gt 0 ]
			then
				$script_path/filtervcf.sh $out/$group.variants.final.vcf $run_info 
			fi
        fi
        cat $out/$group.variants.final.vcf | sed -e 's/-1/\./g' > $out/$group.variants.final.vcf.tmp.vcf
        mv $out/$group.variants.final.vcf.tmp.vcf $out/$group.variants.final.vcf  
	
    if [ ! -s $out/$group.variants.final.vcf ]
    then
        $script_path/errorlog.sh $out/$group.variants.final.vcf merge_variant_group.sh ERROR "does not exist"
        exit 1
    else
		for chr in $chrs
		do
			cat $out/$group.variants.final.vcf | awk -v num=chr${chr} '$0 ~ /^#/ || $1 == num' > $input/$group/$group.variants.chr$chr.final.vcf 
		done
    fi
	echo `date`
fi