#!/bin/sh

########################################################
###### 	Merges variants from vcf files by chromosome

######		Program:			merge_variant_group.sh
######		Date:				12/13/2011
######		Summary:			Using PICARD to sort and mark duplicates in bam 
######		Input files:		$1	=	/path/to/input directory
######					$2	=	group name
######					$3	=	/path/to/run_info.txt
########################################################

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
	
########################################################	
######		Reading run_info.txt and assigning to variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" " " )
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
    output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    filter_variants=$( cat $tool_info | grep -w '^VARIANT_FILTER' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    
########################################################	

    inputargs=""
    for i in $chrs
    do
        inputfile=$input/$group/MergeAllSamples.chr$i.raw.vcf 
        if [ ! -s $inputfile ]
        then		
            echo "ERROR :merge_variant_group_chr. Somatic variant file for group $group, chromosome $i: $inputfile does not exist "
            exit 1
        else
            inputargs="-V $inputfile "$inputargs
        fi
    done

    $script_path/combinevcf.sh "$inputargs" $out/$group.somatic.variants.raw.vcf $run_info no

    if [ $filter_variants == "YES" ]
    then
        $script_path/filter_variant_vqsr.sh $out/$group.somatic.variants.raw.vcf $out/$group.somatic.variants.filter.vcf BOTH $run_info
    else
        cp $out/$group.somatic.variants.raw.vcf $out/$group.somatic.variants.filter.vcf
    fi
    
    if [ ! -s $out/$group.somatic.variants.filter.vcf ]
    then
	echo "ERROR: $out/$group.somatic.variants.filter.vcf not exist"
    else
	for chr in $chrs
	do
	    cat $out/$group.somatic.variants.filter.vcf | awk -v num=chr${chr} '$0 ~ /^#/ || $1 == num' > $input/$group/$group.somatic.variants.chr$chr.filter.vcf 
	done
    fi
			
    inputargs=""
    for i in $chrs
    do
        inputfile=$input/$group/variants.chr$i.raw.vcf 
        if [ ! -s $inputfile ]
        then		
            echo "ERROR :merge_variant_group_chr. Variant file for group $group, chromosome $i: $inputfile does not exist "
            exit 1
        else
            inputargs="-V $inputfile "$inputargs
        fi
    done

    $script_path/combinevcf.sh "$inputargs" $out/$group.variants.raw.vcf $run_info no
	
    if [ $filter_variants == "YES" ]
    then
        $script_path/filter_variant_vqsr.sh $out/$group.variants.raw.vcf $out/$group.variants.filter.vcf BOTH $run_info
    else
        cp $out/$group.variants.raw.vcf $out/$group.variants.filter.vcf
    fi    
    
    if [ ! -s $out/$group.variants.filter.vcf ]
    then
        echo "ERROR: $out/$group.variants.filter.vcf failed to generate the filterd vcf for $sample"
        exit 1
    else
	for chr in $chrs
        do
            cat $out/$group.variants.filter.vcf | awk -v num=chr${chr} '$0 ~ /^#/ || $1 == num' > $input/$group/$group.variants.chr$chr.filter.vcf 
        done
    fi
    echo `date`
fi