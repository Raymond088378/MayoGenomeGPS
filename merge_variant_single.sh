#!/bin/sh

########################################################
###### 	Merges variants from vcf files by chromosome

######		Program:			merge_variant_group.sh
######		Date:				12/13/2011
######		Summary:			Using PICARD to sort and mark duplicates in bam 
######		Input files:		$1	=	/path/to/input directory
######					$2	=	sample name
######					$3	=	/path/to/run_info.txt
########################################################

if [ $# != 4 ];
then
    echo "Usage: </path/to/input directory> <sample name> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    input=$1
    sample=$2
    out=$3
    run_info=$4
    
########################################################	
######		Reading run_info.txt and assigning to variables
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
    all_sites=$( cat $tool_info | grep -w '^EMIT_ALL_SITES' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
    javahome=$( cat $tool_info | grep -w '^JAVA_HOME' | cut -d '=' -f2 )
    threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    depth=$( cat $tool_info | grep -w '^T_DEPTH_FILTER' | cut -d '=' -f2 )
    perllib=$( cat $tool_info | grep -w '^PERLLIB' | cut -d '=' -f2)
    export PERL5LIB=$perllib:$PERL5LIB
    export PATH=$PERL5LIB:$PATH
    export JAVA_HOME=$javahome
    export PATH=$javahome/bin:$PATH
	
    
    inputargs=""
    inputargs_multi=""
    for i in $chrs
    do
		inputfile=$input/$sample/$sample.variants.chr$i.raw.vcf 
		multi=$input/$sample/$sample.variants.chr$i.raw.multi.vcf 
		if [ ! -s $inputfile ]
		then	
			$script_path/errorlog.sh $inputfile merge_variant_single.sh ERROR "not exist"
			exit 1;
		else
			inputargs=$inputargs"$inputfile "
			inputargs_multi=$inputargs_multi"$multi "
		fi
    done
	
    $script_path/concatvcf.sh "$inputargs" $out/$sample.variants.raw.vcf $run_info no
	$script_path/concatvcf.sh "$inputargs_multi" $out/$sample.variants.raw.multi.vcf $run_info yes
	
	
    ### filter the variant calls
    if [ $filter_variants == "YES" ]
    then
        $script_path/filter_variant_vqsr.sh $out/$sample.variants.raw.vcf $out/$sample.variants.filter.vcf BOTH $run_info
    else
        cp $out/$sample.variants.raw.vcf $out/$sample.variants.filter.vcf
    fi
	
	### Filter the variants using total depth 
	### use GATK variant filter to filter using DP
	if [[ $tool == "exome" && $filter_variants == "YES" ]]
	then
		$java/java -Xmx1g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
		-R $ref \
		-et NO_ET \
		-K $gatk/Hossain.Asif_mayo.edu.key \
		-l INFO \
		-T VariantFiltration \
		-V $out/$sample.variants.filter.vcf \
		-o $out/$sample.variants.filter.vcf.tmp \
		--filterExpression "DP < $depth" --filterName DPFilter 
		mv $out/$sample.variants.filter.vcf.tmp $out/$sample.variants.filter.vcf
		mv $out/$sample.variants.filter.vcf.tmp.idx $out/$sample.variants.filter.vcf.idx
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
