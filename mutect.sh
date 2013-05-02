#!/bin/bash
if [ $# != 8 ]
then
	echo -e "script to run somatic caller called mutect\nUsage: ./mutect.sh <normal bam> <tumor bam > <output dir> <chromosome> <tumor sample name> <normal sample name> <output vcf file name> <run info>"
	exit 1;
fi	
	set -x
	echo `date`
	normal_bam=$1
	tumor_bam=$2
	output=$3
	chr=$4
	tumor_sample=$5
	normal_sample=$6
	output_file=$7
	run_info=$8
    
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	mutect=$( cat $tool_info | grep -w '^MUTECT' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
	threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2)
	TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	command_line_params=$( cat $tool_info | grep -w '^MUTECT_params' | cut -d '=' -f2 )
	memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	mem=$( cat $memory_info | grep -w '^MuTecT_JVM=' | sed -e '/MuTecT_JVM=/s///g')
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2) 
	export PATH=$java:$PATH  
    if [ $only_ontarget == "YES" ]
    then
        len=`cat $output/chr$chr.target.bed |wc -l`
        if [ $len -gt 0 ]
        then
            param="-L $output/chr$chr.target.bed"
        else
            param="-L chr$chr"
        fi    
    else
        param="-L chr$chr"
    fi
    
    let check=0
    if [ ! -d $output/temp ]
	then
		mkdir -p $output/temp
	fi
	gatk_params="-R $ref -et NO_ET-k $gatk/Hossain.Asif_mayo.edu.key "
	while [[ $check -eq 0 ]]
    do
		$java/java $mem -jar $mutect \
        -T MuTect \
        -R $ref \
        -nt $threads \
        -I:tumor $tumor_bam --tumor_sample_name $tumor_sample \
        -I:normal $normal_bam --normal_sample_name $normal_sample \
        -B:dbsnp,VCF $dbSNP \
        -et NO_ET \
        --out $output/$output_file.txt -vcf $output/$output_file $gatk_params $param $command_line_params
        rm $output/$output_file.txt
       
        mutectpid=$!
        
        if [[  `find . -name '*.log'` ]]
		then
			len=`cat *.log | grep $output/$output_file | wc -l`
			check=`[ $len -gt 0 ] && echo "0" || echo "1"`
			if [ $check -eq 0 ]
			then
				if [ `grep -l $output/$output_file *.log` ]
				then
					rm `grep -l $output/$output_file *.log`
					rm core.*
				fi
			fi
		else
            let check=1    
        fi
        
        if [[ $check -eq 0 ]]
        then
        	kill -9 $mutectpid
    	fi
    	
    done    
    cat $output/$output_file | sed -e 's/SS:/MUTX_SS:/g' | sed -e 's/SS=/MUTX_SS=/g' > $output/$output_file.tmp
    mv $output/$output_file.tmp $output/$output_file
	if [ ! -s  $output/$output_file ]
	then
		$script_path/errorlog.sh $output/$output_file mutect.sh ERROR "failed to create"
		exit 1;
	fi
	echo `date`
    
