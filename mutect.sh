#!/bin/sh
if [ $# != 8 ]
then
	echo "Usage: <normal bam> <tumor bam > <output dir> <chromosome> <tumor sample name> <normal sample name> <output vcf file name> <run info>"
else
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
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	javahome=$( cat $tool_info | grep -w '^JAVA_HOME' | cut -d '=' -f2 )
	command_line_params=$( cat $tool_info | grep -w '^MUTECT_params' | cut -d '=' -f2 )
	
	export JAVA_HOME=$javahome
	export PATH=$javahome/bin:$PATH
	
        
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
    
	 #param="-L chr$chr"
    check=0
    if [ ! -d $output/temp ]
	then
		mkdir $output/temp
	fi
	
	while [[ $check -eq 0 ]]
    do
        if [ `grep -l $output/$output_file *.log` ]
        then
            rm `grep -l $output/$output_file *.log` 
        fi
        $java/java -XX:MaxPermSize=128M -Xmx6g -Xms512m -jar $mutect/muTect-1.0.27783.jar \
        -T MuTect \
        --reference_sequence $ref \
        $param -nt $threads \
        --input_file:tumor $tumor_bam \
        --input_file:normal $normal_bam \
        -B:dbsnp,VCF $dbSNP \
        -et NO_ET \
        --out $output/$output_file $command_line_params
        len=`cat *.log | grep $output/$output_file | wc -l`
        check=` [ $len -gt 0 ] && echo "0" || echo "1"`
        if [ $check -eq 0 ]
        then
            rm core.*
        fi
    done    
    
	perl $script_path/mutect2vcf.pl -i $output/$output_file -o $output/$output_file.temp -ns $normal_sample -ts $tumor_sample
	mv $output/$output_file.temp  $output/$output_file
	cat $output/$output_file | awk '$0 ~ /^#/ || $5 ~ /,/' > $output/$output_file.multi.vcf
	cat $output/$output_file | awk '$0 ~ /^#/ || $5 !~ /,/' > $output/$output_file.tmp
	mv $output/$output_file.tmp  $output/$output_file
	
	if [ ! -s  $output/$output_file ]
	then
		$script_path/errorlog.sh $output/$output_file mutect.sh ERROR "failed to create"
		exit 1;
	fi
	echo `date`
fi    
