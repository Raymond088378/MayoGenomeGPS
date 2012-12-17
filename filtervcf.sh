#!/bin/bash

if [ $# != 2 ]
then
	echo -e "script to filter the vcf using the expression\nUsage: ./filtervcf.sh <inputvcf><runinfo file>"
else
	set -x
	echo `date`
	inputvcf=$1
	run_info=$2
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)	
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	depth=$( cat $tool_info | grep -w '^T_DEPTH_FILTER' | cut -d '=' -f2 )
	
	mem=$( cat $memory_info | grep -w '^VariantFiltration_JVM' | cut -d '=' -f2)
	$java/java $mem -jar $gatk/GenomeAnalysisTK.jar \
	-R $ref \
	-et NO_ET \
	-K $gatk/Hossain.Asif_mayo.edu.key \
	-l INFO \
	-T VariantFiltration \
	-V $inputvcf \
	-o $inputvcf.tmp --filterExpression "DP < $depth" --filterName DPFilter 
	
	if [ -s $inputvcf.tmp.idx ]
	then
		mv $inputvcf.tmp $inputvcf
		mv $inputvcf.temp.idx $inputvcf.idx
	fi
	echo `date`	
fi