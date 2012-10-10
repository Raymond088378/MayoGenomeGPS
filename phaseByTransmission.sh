#!/bin/bash

if [ $# != 2 ]
then
    echo -e "script to run phaseByTransmission\nUsage: ./phaseByTransmission.sh <vcf input> <vcf output> <run info file>"
else
    set -x
    echo `date`
    input_vcf=$1
	output_vcf=$2
    run_info=$3

    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    ped=$( cat $tool_info | grep -w '^PEDIGREE' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	mem=$( cat $memory_info | grep -w '^PhaseByTransmission_JVM' | cut -d '=' -f2)
	export PATH=$java:$PATH
    check=`[ -s $output_vcf.idx ] && echo "1" || echo "0"`
    while [ $check -eq 0 ]
    do
		$java/java $mem -jar $gatk/GenomeAnalysisTK.jar \
		-R $ref \
		-T PhaseByTransmission \
		-v $input_vcf \
		--ped $ped \
		--out $output_vcf
		check=`[ -s $output_vcf.idx ] && echo "1" || echo "0"`
    done

    if [ ! -s $output_vcf ]
    then
		$script_path/errorlog.sh $output_vcf phaseByTransmission.sh ERROR "failed to create"
        exit 1;
	else
		echo "Replacing PhaseByTransmisson vcf with the one from UnifiedGenoType"
		mv $output_vcf $input_vcf
    fi
    echo `date`
fi
