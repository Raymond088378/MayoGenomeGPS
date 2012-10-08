#!/bin/bash

if [ $# != 6 ]
then
	echo -e "script to run snvmix\nUsage: <sample> <bam file><output vcf> <mode><bed file><run info> "
else
    set -x
    echo `date`
    sample=$1
    bam=`echo $2 |  sed -e '/-I/s///g' | sed -e "s/^ *//"`
    output=$3
    mode=`echo $4 | tr "[A-Z]" "[a-z]"`
    kit=$5
    run_info=$6
	
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    snvmix=$( cat $tool_info | grep -w '^SNVmix' | cut -d '=' -f2)
    only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    filter=$( cat $tool_info | grep -w '^SNVMIX2_Filter' | cut -d '=' -f2)
    TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	command_line_params=$( cat $tool_info | grep -w '^SNVMIX2_params' | cut -d '=' -f2 )
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    temp=`echo $output | sed -e '/.vcf/s///g'`
	
	$samtools/samtools view -H $bam 1>$bam.snvmix.header 2>$bam.snvmix.fix.log
	if [ `cat $bam.snvmix.fix.log | wc -l` -gt 0 ]
	then
		$script_path/email.sh $bam "bam is truncated or corrupt" $run_info
		while [ -f $bam.snvmix.fix.log ]
		do
			echo "waiting for the $bam to be fixed"
			sleep 2m
		done
	else
		rm $bam.snvmix.fix.log
	fi		
	rm $bam.snvmix.header
	mkfifo $bam.pileup 
	$samtools/samtools mpileup -A -s -f $ref $bam > $bam.pileup &
	pileup=$bam.pileup

    if [ $mode == "all" ]
    then	
        $snvmix/SNVMix2 -i $pileup -f -m $snvmix/Mu_pi.txt -o $temp "$command_line_params"
    else
        $snvmix/SNVMix2 -i $pileup -m $snvmix/Mu_pi.txt -o $temp "$command_line_params"
    fi

    if [[ $only_ontarget == "YES" && $tool == "exome" ]]
    then
        $script_path/snvmix_to_vcf.pl -i $temp -o $output -s $sample $filter
        if [ -s $temp_kit ]
		then
            $bedtools/intersectBed -a $output -b $TargetKit -wa -header > $output.i
        else
            cp $output $output.i	
        fi	
        mv $output.i $output
    else
        $script_path/snvmix_to_vcf.pl -i $temp -o $output -s $sample
    fi
    if [ $mode != "all" ]
    then
        cat $output | grep -v "0\/0" > $output.temp
        mv $output.temp $output
    fi    
	
	cat $output | awk '$0 ~ /^#/ || $5 ~ /,/' > $output.multi.vcf
	cat $output | awk '$0 ~ /^#/ || $5 !~ /,/' > $output.temp
	mv $output.temp $output	
	
    if [ -s $output ]
    then
        rm $temp
        rm $pileup
    else
		$script_path/errorlog.sh $output snvmix2.sh ERROR "failed to create"	
        exit 1;
    fi	
    echo `date`
fi	    