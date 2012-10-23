#!/bin/bash

if [ $# -le 2 ]
then
    echo -e "script to split the bam file per chromosome (assuimg the file name as <sample>.sorted.bam)\nUsage: ./split_bam_chr.sh </path/to/input directory> <sample name> </path/to/run_info.txt><SGE_TASK_ID(optional)>";
else
    set -x
    echo `date`
    input=$1
    sample=$2
    run_info=$3
	if [ $4 ]
	then
		SGE_TASK_ID=$4
	fi	
    bam=$input/$sample.sorted.bam
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    tool_info=$(cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2)
	if [ $multi_sample == "YES" ]
	then
		previous="reformat_pairBAM.sh"
	else
		previous="reformat_BAM.sh"
	fi	
    $samtools/samtools view -H $bam 1> $bam.$chr.header 2> $bam.$chr.sbc.fix.log
	if [[ `cat $bam.$chr.sbc.fix.log | wc -l` -gt 0 || `cat $bam.$chr.header |wc -l` -le 0 ]]
	then
		$script_path/email.sh $bam "bam is truncated or corrupt" $previous $run_info
		$script_path/wait.sh $bam.$chr.sbc.fix.log
	else
		rm $bam.$chr.sbc.fix.log
	fi	
	rm $bam.$chr.header	
	$samtools/samtools view -b $bam chr${chr} > $input/chr${chr}.cleaned.bam
    $samtools/samtools index $input/chr${chr}.cleaned.bam
    $samtools/samtools flagstat $input/chr${chr}.cleaned.bam > $input/chr${chr}.flagstat
    if [ ! -s $input/chr${chr}.cleaned.bam.bai ]
	then
		$script_path/errorlog.sh $input/chr${chr}.cleaned.bam split_bam_chr.sh ERROR "not created"  
		exit 1;
	else
		$script_path/filesize.sh SplittingBam $sample $input chr${chr}.cleaned.bam $run_info
	fi	
	echo `date`
fi	
	