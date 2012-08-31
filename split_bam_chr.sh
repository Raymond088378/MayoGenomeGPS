#!/bin/bash
if [ $# -le 2 ];
then
    echo -e "Usage: wrapper to merge bam files and validate the bam for downstream analysis \n merge_align.bam.sh </path/to/input directory> <name of BAM to sort> <sample name> </path/to/run_info.txt>";
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
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	
    $samtools/samtools view -H $bam 2> $bam.fix.log
	if [ `cat $bam.fix.log | wc -l` -gt 0 ]
	then
		$script_path/email.sh $bam "bam is truncated or corrupt" $JOB_NAME $JOB_ID $run_info
		while [ -f $bam.fix.log ]
		do
			echo "waiting for the $bam to be fixed"
			sleep 2m
		done
	else
		rm $bam.fix.log
	fi		
	$samtools/samtools view -b $bam chr${chr} > $input/chr${chr}.cleaned.bam
    $samtools/samtools index $input/chr${chr}.cleaned.bam
    $samtools/samtools flagstat $input/chr${chr}.cleaned.bam > $input/chr${chr}.flagstat
    echo `date`
fi	
	