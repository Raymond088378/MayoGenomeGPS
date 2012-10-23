#!/bin/bash
##	INFO
#	To Intersect bam with OnTarget Kit by splitting the bam file into 200 files

######################################
#		$1		=	input folder (realignment sample folder)
#		$2		=	chromsome index
#		$3		=	Ontarget output folder
#		$4		=	sample name
#		$5		=	run info file
#########################################

if [ $# -le 3 ];
then
    echo -e "Usage : SCRIPT to get Ontaget reads\nUsage: ./OnTarget_BAM.sh </path/to/input sample realignment></path/to/output Ontarget><sample></path/to/run info><SGE_TASK_ID(optional)>";
else	
    set -x 
    echo `date`
    input=$1
    output=$2
    sample=$3
    run_info=$4
	if [ $5 ]
    then
		SGE_TASK_ID=$5
    fi
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    bed=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    CaptureKit=$( cat $tool_info | grep -w '^CAPTUREKIT' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    gene_body=$( cat $tool_info | grep -w '^MATER_GENE_BODY' | cut -d '=' -f2 )
    multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]")	
	
    if [ $tool == "whole_genome" ]
    then
        kit=$gene_body
    else
        kit=$CaptureKit
    fi    
	
    if [ $multi != "YES" ]
	then
		bam=$input/chr$chr.cleaned.bam
		$samtools/samtools view -H $bam 1>$bam.OnTarget_BAM.header 2> $bam.fix.OnTarget_BAM.log
		if [[ `cat $bam.fix.OnTarget_BAM.log | wc -l` -gt 0 || `cat $bam.OnTarget_BAM.header | wc -l` -le 0 ]]
		then
			$script_path/email.sh $bam "bam is truncated or corrupt" realign_recal.sh $run_info
			$script_path/wait.sh $bam.fix.OnTarget_BAM.log
		else
			rm $bam.fix.OnTarget_BAM.log
		fi
		rm $bam.OnTarget_BAM.header	
	fi
	
    if [ $multi == "YES" ]
    then
        pair=$( cat $sample_info | grep -w "^$sample" | cut -d '=' -f2 | tr "\t" " ")
        for i in $pair
        do
            $bed/intersectBed -abam $input/$sample.$i.chr$chr.bam -b $kit | $samtools/samtools view -  | wc -l > $output/$sample.$i.chr$chr.bam.i.out
        if [ ! -s $output/$sample.$i.chr$chr.bam.i.out ]
	then
		$script_path/errorlog.sh $output/$sample.$i.chr$chr.bam.i.out OnTarget_BAM.sh ERROR "failed to create"
		exit 1;
	fi	
        done 
    else   
        $bed/intersectBed -abam $bam -b $kit | $samtools/samtools view - | wc -l > $output/$sample.chr$chr.bam.i.out
        if [ ! -s $output/$sample.chr$chr.bam.i.out ]
	then
		$script_path/errorlog.sh $output/$sample.chr$chr.bam.i.out OnTarget_BAM.sh ERROR "failed to create"
		exit 1;
	fi	
    fi
    
	echo `date`
fi	
	
### END OF SCRIPT	