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
    echo -e "Usage : SCRIPT to get Ontaget reads  \n<input sample realignment><chromsome><output Ontarget><sample><run info>";
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
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    gene_body=$( cat $tool_info | grep -w '^MATER_GENE_BODY' | cut -d '=' -f2 )
    multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    	
    if [ $tool == "whole_genome" ]
    then
        kit=$gene_body
    else
        kit=$CaptureKit
    fi    
    
    
	$samtools/samtools view -H $bam 2> $bam.fix.log
	if [ `cat $bam.fix.log | wc -l` -gt 0 ]
	then
		echo "$bam : bam file is truncated or corrupt"
		$script_path/email.sh $bam "bam is truncated or corrupt" $JOB_NAME $JOB_ID $run_info
		while [ -f $bam.fix.log ]
		do
			echo "waiting for the new and fixed bam file"
			sleep 2m
		done
	else
		rm $bam.fix.log
	fi	

    if [ $multi == "YES" ]
    then
        pair=$( cat $sample_info | grep -w "^$sample" | cut -d '=' -f2 | tr "\t" " ")
        for i in $pair
        do
            $bed/intersectBed -abam $input/$sample.$i.chr$chr.bam -b $kit | $samtools/samtools view -  | wc -l > $output/$sample.$i.chr$chr.bam.i.out  
        done
    else   
        bam=$input/chr$chr.cleaned.bam
		#intersect with the target kit
        $bed/intersectBed -abam $bam -b $kit | $samtools/samtools view - | wc -l > $output/$sample.chr$chr.bam.i.out
    fi
    echo `date`
fi	
	
### END OF SCRIPT	