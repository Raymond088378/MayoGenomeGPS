#!/bin/sh
##	INFO
#	To Intersect bam with OnTarget Kit by splitting the bam file into 200 files

######################################
#		$1		=	input folder (realignment sample folder)
#		$2		=	chromsome index
#		$3		=	Ontarget output folder
#		$4		=	sample name
#		$5		=	run info file
#########################################

if [ $# != 4 ];
then
    echo -e "Usage : SCRIPT to get Ontaget reads  \n<input sample realignment><chromsome><output Ontarget><sample><run info>";
else	
    set -x 
    echo `date`
    input=$1
    output=$2
    sample=$3
    run_info=$4
   # SGE_TASK_ID=2
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    bed=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    CaptureKit=$( cat $tool_info | grep -w '^CAPTUREKIT' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
    out=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    out_dir=$out/$PI/$tool/$run_num
    multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    
	PATH=$bed/:$PATH	
    if [ $tool == "whole_genome" ]
    then
        kit=$out_dir/bed_file.bed
    else
        kit=$CaptureKit
    fi    
    
    if [ ! -s $bam ]
    then
        echo " ERROR: OnTarget.BAM.sh $bam not found"
    fi
    
    if [ $multi == "YES" ]
    then
        pair=$( cat $sample_info | grep -w "$sample" | cut -d '=' -f2)
        for i in $pair
        do
            $bed/intersectBed -abam $input/$sample.$i.chr$chr.bam -b $kit | $samtools/samtools view -  | wc -l > $output/$i.chr$chr.bam.i.out
            if [ ! -s $output/$i.chr$chr.bam.i.out ]
            then
                echo "ERROR: $output/$i.chr$chr.bam.i.out is empty"
            fi    
        done
    else   
        bam=$input/chr$chr.cleaned.bam
		#intersect with the target kit
        $bed/intersectBed -abam $bam -b $kit | $samtools/samtools view -  | wc -l > $output/$sample.chr$chr.bam.i.out
        if [ ! -s $output/$sample.chr$chr.bam.i.out ]
        then
            echo "ERROR : $output/$sample.chr$chr.bam.i.out is empty"
        fi    
    fi
    echo `date`
fi	
	
### END OF SCRIPT	