#!/bin/bash

if [ $# -le 2 ];
then
    echo -e "wrapper script to run the alignment using BWA\nUsage ./align_bwa.sh <sample name> </path/to/output_dir> </path/to/run_info.txt> <SGE TASK ID (optional)";
else	
    set -x 
    echo `date`
    sample=$1
    output_dir=$2
    run_info=$3
    if [ $4 ]
    then
        SGE_TASK_ID=$4
    fi    
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	genome_bwa=$( cat $tool_info | grep -w '^BWA_REF' | cut -d '=' -f2)
	bwa=$( cat $tool_info | grep -w '^BWA' | cut -d '=' -f2)
	center=$( cat $tool_info | grep -w '^CENTER' | cut -d '=' -f2 )
	platform=$( cat $tool_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
	GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	paired=$( cat $run_info | grep -w '^PAIRED' | cut -d '=' -f2)
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
	
	output_dir_sample=$output_dir/alignment/$sample
    fastq=$output_dir/fastq
	
	if [ $paired == 1 ]
	then
		let fidx=($SGE_TASK_ID*2)-1 
		let sidx=($SGE_TASK_ID*2)
		R1=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
        R2=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $sidx | tail -n 1`
		$script_path/sampe.sh $sample $output_dir_sample/$sample.$SGE_TASK_ID.R1.sai $output_dir_sample/$sample.$SGE_TASK_ID.R2.sai \
			$fastq/$R1 $fastq/$R2 $output_dir_sample/$sample.$SGE_TASK_ID.bam $tool_info
	else
		let fidx=($SGE_TASK_ID*2)-1 
		R1=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
		$script_path/samse.sh $sample $output_dir_sample/$sample.$SGE_TASK_ID.R1.sai $fastq/$R1\
			$output_dir_sample/$sample.$SGE_TASK_ID.bam $tool_info	
	fi
		$script_path/filesize.sh alignment.out $sample $output_dir_sample $sample.$SGE_TASK_ID.bam $run_info
 
	$samtools/samtools view -H $output_dir_sample/$sample.$SGE_TASK_ID.bam 1>$output_dir_sample/$sample.$SGE_TASK_ID.bam.header 2> $output_dir_sample/$sample.$SGE_TASK_ID.bam.log
	if [[ `cat $output_dir_sample/$sample.$SGE_TASK_ID.bam.log | wc -l` -gt 0 || `cat $output_dir_sample/$sample.$SGE_TASK_ID.bam.header | wc -l` -le 0 ]]	
	then
		$script_path/errorlog.sh $output_dir_sample/$sample.$SGE_TASK_ID.bam align_bwa.sh ERROR "truncated or corrupt"
		exit 1;
	else
		rm $output_dir_sample/$sample.$SGE_TASK_ID.bam.log
                if [ $paired == 0 ]
        then
            rm $fastq/$R1
            rm $output_dir_sample/$sample.$SGE_TASK_ID.R1.sai
        else
            rm $fastq/$R1 $fastq/$R2
            rm $output_dir_sample/$sample.$SGE_TASK_ID.R1.sai $output_dir_sample/$sample.$SGE_TASK_ID.R2.sai
        fi    
	fi	
	rm $output_dir_sample/$sample.$SGE_TASK_ID.bam.header
    echo `date`
fi	
