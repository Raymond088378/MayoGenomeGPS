#!/bin/bash

########################################################
###### 	ALIGNMENT SCRIPT FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			bwa.splitted.PR.sh
######		Date:				07/25/2011
######		Summary:			Alignment done using BWA on PE FASTQ files and realignment of the
######                                          soft masked reads with novoalign
######							$1	=	/path/to/output directory
######							$2	=	sample name
######							$3	=	/path/to/run_info.txt
######		Output files:		BAM files
######          now checking the quality scores for fastq's and using the specific paramter for illumina and sanger quality in nova align
########################################################

if [ $# -le 2 ]
then
    echo -e "wrapper script to run the alignment using NOVO ALIGN \
		\nUsage: ./align_novo.sh <sample name> </path/to/output_dir> </path/to/run_info.txt> <SGE TASK ID (optional)>";
	exit 1;	
fi
	
    set -x 
    echo `date`
    sample=$1
    output_dir=$2
    run_info=$3
	if [ $4 ]
	then
		SGE_TASK_ID=$4
	fi	
    
########################################################	
######	Reading run_info.txt and assigning to variables
    seq_file=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    center=$( cat $tool_info | grep -w '^CENTER' | cut -d '=' -f2 )
    platform=$( cat $tool_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    paired=$( cat $run_info | grep -w '^PAIRED' | cut -d '=' -f2)
	FOLDER_FASTQC=$( cat $run_info | grep -w '^FOLDER_FASTQC' | cut -d '=' -f2 )
########################################################	
######		Check FASTQ for Illumina or Sanger quality scrore
    
    output_dir_sample=$output_dir/alignment/$sample
    fastq=$output_dir/fastq
    fastqc=$output_dir/fastqc
    
	$script_path/dashboard.sh $sample $run_info Alignment started $SGE_TASK_ID
   
    if [[ $paired == 1 || $paired == "YES" ]]
    then
        let fidx=($SGE_TASK_ID*2)-1 
        let sidx=$SGE_TASK_ID*2
        R1=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
        R2=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $sidx | tail -n 1`
        for i in $R1 $R2
        do
			if [ ! -s $seq_file/$i ]
			then
				$script_path/errorlog.sh align_novo.sh $seq_file/$i ERROR "not found"
				exit 1;
			fi	
			ln -s $seq_file/$i $fastq/$i
			$script_path/fastqc.sh $fastqc $fastq/$i $tool_info $FOLDER_FASTQC
        done
        sequences="$R1:$R2"
    elif [[ $paired == 0 || $paired == "NO" ]]
    then
        let fidx=$SGE_TASK_ID
        R1=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
        if [ ! -s $seq_file/$R1 ]
		then
			$script_path/errorlog.sh align_novo.sh $seq_file/$R1 ERROR "not found"
			exit 1;
		fi	
		ln -s $seq_file/$R1 $fastq/$R1
		$script_path/fastqc.sh $fastqc $fastq/$R1 $tool_info $FOLDER_FASTQC
		sequences="$R1"
    fi  

########################################################	
######		Run novoalign for PE or SR	
	$script_path/novoalign.sh $fastq "$sequences" $sample $output_dir_sample/$sample.$SGE_TASK_ID.bam $tool_info 
    
    echo `date`	
