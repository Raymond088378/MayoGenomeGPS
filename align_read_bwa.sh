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

if [ $# -le 3 ]
then
    echo -e "wrapper script to run the alignment using NOVO ALIGN\nUsage:align_read_bwa.sh <sample name> </path/to/output_dir> </path/to/run_info.txt> <SGE TASK ID (optional)>";
else	
    set -x 
    echo `date`
    sample=$1
    output_dir=$2
    read=$3
    run_info=$4
	if [ $5 ]
	then
		SGE_TASK_ID=$5
	fi	
    
########################################################	
######	Reading run_info.txt and assigning to variables
    seq_file=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    genome_bwa=$( cat $tool_info | grep -w '^BWA_REF' | cut -d '=' -f2)
    bwa=$( cat $tool_info | grep -w '^BWA' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    parameters=$( cat $tool_info | grep -w '^BWA_params' | cut -d '=' -f2 )
########################################################	
######		Check FASTQ for Illumina or Sanger quality scrore
    
    output_dir_sample=$output_dir/alignment/$sample
    fastq=$output_dir/fastq
    fastqc=$output_dir/fastqc
    
    if [ $read -eq 1 ]
    then
        $script_path/dashboard.sh $sample $run_info Alignment started $SGE_TASK_ID
    fi
    
    if [ $read -eq 1 ]
    then
		let fidx=($SGE_TASK_ID*2)-1 
    else
		let fidx=($SGE_TASK_ID*2)
    fi	
	
    R1=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
    if [ ! -s $seq_file/$fastq ]
	then
		$script_path/errorlog.sh align_read_bwa.sh $seq_file/$fastq ERROR "not found"
		exit 1;
	fi	
	$script_path/fastq.sh $R1 $seq_file $fastq $run_info $fastqc
    ILL2SANGER1=`perl $script_path/checkFastqQualityScores.pl $fastq/$R1 10000`

########################################################	
######		Run bwa alignemnt module
    if [ $ILL2SANGER1 -gt 65 ]
    then
        $bwa/bwa aln $parameters -I $genome_bwa $fastq/$R1 > $output_dir_sample/$sample.$SGE_TASK_ID.R$read.sai
    else
        $bwa/bwa aln $parameters $genome_bwa $fastq/$R1 > $output_dir_sample/$sample.$SGE_TASK_ID.R$read.sai
    fi        
  
    if [ ! -s $output_dir_sample/$sample.$SGE_TASK_ID.R$read.sai ]
    then
        $script_path/errorlog.sh $output_dir_sample/$sample.$SGE_TASK_ID.R$read.sai align_read_bwa.sh ERROR "not created"
        exit 1; 
    fi  
	echo `date`
fi
	
   
