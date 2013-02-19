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
    echo -e "wrapper script to run the alignment using NOVO ALIGN\nUsage: ./align_novo.sh <sample name> </path/to/output_dir> </path/to/run_info.txt> <SGE TASK ID (optional)>";
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
    
########################################################	
######	Reading run_info.txt and assigning to variables
    seq_file=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    center=$( cat $tool_info | grep -w '^CENTER' | cut -d '=' -f2 )
    platform=$( cat $tool_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
    fastqc=$( cat $tool_info | grep -w '^FASTQC' | cut -d '=' -f2)
    genome_novo=$( cat $tool_info | grep -w '^NOVO_REF' | cut -d '=' -f2)    
    novoalign=$( cat $tool_info | grep -w '^NOVOALIGN' | cut -d '=' -f2)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    paired=$( cat $run_info | grep -w '^PAIRED' | cut -d '=' -f2)
	paramaters=$( cat $tool_info | grep -w '^NOVO_params' | cut -d '=' -f2)
########################################################	
######		Check FASTQ for Illumina or Sanger quality scrore
    
    output_dir_sample=$output_dir/alignment/$sample
    fastq=$output_dir/fastq
    fastqc=$output_dir/fastqc
    
	$script_path/dashboard.sh $sample $run_info Alignment started $SGE_TASK_ID
   
    if [ $paired == 1 ]
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
			$script_path/fastq.sh $i $seq_file $fastq $run_info $fastqc
	    	$script_path/filesize.sh alignment $sample $fastq $i $run_info
        done
    elif [ $paired == 0 ]
    then
        let fidx=$SGE_TASK_ID
        R1=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
        if [ ! -s $seq_file/$R1 ]
		then
			$script_path/errorlog.sh align_novo.sh $seq_file/$R1 ERROR "not found"
			exit 1;
		fi	
		$script_path/fastq.sh $R1 $seq_file $fastq $run_info $fastqc
		$script_path/filesize.sh alignment $sample $fastq $R1 $run_info
    fi  
    ILL2SANGER1=`perl $script_path/checkFastqQualityScores.pl $fastq/$R1 10000`

########################################################	
######		Run novoalign for PE or SR	
	if [ $ILL2SANGER1 -gt 65 ] 
	then
		qual="-F ILMFQ"
    else
        qual="-F STDFQ"
    fi
    if [ $paired == 1 ]
    then
        $novoalign $paramaters -d $genome_novo $qual -f $fastq/$R1 $fastq/$R2 \
            -o SAM "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:$platform\tCN:$center" |  $samtools/samtools view -bS - > $output_dir_sample/$sample.$SGE_TASK_ID.bam
    else
        $novoalign $paramaters -d $genome_novo $qual -f $fastq/$R1 \
            -o SAM "@RG\tID:$sample\tSM:$sample\tLB:$GenomeBuild\tPL:$platform\tCN:$center"  | $samtools/samtools view -bS - > $output_dir_sample/$sample.$SGE_TASK_ID.bam      
    fi    
    $script_path/filesize.sh alignment.out $sample $output_dir_sample $sample.$SGE_TASK_ID.bam $run_info
    $samtools/samtools view -H $output_dir_sample/$sample.$SGE_TASK_ID.bam 1>$output_dir_sample/$sample.$SGE_TASK_ID.bam.header 2> $output_dir_sample/$sample.$SGE_TASK_ID.bam.fix.log
    if [[ `cat $output_dir_sample/$sample.$SGE_TASK_ID.bam.fix.log | wc -l` -gt 0 || `cat $output_dir_sample/$sample.$SGE_TASK_ID.bam.header | wc -l` -le 0 ]]	
    then
            $script_path/errorlog.sh $output_dir_sample/$sample.$SGE_TASK_ID.bam align_novo.sh ERROR "truncated or corrupt"
            exit 1;
    else
            rm $output_dir_sample/$sample.$SGE_TASK_ID.bam.fix.log
            if [ $paired == 0 ]
            then
                rm $fastq/$R1
            else
                rm $fastq/$R1 $fastq/$R2
            fi    
    fi	
    rm $output_dir_sample/$sample.$SGE_TASK_ID.bam.header
    
    
	#############################################	
	### Sort BAM, adds RG & remove duplicates ###
	#############################################

	### TODO: Break Out Convert_Bam.SH and run as a separate task to avail different memory settings
	### than the novoalign or bwa process (for sorting)

    $script_path/convert_bam.sh $output_dir_sample $sample.$SGE_TASK_ID.bam $sample.$SGE_TASK_ID $SGE_TASK_ID $run_info
	if [ ! -s $output_dir_sample/$sample.$SGE_TASK_ID.flagstat ]
    then
        $script_path/errorlog.sh align_novo.sh $output_dir_sample/$sample.$SGE_TASK_ID.flagstat ERROR "empty"
		exit 1;
	fi
	echo `date`	
fi