#!/bin/bash

########################################################
###### 	SORT BAM, adds RG & REMOVE DUPLICATES SCRIPT FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			convert.bam.sh
######		Date:				07/27/2011
######		Summary:			Using PICARD to sort and mark duplicates in bam 
######		Input files:		$1	=	/path/to/input directory
######							$2	=	name of BAM to sort
######							$3	=	sample name
######							$4	=	/path/to/run_info.txt
######		Output files:		Sorted and clean BAM
######		TWIKI:				http://bioinformatics.mayo.edu/BMI/bin/view/Main/BioinformaticsCore/Analytics/WholeGenomeWo
########################################################

if [ $# != 5 ];
then
    echo -e "Usage: wrapper to add read group and sort the bam </path/to/input directory> <name of BAM to sort> <sample name> <sge task id> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    input=$1
    input_bam=$2
    sample=$3
	id=$4
    run_info=$5
	
########################################################	
######		Reading run_info.txt and assigning to variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    dup=$( cat $run_info | grep -w '^MARKDUP' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    reorder=$( cat $run_info | grep -w '^REORDERSAM' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    ref_path=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    center=$( cat $run_info | grep -w '^CENTER' | cut -d '=' -f2 )
    version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
    platform=$( cat $run_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
    output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
    max_files=$( cat $tool_info | grep -w '^MAX_FILE_HANDLES' | cut -d '=' -f2 )
    max_reads=$( cat $tool_info | grep -w '^MAX_READS_MEM_SORT' | cut -d '=' -f2 )
    
########################################################	
######		PICARD to sort raw BAM file

    ## check if BAM is sorted
    $samtools/samtools view -H $input/$input_bam 2> $input/$input_bam.log
	if [ `cat $input/$input_bam.log | wc -l` -gt 0 ]
	then
		echo "$input/$input_bam : bam is corruped or truncated"
		exit 1;
	else
		rm $input/$input_bam.log
	fi	
	SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $input/$input_bam -s $samtools`
    if [ $SORT_FLAG == 1 ]
    then
        ln -s $input/$input_bam $input/$sample.sorted.bam
    else
        $script_path/sortbam.sh $input/$input_bam $input/$sample.sorted.bam $input coordinate $max_reads true $run_info
    fi

#############################################################	
######		PICARD to check availability of ReadGroup and platform info

    RG_ID=`$samtools/samtools view -H $input/$sample.sorted.bam | grep "^@RG" | tr '\t' '\n' | grep "^ID"| cut -f 2 -d ":"`
    sam=`echo $sample | awk -F'.' '{print $1}'`
    if [ "$RG_ID" != "$sam" ]
    then
        $script_path/addReadGroup.sh $input/$sample.sorted.bam $input/$sample.sorted.rg.bam $input $run_info $sample
    fi
    
    ### flagstat on each bam file
    $samtools/samtools index $input/$sample.sorted.bam
	$samtools/samtools flagstat $input/$sample.sorted.bam > $input/$sample.flagstat
    if [ ! -s $input/$sample.flagstat ]
    then
        $script_path/errorlog.sh convert.bam.sh $input/$sample.flagstat ERROR empty
		exit 1;
	fi
    ## update secondary dahboard
    size=`du -b $input/$sample.sorted.bam | sed 's/\([0-9]*\).*/\1/'`
	$script_path/filesize.sh alignment $sam $sample.sorted.bam $JOB_ID $size $run_info
	$script_path/dashboard.sh $sample $run_info Alignment complete $id
    echo `date`
fi
