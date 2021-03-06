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

if [ $# -le 3 ];
then
    echo -e "wrapper to add read group and sort and reorder the bam \
		\nUsage: ./convert_sam_bam.sh </path/to/input directory> <input name of BAM >\
			 <sample name> </path/to/run_info.txt> <sge task id>";
	exit 1;
fi
    set -x
    echo `date`
    input=$1
    input_bam=$2
    sample=$3
	run_info=$4
	if [ $5 ]
	then
		id=$5
	fi
	
########################################################	
######		Reading run_info.txt and assigning to variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    ref_path=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	aligner=$( cat $run_info | grep -w '^ALIGNER' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]")
	usenovosort=$( cat $tool_info | grep -w '^USENOVOSORT' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]")
	novosort=$( cat $tool_info | grep -w '^NOVOSORT' | cut -d '=' -f2)
	
	if [ $aligner == "bwa" ]
	then
		previous="align_bwa.sh"
	else
		previous="align_novo.sh"
	fi	
    
########################################################	
######		PICARD/novosort to sort raw BAM file

    ## check if BAM is sorted
    $samtools/samtools view -H $input/$input_bam 1> $input/$input_bam.header 2> $input/$input_bam.fix.log
	if [[ `cat $input/$input_bam.fix.log | wc -l` -gt 0 || `cat $input/$input_bam.header | wc -l` -le 0 ]]
	then
		$script_path/email.sh $input/$input_bam "bam is truncated or corrupt" $previous $run_info
		$script_path/wait.sh $input/$input_bam.fix.log
	else
		rm $input/$input_bam.fix.log
	fi	
	rm $input/$input_bam.header
	SORT_FLAG=`$script_path/checkBAMsorted.pl -i $input/$input_bam -s $samtools`
    if [ $SORT_FLAG == 1 ]
    then
        ln -s $input/$input_bam $input/$sample.sorted.bam
    else
    	$script_path/sortbam.sh $input/$input_bam $input/$sample.sorted.bam $input coordinate true $tool_info $memory_info yes no 
    fi

#############################################################	
######		PICARD to check availability of ReadGroup and platform info

    RG_ID=`$samtools/samtools view -H $input/$sample.sorted.bam | grep "^@RG" | tr '\t' '\n' | grep "^ID"| cut -f 2 -d ":"`
    sam=`echo $sample | awk -F'.' '{print $1}'`
    if [ "$RG_ID" != "$sam" ]
    then
        $script_path/addReadGroup.sh $input/$sample.sorted.bam $input/$sample.sorted.rg.bam $input $tool_info $memory_info $sample
    fi
    ### flagstat on each bam file
    $script_path/indexbam.sh $input/$sample.sorted.bam $tool_info
	$script_path/flagstat.sh $input/$sample.sorted.bam $input/$sample.flagstat $tool_info samtools
    if [ ! -s $input/$sample.flagstat ]
    then
        $script_path/errorlog.sh convert_bam.sh $input/$sample.flagstat ERROR "empty"
		exit 1;
	else
		## update secondary dashboard
		$script_path/filesize.sh alignment $sam $input $sample.sorted.bam $run_info
		$script_path/dashboard.sh $sample $run_info Alignment complete $id
    fi
	echo `date`
fi
