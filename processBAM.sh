#!/bin/bash

########################################################
###### 	Merge BAMs for a sample,  FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			merge_align.bam.sh
######		Date:				07/27/2011
######		Summary:			Using PICARD to sort and mark duplicates in bam 
######		Input files:		$1	=	/path/to/input directory
######							$2	=	sample name
######							$3	=	/path/to/run_info.txt
######		Output files:		Sorted and clean BAM
######		TWIKI:				http://bioinformatics.mayo.edu/BMI/bin/view/Main/BioinformaticsCore/Analytics/WholeGenomeWo
########################################################

if [ $# != 3 ];
then
    echo -e "wrapper to merge bam files and validate the bam for downstream analysis\nprocessBAM.sh </path/to/input directory><sample name> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    input=$1
    sample=$2
    run_info=$3
	
########################################################	
######		Reading run_info.txt and assigning to variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2)
    reorder=$( cat $tool_info | grep -w '^REORDERSAM' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	dup_flag=$( cat $tool_info | grep -w '^REMOVE_DUP' | cut -d '=' -f2 )
	dup=$( cat $tool_info | grep -w '^MARKDUP' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
########################################################	

    ### dashboard update
    if [ $analysis == "realign-mayo" ]
	then
		$script_path/dashboard.sh $sample $run_info Alignment started
    fi
	INPUTARGS="";
    files=""
    indexes=""
    cd $input
    for file in $input/*sorted.bam
    do
		base=`basename $file`
		dir=`dirname $file`
        $script_path/filesize.sh processBAM $sample $dir $base $JOB_ID $run_info
		$samtools/samtools view -H $file 1>$file.header 2> $file.fix.log
		if [ `cat $file.fix.log | wc -l` -gt 0 ]
		then
			$script_path/email.sh $file "bam is truncated or corrupt" $run_info
			$script_path/wait.sh $file.fix.log
		else
			rm $file.fix.log 
		fi	
		rm $file.header
		INPUTARGS="INPUT="$file" "$INPUTARGS;
        files=$file" $files";
        indexes=${file}.bai" $indexes"
    done
    
	num_times=`echo $INPUTARGS | tr " " "\n" | grep -c -w 'INPUT'`
	if [ $num_times == 1 ]
	then
		bam=`echo $INPUTARGS | cut -d '=' -f2`
		mv $bam $input/$sample.bam
		if [ -s $bam.bai ]
		then
			rm $bam.bai
		fi	
		SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $input/$sample.bam -s $samtools`
		if [ $SORT_FLAG == 1 ]
		then
			mv $input/$sample.bam $input/$sample.sorted.bam
			$samtools/samtools index $input/$sample.sorted.bam
		else
			$script_path/sortbam.sh $input/$sample.bam $input/$sample.sorted.bam $input coordinate true $run_info
		fi
	else	
	    $script_path/MergeBam.sh "$INPUTARGS" $input/$sample.sorted.bam $input true $run_info
		for i in $indexes
		do
			if [ -s $i ]
			then
				rm $i
			fi
		done			
	fi
	
    ### add read grouup information
    RG_ID=`$samtools/samtools view -H $input/$sample.sorted.bam | grep "^@RG" | tr '\t' '\n' | grep "^ID"| cut -f 2 -d ":"`

    if [ "$RG_ID" == "$sample" ]
    then
		echo " [`date`] no need to convert same read group"
    else	
        $script_path/addReadGroup.sh $input/$sample.sorted.bam $input/$sample.sorted.rg.bam $input $run_info $sample
    fi
    
    if [ $dup == "YES" ]
    then
		DUP_STATUS=`$samtools/samtools view -H $input/$sample.sorted.bam | grep "^@CO" | grep "MarkDuplicates" | wc -l`
		if [ "$DUP_STATUS" -eq 0 ] 
		then
		    $script_path/rmdup.sh $input/$sample.sorted.bam $input/$sample.sorted.rmdup.bam $input/$sample.dup.metrics $input $dup_flag true true $run_info   
		fi
    fi
    ## reorder if required
    if [ $reorder == "YES" ]
    then
        $script_path/reorderBam.sh $input/$sample.sorted.bam $input/$sample.sorted.tmp.bam $input $run_info
    fi
    if [ $analysis == "realignment" -o $analysis == "realign-mayo" ]
    then
        $samtools/samtools flagstat $input/$sample.sorted.bam > $input/$sample.flagstat
    fi
    ### index the bam again to maintain the time stamp for bam and index generation for down stream tools
    if [ $input/$sample.sorted.bam -nt $input/$sample.sorted.bam.bai ]
    then
        $samtools/samtools index $input/$sample.sorted.bam
    fi
    ## dashboard
    $script_path/dashboard.sh $sample $run_info Alignment complete
	## size of the bam file
	$script_path/filesize.sh Alignment.out $sample $input $sample.sorted.bam $JOB_ID $run_info
    echo `date`
fi
