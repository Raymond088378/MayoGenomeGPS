#!/bin/bash
if [ $# != 3 ];
then
    echo -e "Usage: wrapper to merge bam files and validate the bam for downstream analysis \n merge_align.bam.sh </path/to/input directory> <name of BAM to sort> <sample name> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    input=$1
    group=$2
    run_info=$3
	
########################################################	
######		Reading run_info.txt and assigning to variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    reorder=$( cat $tool_info | grep -w '^REORDERSAM' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    ref_path=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
    max_files=$( cat $tool_info | grep -w '^MAX_FILE_HANDLES' | cut -d '=' -f2 )
    max_reads=$( cat $tool_info | grep -w '^MAX_READS_MEM_SORT' | cut -d '=' -f2 )
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    samples=$(cat $sample_info | grep -w "^$group" | cut -d '=' -f2 | tr "\t" " ")	
	########################################################	
######		PICARD to merge BAM file
    INPUTARGS="";
    files=""
    cd $input
    for file in $input/*sorted.bam
    do
		$samtools/samtools view -H $file 1>$file.rf.header 2> $file.fix.log
		if [ `cat $file.fix.log | wc -l` -gt 0 ]
		then
			$script_path/email.sh $file "bam is truncated or corrupt" $JOB_NAME $JOB_ID $run_info
			while [ -f $file.fix.log ]
			do
				echo "waiting for the file to be fixed"
				sleep 2m
			done
		else
			rm $file.fix.log 
		fi	
		rm $file.rf.header
		INPUTARGS="INPUT="$file" "$INPUTARGS;
		files=$file" $files";
    done
    
    num_times=`echo $INPUTARGS | tr " " "\n" | grep -c -w 'INPUT'`
    if [ $num_times == 1 ]
    then
		bam=`echo $INPUTARGS | cut -d '=' -f2`
		mv $bam $input/$group.bam
        SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $input/$group.bam -s $samtools`
        if [ $SORT_FLAG == 1 ]
        then
            mv $input/$group.bam $input/$group.sorted.bam
            $samtools/samtools index $input/$group.sorted.bam
        else
            $script_path/sortbam.sh $input/$group.bam $input/$group.sorted.bam $input coordinate $max_reads $run_info
        fi
    else	
        $script_path/MergeBam.sh $INPUTARGS $input/$group.sorted.bam $input true $run_info
    fi

    ### add read grouup information
    RG_ID=`$samtools/samtools view -H $input/$group.sorted.bam | grep "^@RG" | tr '\t' '\n' | grep "^ID"| cut -f 2 -d ":"`
    
    check=0
    for sample in $samples
    do
        value=`echo $RG_ID | grep -w "$sample" | wc -l`;
        check=`if [ $value -ge 1 ]; then echo "1"; else echo "0"; fi`
    done 
    
    if [ $check == 0 ]
    then
        $script_path/errorlog.sh $input/$group.sorted.bam reformat_pairBAM.sh ERROR "Read group Not available"
        exit 1;
    fi        
	
    if [ $reorder == "YES" ]
    then
        $script_path/reorderBam.sh $input/$group.sorted.bam $input/$group.sorted.tmp.bam $input $run_info  
    fi
	
    echo `date`
fi	
	