#!/bin/bash

if [ $# != 3 ]
then
    echo -e "wrapper to merge bam files and validate the bam for downstream analysis\nUsage: ./reformat_BAM.sh </path/to/input directory> <sample name> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    input=$1
    sample=$2
    run_info=$3
	
########################################################	
######		Reading run_info.txt and assigning to variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    reorder=$( cat $tool_info | grep -w '^REORDERSAM' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	########################################################	
######		PICARD to merge BAM file
    INPUTARGS="";
    files=""
    cd $input
    for file in $input/*sorted.bam
    do
		$samtools/samtools view -H $file 1>$file.rf.header 2> $file.rf.fix.log
		if [[ `cat $file.rf.fix.log | wc -l` -gt 0 || `cat $file.rf.header | wc -l` -le 0 ]]
		then
			$script_path/email.sh $file "bam is truncated or corrupt" "primary_script" $run_info
			$script_path/wait.sh $file.rf.fix.log 
		else
			rm $file.rf.fix.log 
		fi
		rm $file.rf.header	
		INPUTARGS="INPUT="$file" "$INPUTARGS;
		files=$file" $files";
    done
    
    num_times=`echo $INPUTARGS | tr " " "\n" | grep -c -w 'INPUT'`
    if [ $num_times == 1 ]
    then
		bam=`echo $INPUTARGS | cut -d '=' -f2`
		mv $bam $input/$sample.bam
        SORT_FLAG=`$script_path/checkBAMsorted.pl -i $input/$sample.bam -s $samtools`
        if [ $SORT_FLAG == 1 ]
        then
            mv $input/$sample.bam $input/$sample.sorted.bam
            $samtools/samtools index $input/$sample.sorted.bam
        else
            $script_path/sortbam.sh $input/$sample.bam $input/$sample.sorted.bam $input coordinate $run_info
        fi
    else	
        $script_path/MergeBam.sh "$INPUTARGS" $input/$sample.sorted.bam $input true $run_info
    fi

    ### add read grouup information
    RG_ID=`$samtools/samtools view -H $input/$sample.sorted.bam | grep "^@RG" | tr '\t' '\n' | grep "^ID"| cut -f 2 -d ":"`

    if [ "$RG_ID" == "$sample" ]
    then
		echo " [`date`] no need to convert same read group"
    else	
        $script_path/addReadGroup.sh $input/$sample.sorted.bam $input/$sample.sorted.rg.bam $input $run_info $sample
    fi
	
    if [ $reorder == "YES" ]
    then
        $script_path/reoderBam.sh $input/$sample.sorted.bam $input/$sample.sorted.tmp.bam $input $run_info  
    fi
	if [ -s $input/$sample.sorted.bam ]
	then
		$script_path/filesize.sh Realignment $sample $input $sample.sorted.bam $run_info
	else
		$script_path/errorlog.sh reformat_BAM.sh $input/$sample.sorted.bam ERROR "not found"
		exit 1;
	fi	
    echo `date`
fi	
	
