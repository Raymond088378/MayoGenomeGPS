#!/bin/sh
if [ $# != 3 ];
then
    echo -e "Usage: wrapper to merge bam files and validate the bam for downstream analysis \n merge_align.bam.sh </path/to/input directory> <name of BAM to sort> <sample name> </path/to/run_info.txt>";
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
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    reorder=$( cat $run_info | grep -w '^REORDERSAM' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    ref_path=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    center=$( cat $run_info | grep -w '^CENTER' | cut -d '=' -f2 )
    platform=$( cat $run_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
    max_files=$( cat $tool_info | grep -w '^MAX_FILE_HANDLES' | cut -d '=' -f2 )
    max_reads=$( cat $tool_info | grep -w '^MAX_READS_MEM_SORT' | cut -d '=' -f2 )
	
	########################################################	
######		PICARD to merge BAM file
    INPUTARGS="";
    files=""
    cd $input
    for file in $input/*sorted.bam
    do
      INPUTARGS="INPUT="$file" "$INPUTARGS;
      files=$file" $files";
    done
    
    num_times=`echo $INPUTARGS | tr " " "\n" | grep -c -w 'INPUT'`
    if [ $num_times == 1 ]
    then
	bam=`echo $INPUTARGS | cut -d '=' -f2`
	mv $bam $input/$sample.bam
        SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $input/$sample.bam -s $samtools`
        if [ $SORT_FLAG == 1 ]
        then
            mv $input/$sample.bam $input/$sample.sorted.bam
            $samtools/samtools index $input/$sample.sorted.bam
        else
            $script_path/sortbam.sh $input/$sample.bam $input/$sample.sorted.bam $input coordinate $max_reads $run_info
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
	
    echo `date`
fi	
	