#!/bin/bash


### Called by processBAM.sh
if [ $# != 8 ]
then
    echo -e "script to remove or flag the duplicates from a BAM file dending on the flag passed\nUsage: ./rmdup.sh <input bam> <outputbam><temp dir><max files to split><remove of flag dupluicate(true/false)><assume file is aorted or not(true/false)><do indexing or not(true/false)<run info>"
else
    set -x
    echo `date`
    inbam=$1
    outbam=$2
    metrics=$3
    tmp_dir=$4
    remove=$5
    sorted=$6
    index=`echo $7 | tr "[A-Z]" "[a-z]"`
    run_info=$8
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	max_files=$( cat $tool_info | grep -w '^MAX_FILE_HANDLES' | cut -d '=' -f2 )
    mem=$( cat $memory_info | grep -w '^RMDUP_JVM' | cut -d '=' -f2)
	
	$java/java $mem -jar $picard/MarkDuplicates.jar \
    INPUT=$inbam \
    OUTPUT=$outbam \
    METRICS_FILE=$metrics \
    ASSUME_SORTED=$sorted \
    REMOVE_DUPLICATES=$remove \
    MAX_FILE_HANDLES=$max_files \
    CREATE_INDEX=$index \
	VALIDATION_STRINGENCY=SILENT \
    CO=MarkDuplicates TMP_DIR=$tmp_dir
	file=`echo $outbam | sed -e 's/\(.*\)..../\1/'`
	mv $file.bai $file.bam.bai
	
    if [ -s $outbam ]
    then
       $samtools/samtools view -H $outbam 1>$outbam.rmdup.header 2> $outbam.rmdup.log
		if [[ `cat $outbam.rmdup.log | wc -l` -gt 0 || `cat $outbam.rmdup.header | wc -l` -le 0 ]]
		then
			$script_path/errorlog.sh rmdup.sh $outbam ERROR "truncated or corrupted"
			exit 1;
		else
			rm $outbam.rmdup.log
		fi	
		rm $outbam.rmdup.header
		mv $outbam $inbam
        if [ $index == "false" ]
        then
            rm $file.bam.bai
        else
			mv $outbam.bai $inbam.bai
		fi
    else
        $script_path/errorlog.sh rmdup.sh $outbam ERROR "is empty"
		exit 1;
    fi
    echo `date`
fi	