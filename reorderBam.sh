#!/bin/bash

if [ $# != 4 ]
then
    echo -e "script to reorder the bam file\nUsage: <input bam> <outputbam><temp dir><tool info><memory_info>"
else
    set -x
    echo `date`
    inbam=$1
    outbam=$2
    tmp_dir=$3
    tool_info=$4
    memory_info=$5
    
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	mem=$( cat $memory_info | grep -w '^REORDER_JVM' | cut -d '=' -f2)
    
    $java/java $mem -jar $picard/ReorderSam.jar \
    I=$inbam \
    O=$outbam \
    R=$ref \
    TMP_DIR=$tmp_dir \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
    file=`echo $outbam | sed -e 's/\(.*\)..../\1/'`
    mv $file.bai $file.bam.bai
    
    if [ -s $outbam ]
    then
    	$samtools/samtools view -H $outbam 1> $outbam.ro.header 2> $outbam.fix.ro.log
		if [ `cat $outbam.fix.ro.log | wc -l` -gt 0 ]
		then
			$script_path/errorlog.sh $outbam reorderBam.sh ERROR "truncated or corrupted "
			exit 1;
		else
			rm $outbam.fix.ro.log
        	mv $outbam $inbam
    	fi
    	rm $outbam.ro.header
	else
        $script_path/errorlog.sh $outbam reorderBam.sh ERROR "failed to reorder"
        exit 1;
    fi
    echo `date`
fi	