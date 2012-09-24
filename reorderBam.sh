#!/bin/bash

if [ $# != 4 ]
then
    echo "Usage: <input bam> <outputbam><temp dir><run info>"
else
    set -x
    echo `date`
    inbam=$1
    outbam=$2
    tmp_dir=$3
    run_info=$4
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
    
    $java/java -Xmx6g -Xms512m \
    -jar $picard/ReorderSam.jar \
    I=$inbam \
    O=$outbam \
    R=$ref \
    TMP_DIR=$tmp_dir \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
    file=`echo $outbam | sed -e 's/\(.*\)..../\1/'`
    mv $file.bai $file.bam.bai
    if [ -s $outbam ]
    then
        mv $outbam $inbam
    else
        $script_path/errorlog.sh $outbam reorderBam.sh ERROR "failed to reorder"
        exit 1;
    fi
    echo `date`
fi	