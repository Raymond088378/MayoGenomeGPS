#!/bin/sh

if [ $# != 9 ]
then
    echo "Usage: <input bam> <outputbam><temp dir><max files to split><remove of flag dupluicate(true/false)><assume file is aorted or not(true/false)><do indexing or not(true/false)<run info>"
else
    set -x
    echo `date`
    inbam=$1
    outbam=$2
    metrics=$3
    tmp_dir=$4
    files=$5
    remove=$6
    sorted=$7
    index=`echo $8 | tr "[a-z]" "[A-Z]"`
    run_info=$9
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    
    $java/java -Xmx6g -Xms512m \
    -jar $picard/MarkDuplicates.jar \
    INPUT=$inbam \
    OUTPUT=$outbam \
    METRICS_FILE=$metrics \
    ASSUME_SORTED=$sorted \
    REMOVE_DUPLICATES=$remove \
    MAX_FILE_HANDLES=$files \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=$tmp_dir
    
    if [ -s $outbam ]
    then
        mv $outbam $inbam
        if [ $index == "TRUE" ]
        then
            $samtools/samtools index $inbam
        fi
    else
        echo "ERROR: [`date`] remove duplicate failed for $inbam"
    fi
    echo `date`
fi	