#!/bin/sh

if [ $# != 5 ]
then
    echo "Usage: <input bam> <outputbam><temp dir><indexing flag (true/false)><run info>"
else
    set -x
    echo `date`
    inbam=$1
    outbam=$2
    tmp_dir=$3
    index=`echo $4 | tr "[a-z]" "[A-Z]"`
    run_info=$5
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    reads=$( cat $tool_info | grep -w '^MAX_READS_MEM_SORT' | cut -d '=' -f2 )
    
    $java/java -Xmx6g -Xms512m \
    -jar $picard/MergeSamFiles.jar \
    $inbam \
    MAX_RECORDS_IN_RAM=$reads \
    OUTPUT=$outbam \
    TMP_DIR=$tmp_dir \
    CREATE_INDEX=false \
    VALIDATION_STRINGENCY=SILENT
    
    if [ ! -s $outbam ]
    then
        echo "ERROR: merging bams failed $outbam"
    else
        files=`echo $inbam | sed -e '/INPUT=/s///g'`
        indexes=`echo $inbam | sed -e '/INPUT=/s///g' | tr " " "\n" | awk '{print $0".bai"}'`
		rm $files $indexes
        #if [ $index == "TRUE" ]
        #then
            $samtools/samtools index $outbam
        #fi
    fi
    echo `date`
fi	