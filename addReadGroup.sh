#!/bin/sh

if [ $# != 5 ]
then
    echo "Usage: <input bam> <outputbam><temp dir><run info><sample>"
else
    set -x
    echo `date`
    inbam=$1
    outbam=$2
    tmp_dir=$3
    run_info=$4
    sample=$5
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    platform=$( cat $run_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
    center=$( cat $run_info | grep -w '^CENTER' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
                    
    $java/java -Xmx6g -Xms512m \
    -jar $picard/AddOrReplaceReadGroups.jar \
    INPUT=$inbam \
    OUTPUT=$outbam \
    PL=$platform SM=$sample CN=$center ID=$sample PU=$sample LB=$GenomeBuild \
    TMP_DIR=$tmp_dir \
    VALIDATION_STRINGENCY=SILENT		
    
    if [ -s $outbam ]
    then
        mv $outbam $inbam
        $samtools/samtools index $inbam
    else
        echo "ERROR: add read group failed for $inbam"
    fi
    echo `date`
fi	