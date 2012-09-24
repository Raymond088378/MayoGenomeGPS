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
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	
    $java/java -Xmx5g -Xms512m \
    -jar $picard/MarkDuplicates.jar \
    INPUT=$inbam \
    OUTPUT=$outbam \
    METRICS_FILE=$metrics \
    ASSUME_SORTED=$sorted \
    REMOVE_DUPLICATES=$remove \
    MAX_FILE_HANDLES=$files \
    CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
    CO=MarkDuplicates \
    TMP_DIR=$tmp_dir
	file=`echo $outbam | sed -e 's/\(.*\)..../\1/'`
	mv $file.bai $file.bam.bai
	
    if [ -s $outbam ]
    then
       $samtools/samtools view -H $outbam 2> $outbam.log
		if [ `cat $outbam.log | wc -l` -gt 0 ]
		then
			$script_path/errorlog.sh rmdup.sh $outbam ERROR "truncated or corrupted"
			exit 1;
		else
			rm $outbam.log
		fi	
		
		mv $outbam $inbam
        if [ $index == "FALSE" ]
        then
            rm $file.bam.bai
        else
			mv $outbam.bai $inbam.bai
		fi
    else
        $script_path/errorlog.sh rmdup.sh $outbam ERROR empty
		exit 1;
    fi
    echo `date`
fi	