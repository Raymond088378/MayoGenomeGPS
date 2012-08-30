#!/bin/sh

if [ $# != 7 ]
then
    echo "Usage: <input bam> <outputbam><temp dir><sorting order><# of reads in RAM<run info>"
else
    set -x
    echo `date`
    inbam=$1
    outbam=$2
    tmp_dir=$3
    order=$4
    reads=$5
    index=`echo $6 | tr "[a-z]" "[A-Z]"`
    run_info=$7
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )	
    
    $java/java -Xmx6g -Xms512m \
    -jar $picard/SortSam.jar \
    INPUT=$inbam \
    OUTPUT=$outbam \
    MAX_RECORDS_IN_RAM=$reads \
    SO=$order \
    TMP_DIR=$tmp_dir \
    VALIDATION_STRINGENCY=SILENT	
    
    
    if [ ! -s $outbam ]
    then
        $script_path/errorlog.sh sortbam.sh $outbam ERROR empty
        exit 1;
    else
        $samtools/samtools view -H $outbam 2> $outbam.log
		if [ `cat $outbam.log | wc -l` -gt 0 ]
		then
			echo "$outbam :bam is truncated or corrupted "
			exit 1;
		else
			rm $outbam.log
		fi	
		rm $inbam
        if [ $index == "TRUE" ]
        then
            $samtools/samtools index $outbam
        fi    
    fi
    echo `date`
fi
	