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
    index=`echo $4 | tr "[A-Z]" "[a-z]" `
    run_info=$5
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    reads=$( cat $tool_info | grep -w '^MAX_READS_MEM_SORT' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
	
    $java/java -Xmx3g -Xms512m \
    -jar $picard/MergeSamFiles.jar \
    $inbam \
    MAX_RECORDS_IN_RAM=$reads \
    OUTPUT=$outbam \
    TMP_DIR=$tmp_dir \
    CREATE_INDEX=$index \
    VALIDATION_STRINGENCY=SILENT
    file=`echo $outbam | sed -e 's/\(.*\)..../\1/'`
	mv $file.bai $file.bam.bai
    
	if [ ! -s $outbam ]
    then
        $script_path/errorlog.sh $outbam MergeBam.sh ERROR empty
		exit 1;
	else
		$samtools/samtools view -H $outbam 1>$outbam.header 2>$outbam.fix.log
		if [ `cat $outbam.fix.log | wc -l` -gt 0 ]
		then
			$script_path/errorlog.sh $outbam MergeBam.sh ERROR "truncated or corrupt"
			exit 1;
		else
			rm $outbam.fix.log
		fi
		rm $outbam.header
        files=`echo $inbam | sed -e '/INPUT=/s///g'`
        indexes=`echo $inbam | sed -e '/INPUT=/s///g' | tr " " "\n" | awk '{print $0".bai"}'`
		for i in $files 
		do
			if [ -s $i ]
			then
				rm $i
			fi
		done
		for i in $indexes
		do
			if [ -s $i ]
			then
				rm $i 
			fi
		done	
    fi
    echo `date`
fi	