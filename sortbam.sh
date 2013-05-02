#!/bin/bash

### Called by ProcessBAM.sh, 
if [ $# != 9 ]
then
    echo -e "Script to sort the bam file using picard samtools \
		\nUsage: ./sortbam.sh <input bam> <outputbam></path/to/temp dir><sorting order> \
			<flag for indexing(true/false)></path/to/tool info></path/to/memory info><flag to remove inputbam(yes/no)><flag to mention if BAM is already sorted(yes/no)>"
else
    set -x
    echo `date`
    inbam=$1
    outbam=$2
    tmp_dir=$3
    order=$4
    index=`echo $5 | tr "[A-Z]" "[a-z]"`
    tool_info=$6
    memory_info=$7
    removeflag=`echo $8 | tr "[A-Z]" "[a-z]"`
    assume_sorted=`echo $9 | tr "[A-Z]" "[a-z]"`
    
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )	
    reads=$( cat $tool_info | grep -w '^MAX_READS_MEM_SORT' | cut -d '=' -f2 )
    mem=$( cat $memory_info | grep -w '^SORT_JVM' | cut -d '=' -f2)
    ### Included for future use 
    usenovosort=$( cat $tool_info | grep -w '^USENOVOSORT' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]")
	novosort=$( cat $tool_info | grep -w '^NOVOSORT' | cut -d '=' -f2)
	novosortopt=$( cat $tool_info | grep -w '^NOVOSORTPBOPT' | cut -d '=' -f2)
	
	if [ $usenovosort == "yes" ]
	then
		### NOVOSORT INJECTION
		if [ $assume_sorted == "yes" ]
		then
			if [ $order == "coordinate" ]
			then
				$novosort $novosortopt --assumesorted --index --tmpdir=$tmp_dir $inbam -o $outbam
			else
				$novosort $novosortopt --namesort --assumesorted --index --tmpdir=$tmp_dir $inbam -o $outbam
			fi	
		else
			if [ $order == "coordinate" ]
			then
				$novosort $novosortopt --index --tmpdir=$tmp_dir $inbam -o $outbam
			else
				$novosort $novosortopt --namesort --index --tmpdir=$tmp_dir $inbam -o $outbam
			fi		
		fi
	else	
		$java/java $mem -jar $picard/SortSam.jar \
		INPUT=$inbam \
		OUTPUT=$outbam \
		MAX_RECORDS_IN_RAM=$reads \
		SO=$order \
		TMP_DIR=$tmp_dir \
		CREATE_INDEX=$index \
		VALIDATION_STRINGENCY=SILENT	
		file=`echo $outbam | sed -e 's/\(.*\)..../\1/'`
		mv $file.bai $file.bam.bai
    fi
    if [ ! -s $outbam ]
    then
        $script_path/errorlog.sh sortbam.sh $outbam ERROR "is empty"
        exit 1;
    else
        $samtools/samtools view -H $outbam 1>$outbam.sort.header 2> $outbam.sort.log
		if [[ `cat $outbam.sort.log | wc -l` -gt 0 || `cat $outbam.sort.header| wc -l` -le 0 ]]
		then
			$script_path/errorlog.sh rmdup.sh $outbam ERROR "truncated or corrupted"
			exit 1;
		else
			rm $outbam.sort.log
		fi	
		if [ $removeflag == "yes" ]
		then
			rm $inbam 
		fi
		rm $outbam.sort.header
        if [ $index == "true" ]
        then
            echo "index already created"
        else
			rm $file.bam.bai
		fi    
    fi
    echo `date`
fi
	