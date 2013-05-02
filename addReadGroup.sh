#!/bin/bash

### Called by processBAM.sh
if [ $# != 6 ]
then
    echo -e "script to add read group and platform information to a BAM file\
		\nUsage: ./addReadGroup.sh </fullpath/to/input bam> <fullpath/to/outputbam></path/to/temp dir></path/to/tool info></path/to/memory info><sample name>"
else
    set -x
    echo `date`
    inbam=$1
    outbam=$2
    tmp_dir=$3
    tool_info=$4
    memory_info=$5
    sample=$6
    
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	params=$( cat $tool_info | grep -w '^PICARD_ReadGroup_params' |sed -e '/PICARD_ReadGroup_params=/s///g')
	mem=$( cat $memory_info | grep -w '^READGROUP_JVM' | cut -d '=' -f2)
	
	
    $java/java $mem -jar $picard/AddOrReplaceReadGroups.jar \
    INPUT=$inbam \
    OUTPUT=$outbam \
    SM=$sample ID=$sample PU=$sample \
    TMP_DIR=$tmp_dir \
    VALIDATION_STRINGENCY=SILENT $params
    file=`echo $outbam | sed -e 's/\(.*\)..../\1/'`
	mv $file.bai $file.bam.bai
    
    if [ -s $outbam ]
    then
        $samtools/samtools view -H $outbam 1> $outbam.rg.header 2> $outbam.fix.rg.log
		if [[ `cat $outbam.fix.rg.log | wc -l` -gt 0 || `cat $outbam.rg.header | wc -l` -le 0 ]]
		then
			$script_path/errorlog.sh $outbam addReadGroup.sh ERROR "truncated or corrupted "
			exit 1;
		else
			rm $outbam.fix.rg.log
			mv $outbam $inbam
			mv $outbam.bai $inbam.bai
		fi	
		rm $outbam.rg.header
    else
        $script_path/errorlog.sh $outbam convert.bam.sh ERROR "is empty"
		exit 1;
	fi
    echo `date`
fi	