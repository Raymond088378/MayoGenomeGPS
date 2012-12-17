#!/bin/bash
#	INFO
#	wrapper for 2A and 2B

if [ $# -le 5 ]
then
    echo -e "wrapper script to do realignment and variant calling\nUsage: ./realign_recal.sh Multi-Samples\n<input ':' sep> <bam ':' sep[normal:tumor1:tumor2:tumorN]> <samples ':' sep[normal:tumor1:tumor2:tumorN]> <outputdir bams> <run_info><1 for realign-recal or 0 for recal-realign><SGE_TASK_ID (optional)>\nelse\n<input> <bam > <samples> <outputdir bams><run_info><1 for realign-recal or 0 for recal-realign><SGE_TASK_ID (optional)>"
else
    set -x
    echo `date`
    input=$1    
    bam=$2
    samples=$3
    output_bam=$4
    run_info=$5
    flag=$6
	if [ $7 ]
	then
		SGE_TASK_ID=$7
	fi	
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
	export PATH=$java:$PATH

	
    ### update dash board    
    if [ $SGE_TASK_ID == 1 ]
    then
        for i in `echo $samples | tr ":" " "`
		do
			$script_path/dashboard.sh $i $run_info Realignment started
			id=`echo $samples | awk -v sample=$i -F ':' '{ for(i=1;i<=NF;i++){ if ($i == sample) {print i} } }'`
			in=`echo $input | cut -d ":" -f "$id"`
			bb=`echo $bam | cut -d ":" -f "$id"`
			$script_path/filesize.sh Realignment $i $in $bb $run_info
		done
	fi    
	
	### check and validate the input files
	for i in `echo $samples | tr ":" " "`
	do
		id=`echo $samples | awk -v sample=$i -F ':' '{ for(i=1;i<=NF;i++){ if ($i == sample) {print i} } }'`
		in=`echo $input | cut -d ":" -f "$id"`
		bb=`echo $bam | cut -d ":" -f "$id"`
		$samtools/samtools view -H $in/$bb 1>$in/$bb.rr.$chr.header 2> $in/$bb.fix.rr.$chr.log
		if [[ `cat $in/$bb.fix.rr.$chr.log | wc -l` -gt 0 || `cat $in/$bb.rr.$chr.header | wc -l` -le 0 ]]
		then
			$script_path/email.sh $in/$bb "bam is truncated or corrupt" processBam.sh $run_info
			$script_path/wait.sh $in/$bb.fix.rr.$chr.log
		else
			rm $in/$bb.fix.rr.$chr.log
		fi
		rm $in/$bb.rr.$chr.header
	done		
	
    if [ $flag == 1 ]
    then
        $script_path/realign_per_chr.sh $input $bam $output_bam $run_info 0 1 $samples $chr
        $script_path/recal_per_chr.sh $output_bam chr${chr}.realigned.bam $output_bam $run_info 1 0 multi $chr
    else
        $script_path/recal_per_chr.sh $input $bam $output_bam $run_info 0 1 $samples $chr
        $script_path/realign_per_chr.sh $output_bam chr${chr}.recalibrated.bam $output_bam $run_info 1 0 multi $chr
    fi
	
    ## update the dash board
    if [ $SGE_TASK_ID == 1 ]
    then
		 for i in `echo $samples | tr ":" " "`
		do
			$script_path/dashboard.sh $i $run_info Realignment complete
		done
    fi    
	### file name will be chr*.cleaned.bam
	$samtools/samtools view -H $output_bam/chr$chr.cleaned.bam 1>$output_bam/chr$chr.cleaned.bam.rr.$chr.header 2>$output_bam/chr$chr.cleaned.bam.fix.rr.$chr.log
	if [[ `cat $output_bam/chr$chr.cleaned.bam.fix.rr.$chr.log | wc -l` -gt 0  || `cat $output_bam/chr$chr.cleaned.bam.rr.$chr.header | wc -l` -le 0 ]]
	then
		$script_path/errorlog.sh $output_bam/chr$chr.cleaned.bam realign_recal.sh ERROR "is truncated or corrupt"
		exit 1;
	else
		rm $output_bam/chr$chr.cleaned.bam.fix.rr.$chr.log
	fi	
	rm $output_bam/chr$chr.cleaned.bam.rr.$chr.header
	if [ `echo $samples | tr ":" "\n" | wc -l` -gt 1 ]
	then
		$script_path/filesize.sh Realignment multi_sample $output_bam chr$chr.cleaned.bam $run_info
	else
		$script_path/filesize.sh Realignment $samples $output_bam chr$chr.cleaned.bam $run_info
	fi
	echo `date`
fi	