#!/bin/bash

if [ $# -le 4 ]
then
	echo -e "Usage: SCRIPT to split the bam uisng read group information \n</path/to/realign dir> </path/to/output folder> <sample> </path/to/alignment folder></path/to/run info>";
else
    set -x
    echo `date`
    input=$1
    output=$2
    sample=$3
    alignment=$4
    run_info=$5
    if [ $6 ]
	then
		SGE_TASK_ID=$6
	fi	
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    cd $input/$sample
    pair=$( cat $sample_info | grep -w "^$sample" | cut -d '=' -f2 | tr "\t" " ")
    if [ -f $input/$sample/chr$chr.cleaned.bam ]
    then
        $samtools/samtools view -H chr$chr.cleaned.bam 1>$output/$sample.chr$chr.header.sam 2>$output/$sample.chr$chr.cleaned.bam.fix.ssp.log
 		if [ `cat $output/$sample.chr$chr.cleaned.bam.fix.ssp.log | wc -l` -gt 0 ]
		then
			$script_path/email.sh $input/$sample/chr$chr.cleaned.bam "bam is truncated or corrupt" $JOB_NAME $JOB_ID $run_info
			while [ -f $output/$sample.chr$chr.cleaned.bam.fix.ssp.log ]
			do
				echo "waiting for the $input/$sample/chr$chr.cleaned.bam to be fixed"
				sleep 2m
			done
		else
			rm $output/$sample.chr$chr.cleaned.bam.fix.ssp.log
		fi
        for i in $pair
        do
            sam=`echo $pair | tr " " "\n" | grep -v $i | tr "\n" " "`
            gr=""
            for s in $sam
            do
                a="ID:$s|";
                gr="$gr$a"
            done
            gr=`echo $gr |  sed "s/|$//"`
            cat $output/$sample.chr$chr.header.sam | grep -w -E -v "$gr" > $output/$sample.chr$chr.$i.header.sam
            $samtools/samtools view -b -r $i $input/$sample/chr$chr.cleaned.bam > $output/$sample.$i.chr$chr.bam
            $samtools/samtools reheader $output/$sample.chr$chr.$i.header.sam $output/$sample.$i.chr$chr.bam > $output/$sample.$i.chr$chr.re.bam
            mv $output/$sample.$i.chr$chr.re.bam $output/$sample.$i.chr$chr.bam
            rm $output/$sample.chr$chr.$i.header.sam
        done
        rm $output/$sample.chr$chr.header.sam
    fi
    echo `date`
fi