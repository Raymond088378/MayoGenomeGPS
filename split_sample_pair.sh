#!/bin/sh

if [ $# != 5 ]
then
	echo -e "Usage: SCRIPT to create IGV BAM \n</path/to/realign dir> </path/to/output folder> <sample> </path/to/alignment folder><run ifno>";
else
	set -x
	echo `date`
	input=$1
    output=$2
	sample=$3
	alignment=$4
    run_info=$5
#SGE_TASK_ID=1
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )

	cd $input/$sample
    pair=$( cat $sample_info | grep -w "$sample" | cut -d '=' -f2)
	if [ -f $input/$sample/chr$chr.cleaned.bam ]
	then
		$samtools/samtools view -H chr$chr.cleaned.bam > $output/$sample.chr$chr.header.sam
		for i in $pair
		do
			sam=`echo $pair | tr " " "\n" | grep -v $i | tr "\n" " "`
			gr=""
			for s in $sam
			do
				a="ID:$s|";
				gr="$gr $a"
			done
			gr=`echo $gr |  sed "s/|$//"`
			cat $output/$sample.chr$chr.header.sam |grep -E -v '$gr' > $output/$sample.chr$chr.$i.header.sam
			$samtools/samtools view -b -r $i $input/$sample//chr$chr.cleaned.bam > $output/$sample.$i.chr$chr.bam
			$samtools/samtools reheader $output/$sample.chr$chr.$i.header.sam $output/$sample.$i.chr$chr.bam > $output/$sample.$i.chr$chr.re.bam
			mv $output/$sample.$i.chr$chr.re.bam $output/$sample.$i.chr$chr.bam
			rm $output/$sample.chr$chr.$i.header.sam
		done
		rm $output/$sample.chr$chr.header.sam
	fi
	echo `date`
fi

		
            



