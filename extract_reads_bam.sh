#!/bin/bash

if [ $# -le 3 ]
then
	echo -e "Usage: script to extract reads not used for downstream preocessing \n <bam file><input directory><run info><igv folder><single/pair>"
else
	set -x
	echo `date`	
	output=$1
	bam=$2
	run_info=$3
	igv=$4
	if [ $5 ]
	then
		group=$5
	fi
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	chrindex=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | awk '{print "chr"$0}' )
	delivery_folder=$( cat $run_info | grep -w '^DELIVERY_FOLDER' | cut -d '=' -f2)
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]" )
	out=$delivery_folder/IGV_BAM
	chrs=`cat $ref.fai | cut -f1 | tr ":" "\n"`
	i=1
	for chr in $chrs
	do
		if [ `echo $chrindex | grep -w "$chr" | wc -l` -eq 0 ]
		then
			chrArray[$i]=$chr
			let i=i+1
		fi
	done

	$samtools/samtools view -H $output/$bam 2> $output/$bam.fix.log
	if [ `cat $output/$bam.fix.log | wc -l` -gt 0 ]
	then
		echo " $output/$bam :bam is truncated or corrupted"
		exit 1;
	else
		rm $output/$bam.fix.log
	fi
	
	## extract read for specfic chromosome
	if [ ! -s $output/$bam.bai ]
	then
		$samtools/samtools index $output/$bam
	fi	
	
	input=""
	for i in $(seq 1 ${#chrArray[@]})
	do
		chr=${chrArray[$i]}
		$samtools/samtools view -b $output/$bam $chr > $output/$bam.$chr.bam
		$samtools/samtools index $output/$bam.$chr.bam
		input=$input" INPUT=$output/$bam.$chr.bam"
	done

	### extract unmapped reads
	$samtools/samtools view -b -f 12 $output/$bam > $output/$bam.unmapped.bam
	$samtools/samtools index $output/$bam.unmapped.bam
	input=$input" INPUT=$output/$bam.unmapped.bam"

	$script_path/MergeBam.sh "$input" $output/$bam.extra.bam $output yes $run_info
	
	if [ $5 ]
	then
		sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
		samples=$(cat $sample_info | grep -w "^$group" | cut -d '=' -f2 | tr "\t" " ")	
		for sample in $samples
		do	
			sam=`echo $samples | tr " " "\n"| grep -w -v "$sample" | tr "\n" " "`
			gr=""
			for s in $sam
			do
				a="ID:$s|";
				gr="$gr$a"
			done
			gr=`echo $gr |  sed "s/|$//"`
			$samtools/samtools view -b -r $sample $output/$bam.extra.bam > $output/$sample.extra.bam
			$samtools/samtools view -H $output/$sample.extra.bam | grep -w -E -v "$gr" | $samtools/samtools reheader - $output/$sample.extra.bam > $output/$sample.extra.re.bam
			mv $output/$sample.extra.re.bam $output/$sample.extra.bam
			$samtools/samtools index $output/$sample.extra.bam
			if [ $delivery_folder != "NA" ]
			then
				if [ -d $delivery_folder ]
				then
					if [ ! -d $out ]
					then
						mkdir $out
					fi
					mv $output/$sample.extra.bam $out/
					mv $output/$sample.extra.bam.bai $out/
				fi
			else
				mv $output/$sample.extra.bam $igv/
				mv $output/$sample.extra.bam.bai $igv/
			fi		
		done
		rm $output/$bam.extra.bam $output/$bam.extra.bam.bai
	else
		if [ $delivery_folder != "NA" ]
		then
			if [ -d $delivery_folder ]
			then
				if [ ! -d $out ]
				then
					mkdir $out
				fi
				mv $output/$bam.extra.bam $out/
				mv $output/$bam.extra.bam.bai $out/
			fi
		else	
			mv $output/$bam.extra.bam $igv/
			mv $output/$bam.extra.bam.bai $igv/
		fi		
	fi
	echo `date`
fi 
 
	




