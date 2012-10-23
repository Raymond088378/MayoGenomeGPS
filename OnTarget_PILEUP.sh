#!/bin/bash
##	INFO
#	To Intersect pileup with OnTarget Kit by splitting the bam file into 200 files

######################################
#		$1		=	input folder (realignment sample folder)
#		$2		=	chromsome index
#		$3		=	Ontarget output folder
#		$4		=	sample name
#		$5		=	run info file
#########################################

if [ $# -le 3 ];
then
    echo -e "SCRIPT to get Ontarget pileup\nUsage: ./OnTarget_PILEUP.sh </path/to/input dir> </path/to/output Ontarget> <sample> </path/to/run ifno><SGE_TASK_ID(optional)>";
else	
    set -x
    echo `date`
    input=$1
    output=$2
    sample=$3
    run_info=$4
    
    if [ $5 ]
    then
		SGE_TASK_ID=$5
    fi
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    CaptureKit=$( cat $tool_info | grep -w '^CAPTUREKIT' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    gene_body=$( cat $tool_info | grep -w '^MATER_GENE_BODY' | cut -d '=' -f2 )
	multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	mem=$( cat $memory_info | grep -w '^CoverageBySample_JVM' | cut -d '=' -f2)
   
    if [ $tool == "whole_genome" ]
    then
        kit=$gene_body
    else
        kit=$CaptureKit
    fi
    cat $kit | grep -w chr$chr > $output/$sample.chr$chr.bed
	if [ `cat $kit | wc -l` -gt 0 ]
	then
		param="$output/$sample.chr$chr.bed"
	else
		param="chr$chr"
	fi
	
    bam=$input/chr$chr.cleaned.bam 
	$samtools/samtools view -H $bam 1> $bam.OnTarget_PILEUP.header 2> $bam.fix.OnTarget_PILEUP.log
	if [[ `cat $bam.fix.OnTarget_PILEUP.log | wc -l` -gt 0 || `cat $bam.OnTarget_PILEUP.header | wc -l` -le 0 ]]
	then
		$script_path/email.sh $bam "bam is truncated or corrupt" realign_recal.sh $run_info
		$script_path/wait.sh $bam.fix.OnTarget_PILEUP.log
	else
		rm $bam.fix.OnTarget_PILEUP.log
	fi	
	rm $bam.OnTarget_PILEUP.header
    mkdir -p $output/temp
	
    if [ $multi == "YES" ]
    then
        $java/java $mem -Djava.io.tmpdir=$output/temp/ -jar \
		$gatk/GenomeAnalysisTK.jar \
		-et NO_ET \
		-K $gatk/Hossain.Asif_mayo.edu.key \
		-T CoverageBySample  \
		-I $bam -R $ref \
		-L $param -o $output/$sample.chr$chr.txt
		pair=$( cat $sample_info | grep -w "^$sample" | cut -d '=' -f2 | tr "\t" " ")
		for i in $pair
        do
            for((j=0; j<=99; j++))
            do
				a=`cat $output/$sample.chr$chr.txt | grep -w "$i" | awk '$NF>'$j'' | wc -l`
				echo $a >> $output/$sample.$i.chr$chr.pileup.i.out
			done
             if [ ! -s $output/$sample.$i.chr$chr.pileup.i.out ]
            then
                $script_path/errorlog.sh $output/$sample.$i.chr$chr.pileup.i.out OnTarget_PILEUP.sh ERROR "failed to create"
        fi    
        done    
    else	
		$java/java $mem -Djava.io.tmpdir=$output/temp/ -jar \
		$gatk/GenomeAnalysisTK.jar \
		-et NO_ET \
		-K $gatk/Hossain.Asif_mayo.edu.key \
		-T CoverageBySample  \
		-I $bam -R $ref \
		-L $param -o $output/$sample.chr$chr.txt
        for((j=0; j<=99; j++))
        do
			a=`cat $output/$sample.chr$chr.txt | grep -w "$sample" | awk '$NF>'$j''  | wc -l`
			echo $a >> $output/$sample.chr$chr.pileup.i.out
        done    
		if [ ! -s $output/$sample.chr$chr.pileup.i.out ]
		then
			$script_path/errorlog.sh $output/$sample.chr$chr.pileup.i.out OnTarget_PILEUP.sh ERROR "failed to create"
			exit 1;
		fi
	fi
    rm $output/$sample.chr$chr.txt $output/$sample.chr$chr.bed
	echo `date`
fi	
    
	
	
	