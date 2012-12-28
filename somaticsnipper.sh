#!/bin/bash

if [ $# != 8 ]
then
    echo -e "script to run somatic snipper on the tumor normal bam file and outputs a vcf file\nUsage: ./somaticsipper.sh <normal bam> <tumor bam > <output dir> <chromosome> <tumor sample name> <normal sample name> <output vcf file name> <run info>"
else
    set -x
    echo `date`
    normal_bam=$1
    tumor_bam=$2
    output=$3
    chr=$4
    tumor_sample=$5
    normal_sample=$6
    output_file=$7
    run_info=$8
	
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    somatic_sniper=$( cat $tool_info | grep -w '^SOMATIC_SNIPER' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	command_line_params=$( cat $tool_info | grep -w '^SOMATIC_SNIPER_params' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 )
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
	
	snv=$tumor_sample.chr$chr.snv.output
	
	$samtools/samtools view -H $tumor_bam 1> $tumor_bam.$chr.header 2> $tumor_bam.$chr.fix.ss.log
	$samtools/samtools view -H $normal_bam 1>$normal_bam.$chr.header 2> $normal_bam.$chr.fix.ss.log
	if [[ `cat $tumor_bam.$chr.fix.ss.log | wc -l` -gt 0 || `cat $tumor_bam.$chr.header | wc -l` -le 0 ]]
	then
		$script_path/email.sh $tumor_bam "bam is truncated or corrupt" realign_recal.sh $run_info
		$script_path/wait.sh $tumor_bam.$chr.fix.ss.log
	else
		rm $tumor_bam.$chr.fix.ss.log
	fi	
	rm $tumor_bam.$chr.header
	if [[ `cat $normal_bam.$chr.fix.ss.log | wc -l` -gt 0 || `cat $normal_bam.$chr.header | wc -l` -le 0 ]]
	then
		$script_path/email.sh $normal_bam "bam is truncated or corrupt" realign_recal.sh $run_info
		$script_path/wait.sh $normal_bam.$chr.fix.ss.log
	else
		rm $normal_bam.$chr.fix.ss.log
	fi	
	rm $normal_bam.$chr.header
	### removing duplicates from the bam files
	#$samtools/samtools view -b -F 1024 $tumor_bam > $output_file/$tumor_sample.chr$chr.ss.bam
	#$samtools/samtools view -b -F 1024 $normal_bam > $output_file/$normal_sample.chr$chr.ss.bam
	
	#normal_bam=$output_file/$normal_sample.chr$chr.ss.bam
	#tumor_bam=$output_file/$tumor_sample.chr$chr.ss.bam
    $somatic_sniper/bam-somaticsniper $command_line_params -F vcf -f $ref $tumor_bam $normal_bam $output/$snv
    #rm $normal_bam $tumor_bam
	
	cat $output/$snv | awk 'BEGIN {OFS="\t"} {if($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,"PASS",$8,$9,$10,$11;}' | sed -e "/NORMAL/s//$normal_sample/g" | sed -e "/TUMOR/s//$tumor_sample/g"  | awk '$0 ~ /^#/ || $5 !~ /,/' | $script_path/ssniper_vcf_add_AD.pl > $output/$output_file
    cat $output/$snv | awk 'BEGIN {OFS="\t"} {if($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,"PASS",$8,$9,$10,$11;}' | sed -e "/NORMAL/s//$normal_sample/g" | sed -e "/TUMOR/s//$tumor_sample/g"  | awk '$0 ~ /^#/ || $5 ~ /,/' | $script_path/ssniper_vcf_add_AD.pl > $output/$output_file.multi.vcf    
    
    only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
    TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
    
    if [ $only_ontarget == "YES" ]
    then
        len=`cat $output/chr$chr.target.bed |wc -l`
        if [ $len -gt 0 ]
        then
            $bedtools/intersectBed -a $output/$output_file -b $output/chr$chr.target.bed -wa -header > $output/$output_file.i
            mv $output/$output_file.i $output/$output_file
        fi    
    fi
    	
    if [ ! -s $output/$output_file ]
    then
        $script_path/errorlog.sh $output/$output_file somaticsnipper.sh ERROR "failed to create"
		exit 1;
    else
        rm $output/$snv   
    fi
    echo `date`
fi	
