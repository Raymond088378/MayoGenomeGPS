#!/bin/bash

if [ $# != 8 ]
then
    echo -e "script to run somtic indel caller\nUsage: ./somaticindel.sh <normal bam ><tumor bam><chromosome><range><tumor sample><output dir><output file vcf><run info>"
else
    set -x
    echo `date`
    tumor_bam=$1
    normal_bam=$2
    chr=$3
    param=$4
    tumor_sample=$5
    output=$6
    output_file=$7
    run_info=$8
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	command_line_params=$( cat $tool_info | grep -w '^SOMATIC_INDEL_params' | cut -d '=' -f2 )
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	mem=$( cat $memory_info | grep -w '^SomaticIndelDetector_JVM' | cut -d '=' -f2)
	
	export PATH=$java:$PATH
	indel_v=$tumor_sample.chr$chr.indel.txt
	let check=0
	let count=0
	if [ ! -d $output/temp ]
	then
		mkdir -p $output/temp
	fi
	
	$samtools/samtools view -H $tumor_bam 1>$tumor_bam.si.header 2> $tumor_bam.fix.si.log
	$samtools/samtools view -H $normal_bam 1>$normal_bam.si.header 2> $normal_bam.fix.si.log
	if [[ `cat $tumor_bam.fix.si.log | wc -l` -gt 0 || `cat $tumor_bam.si.header | wc -l` le 0 ]]
	then
		$script_path/email.sh $tumor_bam "bam is truncated or corrupt" realign_recal.sh $run_info
		$script_path/wait.sh $tumor_bam.fix.si.log 
	else
		rm $tumor_bam.fix.si.log
	fi	
	rm $tumor_bam.si.header
	if [[ `cat $normal_bam.fix.si.log | wc -l` -gt 0 || `cat $normal_bam.si.header | wc -l` -le 0 ]]
	then
		$script_path/email.sh $normal_bam "bam is truncated or corrupt" realign_recal.sh $run_info
		$script_path/wait.sh $normal_bam.fix.si.log 
	else
		rm $normal_bam.fix.si.log
	fi	
	rm $normal_bam.si.header
	while [[ $check -eq 0 && $count -le 10 ]]
    do
		$java/java $mem -Djava.io.tmpdir=$output/temp/ -jar $gatk/GenomeAnalysisTK.jar \
		-R $ref \
		-et NO_ET \
		-K $gatk/Hossain.Asif_mayo.edu.key \
		-T SomaticIndelDetector \
		$range \
		-o $output/$output_file \
		-verbose $output/$indel_v \
		-I:normal $normal_bam \
		-I:tumor $tumor_bam $command_line_params
		sleep 5
		check=`[ -s $output/$output_file.idx ] && echo "1" || echo "0"`
        if [ $check -eq 0 ]
        then
            if [[  `find . -name '*.log'` ]]
			then
				if [ `grep -l $output/$output_file *.log` ]
				then
					rm `grep -l $output/$output_file *.log`
					rm core.*
				fi	
			fi
		fi
		let count=count+1	
    done 
	
    if [ ! -s $output/$output_file.idx ]
    then
        $script_path/errorlog.sh $output/$output_file somaticindel.sh ERROR "failed to create"
		exit 1;
    else
        if [ `cat $output/$output_file | awk '$0 !~ /^#/' | wc -l ` -eq 0 ]
        then
            $script_path/errorlog.sh $output/$output_file somaticindel.sh WARNING "no calls"
        fi 
		$script_path/convertvcf.pl $output/$output_file > $output/$output_file.tmp
		mv $output/$output_file.tmp $output/$output_file
		$script_path/fixindelAD.pl $output/$output_file $output/$output_file.temp
        mv $output/$output_file.temp $output/$output_file
        rm $output/$indel_v
        if [ `cat $output/$output_file |  awk '$0 ~ /^##FORMAT=<ID=GT/' | wc -l` == 0 ]
        then
            $script_path/add.gt.to.vcf.pl $output/$output_file | awk '$0 ~ /^#/ || $8 ~ /SOMATIC/' > $output/$output_file.tmp
        else
            cat $output/$output_file | awk '$0 ~ /^#/ || $8 ~ /SOMATIC/' > $output/$output_file.tmp
		fi
        cat $output/$output_file.tmp | awk '$0 ~ /^#/ || $5 ~ /,/' > $output/$output_file.multi.vcf
		cat $output/$output_file.tmp | awk '$0 ~ /^#/ || $5 !~ /,/' > $output/$output_file
		rm $output/$output_file.tmp $output/$output_file.idx
    fi
    echo `date`
fi
