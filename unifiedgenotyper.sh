#!/bin/bash

if [ $# != 6 ]
then
    echo -e "script to run unified genotyper\nUsage: ./unifiedgenotyper.sh <bamfile> <vcf output> \
<type of variant> <range of positions> <output mode> <run info file>"
	exit 1;
fi
    set -x
    echo `date`
    bam=$1
    vcf=$2
    type=$3
    range=$4
    mode=$5
    run_info=$6

    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	ped=$( cat $tool_info | grep -w '^PEDIGREE' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	command_line_params=$( cat $tool_info | grep -w '^UnifiedGenotyper_params' | cut -d '=' -f2 )
	memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	mem=$( cat $memory_info | grep -w '^UnifiedGenotyper_JVM' | cut -d '=' -f2)
	 
	export PATH=$java:$PATH
	if [[ ${#ped} -eq 0 && $ped == "NA" ]]
    then
        ped="NA"
    fi
    
    let check=0
    out=`dirname $vcf`
    
    if [ ! -d $out/temp ]
	then
		mkdir -p $out/temp
		sleep 10s
	fi
	let count=0
	gatk_param="-R $ref -et NO_ET -K $gatk/Hossain.Asif_mayo.edu.key "
	while [[ $check -eq 0 && $count -le 10 ]]
    do
		if [ $ped != "NA" ]
		then
			$java/java $mem -Djava.io.tmpdir=$out/temp/ -jar $gatk/GenomeAnalysisTK.jar \
			-T UnifiedGenotyper --output_mode $mode -nt $threads \
			-glm $type $range $bam \
			--ped $ped --out $vcf $command_line_params $gatk_param
		else
			$java/java $mem -Djava.io.tmpdir=$out/temp/ -jar $gatk/GenomeAnalysisTK.jar \
			-T UnifiedGenotyper --output_mode $mode -nt $threads \
			-glm $type $range $bam \
			--out $vcf $command_line_params $gatk_param
		fi
		sleep 5
        check=`[ -s $vcf.idx ] && echo "1" || echo "0"`
        if [ $check -eq 0 ]
        then
			if [[  `find . -name '*.log'` ]]
			then
				if [ `grep -l $vcf *.log` ]
				then
					rm `grep -l $vcf *.log`
					rm core.*
				fi
			fi
		fi 
		let count=count+1	
    done 
	
	$script_path/convertvcf.pl $vcf > $vcf.tmp
	mv $vcf.tmp $vcf
	
    if [ ! -s $vcf.idx ]
    then
        $script_path/errorlog.sh $vcf unifiedgenotyper.sh ERROR "is empty"
        exit 1;
    else
		rm $vcf.idx	
	fi
    echo `date`

