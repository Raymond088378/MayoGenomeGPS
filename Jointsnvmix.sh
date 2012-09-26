#!/bin/sh

if [ $# != 8 ]
then
    echo "Usage: <normal bam> <tumor bam > <output dir> <chromosome> <tumor sample name> <normal sample name ><output file> <run info>"
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
    jointsnvmix=$( cat $tool_info | grep -w '^JOINTSNVMIX' | cut -d '=' -f2)    
    python=$( cat $tool_info | grep -w '^PYTHON' | cut -d '=' -f2) 
    pythonpath=$( cat $tool_info | grep -w '^PYTHONLIB' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2) 
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)	
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
	only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	command_line_params=$( cat $tool_info | grep -w '^JOINTSNVMIX_params' | cut -d '=' -f2 )
	bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
	JSM_Filter=$( cat $tool_info|grep -w 'JSM_Filter'|cut -d  '=' -f2 )
        

    export PYTHONPATH=$pythonpath:$PYTHONPATH
    export PATH=$python:$PYTHONPATH:$PATH
    
    if [ ! -s $normal_bam ]
    then
        $script_path/errorlogs.sh $normal_bam Jointsnvmix.sh ERROR "not exist"
        exit 1;
    else
    	$samtools/samtools view -H $normal_bam 1>$normal_bam.jsm.header 2> $normal_bam.fix.jsm.log
    	if [ `cat $normal_bam.fix.jsm.log | wc -l` -gt 0 ]
		then
			$script_path/email.sh $normal_bam "bam is truncated or corrupt" $JOB_NAME $JOB_ID $run_info
			while [ -f $normal_bam.fix.jsm.log ]
			do
				echo "waiting for the $normal_bam to be fixed"
				sleep 2m
			done
		else
			rm $normal_bam.fix.jsm.log
		fi	
		rm $normal_bam.jsm.header
    fi
    
    if [ ! -s $tumor_bam ]
    then
        $script_path/errorlogs.sh $tumor_bam Jointsnvmix.sh ERROR "not exist"
        exit 1;
    else
    	$samtools/samtools view -H $tumor_bam 1>$tumor_bam.jsm.header 2> $tumor_bam.fix.jsm.log
		if [ `cat $tumor_bam.fix.jsm.log | wc -l` -gt 0 ]
		then
			$script_path/email.sh $tumor_bam "bam is truncated or corrupt" $JOB_NAME $JOB_ID $run_info
			while [ -f $tumor_bam.fix.jsm.log ]
			do
				echo "waiting for the $tumor_bam to be fixed"
				sleep 2m
			done
		else
			rm $tumor_bam.fix.jsm.log
		fi	
		rm $tumor_bam.jsm.header
    fi
    
    ### run joint snvmix classify to call the somatic mutation
	$python/python $jointsnvmix/build/scripts-2.7/jsm.py classify --model snvmix2 $command_line_params --chromosome chr$chr --out_file $output/$output_file.txt --parameters_file $jointsnvmix/config/params.cfg $ref $normal_bam $tumor_bam
	
	### script to convert text output to vcf output 
	perl $script_path/jsm2vcf.pl -i $output/$output_file.txt -o $output/$output_file -ns $normal_sample -ts $tumor_sample $JSM_Filter
	cat $output/$output_file | awk '$0 ~ /^#/ || $5 ~ /,/' > $output/$output_file.multi.vcf
	cat $output/$output_file | awk '$0 ~ /^#/ || $5 !~ /,/' > $output/$output_file.tmp
	mv $output/$output_file.tmp $output/$output_file
	rm $output/$output_file.txt
        
	if [ $only_ontarget == "YES" ]
	then
		len=`cat $output/chr$chr.target.bed | wc -l`
		if [ $len -gt 0 ]
		then
			$bedtools/intersectBed -a $output/$output_file -b $output/chr$chr.target.bed -wa -header > $output/$output_file.i
			mv $output/$output_file.i $output/$output_file   
		fi
	fi
	if [ ! -s $output/$output_file ]
	then
		$script_path/errorlog.sh $output/$output_file Jointsnvmix.sh ERROR "not exist"
		exit 1;
	fi	
    echo `date`
fi
