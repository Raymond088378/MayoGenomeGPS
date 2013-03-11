#!/bin/bash

########################################################
###### 	SV CALLER FOR TUMOR/NORMAL PAIR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			run_crest_multi.sh
######		Date:				09/26/2011
######		Summary:			Calls Crest
######		Input 
######		$1	=	group name
######		$2	=	bam list (first is normal) : separated
######		$3	=	names of the samples : separated
######		$4	=	/path/to/output directory
######		$5	=	/path/to/run_info.txt
########################################################

if [ $# -le 3 ]
then
	echo -e "Script to run crest on a paired sample\nUsage: ./run_crest_multi.sh <group name> </path/to/input directory> </path/to/output directory> </path/to/run_info.txt><SGE_TASK_ID (optional)>"
else
	set -x
	echo `date`
	group=$1
	input=$2
	output_dir=$3
	run_info=$4
	if [ $5 ]
	then
		SGE_TASK_ID=$5
	fi	
	########################################################	
	######		Reading run_info.txt and assigning to variables
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	crest=$( cat $tool_info | grep -w '^CREST' | cut -d '=' -f2 )
	picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
	perllib=$( cat $tool_info | grep -w '^PERLLIB' | cut -d '=' -f2 )
	blat=$( cat $tool_info | grep -w '^BLAT' | cut -d '=' -f2 )
	cap3=$( cat $tool_info | grep -w '^CAP3' | cut -d '=' -f2 )
	blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
	blat_ref=$( cat $tool_info | grep -w '^BLAT_REF' | cut -d '=' -f2 )
	blat_server=$( cat $tool_info | grep -w '^BLAT_SERVER' | cut -d '=' -f2 )
	chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
	ref_genome=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
	samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2)
	min_read=$( cat $tool_info | grep -w '^STRUCT_MIN_SUPPORT' | cut -d '=' -f2)
	min_id=$( cat $tool_info | grep -w '^STRUCT_MIN_IDENTITY' | cut -d '=' -f2)
	blacklist_sv=$( cat $tool_info | grep -w '^BLACKLIST_SV' | cut -d '=' -f2 )
	bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
	crest_params=$( cat $tool_info | grep -w '^CREST_params' | cut -d '=' -f2 )

	########################################################	
	######		

	pid=""
	export PERL5LIB=$perllib
	PATH=$PATH:$blat:$crest:$perllib
	mkdir -p $output_dir/$group/log
	range=20000
	let blat_port+=$RANDOM%range
	status=`$blat/gfServer status localhost $blat_port | wc -l`;
	if [ "$status" -eq 0 ]
	then
		$blat/gfServer start $blat_server $blat_port -log=$output_dir/$group/log/blat.$group.$chr.txt $blat_ref  &
		pid=$!
		sleep 5m
	fi

	let num_tumor=`echo $samples|tr " " "\n"|wc -l`-1
	normal_sample=`echo $samples| tr " " "\n" | head -n 1 `
	tumor_list=`echo $samples | tr " " "\n" | tail -$num_tumor`

	i=2
	for file in $tumor_list
	do
		sample=`echo $samples| tr " " "\n"| head -n $i | tail -1`
		if [ ! -f $output_dir/$group/$sample.chr$chr.cover ]
		then
			$script_path/errorlog.sh $output_dir/$group/$sample.chr$chr.cover run_crest_multi.sh ERROR "not exist" 
			exit 1;
		fi

		if [ ! -f $output_dir/$group/$sample.chr$chr.sclip.txt ]
		then
			$script_path/errorlog.sh $output_dir/$group/$sample.chr$chr.sclip.txt run_crest_multi.sh ERROR "not exist" 
			exit 1;
		fi
		status=`$blat/gfServer status $blat_server $blat_port | wc -l`;
		let count=0
		while [[ "$status" -le 1 && $count -le 5 ]]
		do
			`kill -9 $pid`
			blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
			range=20000
			let blat_port+=$RANDOM%range
			status=`$blat/gfServer status $blat_server $blat_port | wc -l`;
			if [ "$status" -le 1 ]
			then
				rm $output_dir/$sample/log/blat.$sample.$chr.txt
				$blat/gfServer start $blat_server $blat_port -log=$output_dir/$group/log/blat.$group.$chr.txt $blat_ref  &
				pid=$!
				sleep 5m
			fi
			status=`$blat/gfServer status $blat_server $blat_port | wc -l`;
			let count=count+1
		done 	
		
		if [ $count -ge 5 ]
		then
			$script_path/errorlog.sh GFSERVER run_single_crest.sh ERROR "failed to create gfserver"
			exit 1;
		fi	
	
		$crest/CREST.pl -f $output_dir/$group/$sample.chr$chr.cover \
        -d $output_dir/$group/${file}.chr$chr.bam -g $output_dir/$group/${normal_sample}.chr$chr.bam \
		--ref_genome $ref_genome -t $blat_ref \
		--blatport $blat_port -blatserver $blat_server \
		--cap3 $cap3/cap3 -r chr$chr -o $output_dir/$group -p $sample.$chr $crest_params

        if [ -f $output_dir/$group/${file}.$chr.predSV.txt ] 
        then
            rm $output_dir/$group/${file}.chr$chr.bam $output_dir/$group/${file}.chr$chr.bam.bai
        fi
        $script_path/CREST2VCF.pl -i $output_dir/$group/${file}.$chr.predSV.txt -f $ref_genome -o $output_dir/$group/$file.$chr.vcf -s $file -t $samtools
        $script_path/vcfsort.pl ${ref_genome}.fai $output_dir/$group/$file.$chr.vcf > $output_dir/$group/$file.$chr.vcf.sort
        mv $output_dir/$group/$file.$chr.vcf.sort $output_dir/$group/$file.$chr.vcf
		if [ ! -s $output_dir/$group/$file.$chr.vcf.fail ]
        then
            rm $output_dir/$group/$file.$chr.vcf.fail
        fi  
		
		if [ -s $output_dir/$group/${file}.$chr.predSV.txt ]
        then
            awk "((\$10>=$min_read)&&(\$11>=$min_read)&&(\$14>=$min_id)&&(\$16>=$min_id))" $output_dir/$group/${file}.$chr.predSV.txt | awk '{print $1"\t"$2"\t"$2+1"\t"$5"\t"$6"\t"$6+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $output_dir/$sample/$sample.$chr.predSV.txt > $output_dir/$group/${file}.$chr.final.predSV.txt
		else
            touch $output_dir/$group/${file}.$chr.final.predSV.txt
        fi
		$script_path/CREST2VCF.pl -i $output_dir/$group/${file}.$chr.final.predSV.txt -f $ref_genome -o $output_dir/$group/${file}.$chr.final.vcf -s $file -t $samtools
		if [ ! -s $output_dir/$group/${file}.$chr.final.vcf.fail ]
        then
            rm $output_dir/$group/${file}.$chr.final.vcf.fail
        fi  
		$script_path/vcfsort.pl ${ref_genome}.fai $output_dir/$group/${file}.$chr.final.vcf > $output_dir/$group/${file}.$chr.final.vcf.sort
		mv $output_dir/$group/${file}.$chr.final.vcf.sort $output_dir/$group/${file}.$chr.final.vcf
		### vcf converter for CREST output to VCF
		if [ ! -f $output_dir/$group/${file}.$chr.final.predSV.txt ]
		then
			$script_path/errorlog.sh $output_dir/$group/${file}.$chr.final.predSV.txt run_crest_multi.sh ERROR "failed to create" 
			exit 1;
		else
			rm $output_dir/$group/$file.chr$chr.sclip.txt $output_dir/$group/$file.chr$chr.cover
		fi
		let i=i+1
	done
    rm $output_dir/$group/${normal_sample}.chr$chr.bam $output_dir/$group/${normal_sample}.chr$chr.bam.bai
    rm $output_dir/$group/${normal_sample}.chr$chr.sclip.txt $output_dir/$group/${normal_sample}.chr$chr.cover
	`kill -9 $pid`
    echo `date`
fi
