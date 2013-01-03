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

if [ $# -le 4 ]
then
    echo -e "Script to run crest on a paired sample\nUsage: ./run_crest_multi_cover.sh <sample name> <group name> </path/to/input directory> </path/to/output directory> </path/to/run_info.txt>"
else
    set -x
    echo `date`
    sample=$1
    group=$2
    input=$3
    output_dir=$4
    run_info=$5
	if [ $6 ]
	then
		SGE_TASK_ID=$6
	fi	
	########################################################	
    ######		Reading run_info.txt and assigning to variables
    #SGE_TASK_ID=1
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    cap3=$( cat $tool_info | grep -w '^CAP3' | cut -d '=' -f2 )
	crest=$( cat $tool_info | grep -w '^CREST' | cut -d '=' -f2 )
    perllib=$( cat $tool_info | grep -w '^PERLLIB' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    ref_genome=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	somatic_calling=$( cat $tool_info | grep -w '^SOMATIC_CALLING' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	blat=$( cat $tool_info | grep -w '^BLAT' | cut -d '=' -f2 )
	blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
	blat_ref=$( cat $tool_info | grep -w '^BLAT_REF' | cut -d '=' -f2 )
	blat_server=$( cat $tool_info | grep -w '^BLAT_SERVER' | cut -d '=' -f2 )
	crest_params=$( cat $tool_info | grep -w '^CREST_params' | cut -d '=' -f2 )
	min_read=$( cat $tool_info | grep -w '^STRUCT_MIN_SUPPORT' | cut -d '=' -f2)
	min_id=$( cat $tool_info | grep -w '^STRUCT_MIN_IDENTITY' | cut -d '=' -f2)
	blacklist_sv=$( cat $tool_info | grep -w '^BLACKLIST_SV' | cut -d '=' -f2 )
	bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
	
	if [ $analysis == "variant" ]
	then
		input_bam=$input/$group.$sample.chr$chr.bam
		previous="split_sample_pair.sh"
	else
		input_bam=$input/$sample.sorted.bam
		previous="processBAM.sh"
	fi	
    export PERL5LIB=$perllib:$crest
    PATH=$PATH:$blat:$crest:$perllib
	mkdir -p $output_dir/$group/log
	$samtools/samtools view -H $input_bam 1>$input_bam.crest.$sample.$chr.header 2>$input_bam.crest.$sample.$chr.fix.log
	if [ `cat $input_bam.crest.$sample.$chr.fix.log | wc -l` -gt 0 ]
	then
		$script_path/email.sh $input_bam "bam is truncated or corrupt" $previous $run_info
		$script_path/wait.sh $input_bam.crest.$sample.$chr.fix.log
	else
		rm $input_bam.crest.$sample.$chr.fix.log
	fi
	rm $input_bam.crest.$sample.$chr.header
	if [ ! -s ${input_bam}.bai ]
	then
		$samtools/samtools index $input_bam 
	fi	
	$samtools/samtools view -b $input_bam chr$chr >  $output_dir/$group/$sample.chr$chr.bam
	$samtools/samtools index $output_dir/$group/$sample.chr$chr.bam
    file=$output_dir/$group/$sample.chr$chr.bam
    SORT_FLAG=`$script_path/checkBAMsorted.pl -i $file -s $samtools`
    if [ $SORT_FLAG == 0 ]
    then
        $script_path/errorlog.sh $file run_crest_multi_cover.sh ERROR "is not sorted"
		exit 1;
    fi
    # check if BAM has an index
    if [ ! -s $file.bai ]
    then
        $samtools/samtools index $file
    fi
	$crest/extractSClip.pl -i $file -r chr$chr --ref_genome $ref_genome -o $output_dir/$group -p $sample
    if [ $somatic_calling == "NO" ]
	then
		if [ -f $output_dir/$group/$sample.chr$chr.cover ]
		then
			mkdir -p $output_dir/$group/log/$sample/
			export TMPDIR=$output_dir/$sample/log/$sample/
			range=20000
			let blat_port+=$RANDOM%range
			status=`$blat/gfServer status localhost $blat_port | wc -l`;
			let count=0
			if [ "$status" -le 1 ]
			then
				$blat/gfServer start localhost $blat_port -log=$output_dir/$group/log/$sample/blat.$sample.$chr.txt $blat_ref  &
				pid=$!
				sleep 5m
			fi
			status=`$blat/gfServer status localhost $blat_port | wc -l`;
			while [[ "$status" -le 1 && $count -le 5 ]]
			do
				`kill -9 $pid`
				blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
				range=20000
				let blat_port+=$RANDOM%range
				status=`$blat/gfServer status localhost $blat_port | wc -l`;
				if [ "$status" -le 1 ]
				then
					rm $output_dir/$group/log/$sample/blat.$sample.$chr.txt
					$blat/gfServer start localhost $blat_port -log=$output_dir/$group/log/$sample/blat.$sample.$chr.txt $blat_ref  &
					pid=$!
					sleep 5m
				fi
				status=`$blat/gfServer status localhost $blat_port | wc -l`;
				let count=count+1
			done 	
		
			if [ $count -ge 5 ]
			then
				$script_path/errorlog.sh GFSERVER run_single_crest.sh ERROR "failed to create gfserver"
				exit 1;
			fi	
			cd $crest    
			$crest/CREST.pl -f $output_dir/$group/$sample.chr$chr.cover -d $file \
			--ref_genome $ref_genome -t $blat_ref \
			--blatport $blat_port -blatserver localhost \
			--cap3 $cap3/cap3 \
			-o $output_dir/$group/ -p $sample.$chr $crest_params
			if [ -f $output_dir/$group/${sample}.$chr.predSV.txt ] 
			then
				rm $file ${file}.bai
			fi
			$script_path/CREST2VCF.pl -i $output_dir/$group/${sample}.$chr.predSV.txt -f $ref_genome -o $output_dir/$group/${sample}.$chr.raw.vcf -s $sample -t $samtools
			if [ ! -s $output_dir/$group/${sample}.$chr.raw.vcf.fail ]
			then
				rm $output_dir/$group/${sample}.$chr.raw.vcf.fail
			fi  
			$script_path/vcfsort.pl ${ref_genome}.fai $output_dir/$group/${sample}.$chr.raw.vcf > $output_dir/$group/${sample}.$chr.raw.vcf.sort
			mv $output_dir/$group/${sample}.$chr.raw.vcf.sort $output_dir/$group/${sample}.$chr.raw.vcf
			awk "((\$10>=$min_read)&&(\$11>=$min_read)&&(\$14>=$min_id)&&(\$16>=$min_id))" $output_dir/$group/${sample}.$chr.predSV.txt | awk '{print $1"\t"$2"\t"$2+1"\t"$5"\t"$6"\t"$6+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $output_dir/$group/${sample}.$chr.predSV.txt > $output_dir/$group/${sample}.$chr.filter.predSV.txt
			$script_path/CREST2VCF.pl -i $output_dir/$group/${sample}.$chr.filter.predSV.txt -f $ref_genome -o $output_dir/$group/${sample}.$chr.filter.vcf -s $sample -t $samtools
			if [ ! -s $output_dir/$group/${sample}.$chr.filter.vcf.fail ]
			then
				rm $output_dir/$group/${sample}.$chr.filter.vcf.fail
			fi  
			$script_path/vcfsort.pl ${ref_genome}.fai $output_dir/$group/${sample}.$chr.filter.vcf > $output_dir/$group/${sample}.$chr.filter.vcf.sort
			mv $output_dir/$group/${sample}.$chr.filter.vcf.sort $output_dir/$group/${sample}.$chr.filter.vcf
			if [ -f $output_dir/$group/${sample}.$chr.predSV.txt ]
			then
				rm $output_dir/$group/$sample.chr$chr.cover $output_dir/$group/$sample.chr$chr.sclip.txt 
			else
				$script_path/errorlog.sh $output_dir/$group/${sample}.$chr.predSV.txt run_single_crest.sh ERROR "failed to create"
				exit 1;
			fi			
		else
			$script_path/errorlog.sh $output_dir/$group/$sample.chr$chr.cover run_single_crest.sh ERROR "failed to create"
			touch $output_dir/$group/${sample}.$chr.predSV.txt
			$script_path/CREST2VCF.pl -i $output_dir/$group/${sample}.$chr.predSV.txt -f $ref_genome -o $output_dir/$group/${sample}.$chr.raw.vcf -s $sample -t $samtools
			rm $output_dir/$group/${sample}.$chr.raw.vcf.fail
			touch  $output_dir/$group/${sample}.$chr.filter.predSV.txt
			$script_path/CREST2VCF.pl -i $output_dir/$group/${sample}.$chr.filter.predSV.txt -f $ref_genome -o $output_dir/$group/${sample}.$chr.filter.vcf -s $sample -t $samtools
			rm $output_dir/$group/${sample}.$chr.filter.vcf.fail
		fi
		`kill -9 $pid`	
	fi
	echo `date`
fi
