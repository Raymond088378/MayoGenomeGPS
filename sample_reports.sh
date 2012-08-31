#!/bin/bash
##	INFO
##	to merge the scripts to call in one script to control sapnning jobs
if [ $# -le 7 ]
then
	echo "Usage: <run info> <sample> <Temp reports> <Output OnTarget> <sift> <snpeff> <polyphen> <output dir><group optional>";
else
	set -x
	echo `date`			
	run_info=$1 
	sample=$2
	TempReports=$3 
	output_OnTarget=$4 
	sift=$5 
	snpeff=$6
	poly=$7
	output_dir=$8
	if [ $9 ]
	then
		gr=$9
	fi	
	#SGE_TASK_ID=22
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	which_chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
	
	## prepocessing the input file from variant module or user added 
	if [[ $9 ]]
	then
		if [ ! -s $output_OnTarget/$gr.$sample.variants.chr$which_chr.SNV.filter.i.c.vcf  ]
		then
			touch $output_OnTarget/$gr.$sample.variants.chr$which_chr.SNV.filter.i.c.vcf.fix.log
			$script_path/email.sh $output_OnTarget/$gr.$sample.variants.chr$which_chr.SNV.filter.i.c.vcf  "not found" $JOB_NAME $JOB_ID $run_info
			while [ -f $output_OnTarget/$gr.$sample.variants.chr$which_chr.SNV.filter.i.c.vcf.fix.log ]
			do
				echo "waiting for the fix"
				sleep 2m
			done
		fi			
		$script_path/parse.vcf.sh $output_OnTarget/$gr.$sample.variants.chr$which_chr.SNV.filter.i.c.vcf $TempReports/$gr.$sample.chr${which_chr}.snv $run_info SNV
		##INDEL
		if [ ! -s $output_OnTarget/$gr.$sample.variants.chr$which_chr.INDEL.filter.i.c.vcf  ]
		then
			touch $output_OnTarget/$gr.$sample.variants.chr$which_chr.INDEL.filter.i.c.vcf.fix.log
			$script_path/email.sh $output_OnTarget/$gr.$sample.variants.chr$which_chr.INDEL.filter.i.c.vcf  "not found" $JOB_NAME $JOB_ID $run_info
			while [ -f $output_OnTarget/$gr.$sample.variants.chr$which_chr.INDEL.filter.i.c.vcf.fix.log ]
			do
				echo "waiting for the fix"
				sleep 2m
			done
		fi	
		$script_path/parse.vcf.sh $output_OnTarget/$gr.$sample.variants.chr$which_chr.INDEL.filter.i.c.vcf $TempReports/$gr.$sample.chr${which_chr}.indel $run_info INDEL
		sam=$gr.$sample
	else
		if [ ! -s $output_OnTarget/$sample.variants.chr$which_chr.SNV.filter.i.c.vcf  ]
		then
			touch $output_OnTarget/$sample.variants.chr$which_chr.SNV.filter.i.c.vcf.fix.log
			$script_path/email.sh $output_OnTarget/$sample.variants.chr$which_chr.SNV.filter.i.c.vcf  "not found" $JOB_NAME $JOB_ID $run_info
			while [ -f $output_OnTarget/$sample.variants.chr$which_chr.SNV.filter.i.c.vcf.fix.log ]
			do
				echo "waiting for the fix"
				sleep 2m
			done
		fi		
		$script_path/parse.vcf.sh $output_OnTarget/$sample.variants.chr$which_chr.SNV.filter.i.c.vcf $TempReports/$sample.chr${which_chr}.snv $run_info SNV
		##INDEL
		if [ ! -s $output_OnTarget/$sample.variants.chr$which_chr.INDEL.filter.i.c.vcf  ]
		then
			touch $output_OnTarget/$sample.variants.chr$which_chr.INDEL.filter.i.c.vcf.fix.log
			$script_path/email.sh $output_OnTarget/$sample.variants.chr$which_chr.INDEL.filter.i.c.vcf  "not found" $JOB_NAME $JOB_ID $run_info
			while [ -f $output_OnTarget/$sample.variants.chr$which_chr.INDEL.filter.i.c.vcf.fix.log ]
			do
				echo "waiting for the fix"
				sleep 2m
			done
		fi	
		$script_path/parse.vcf.sh $output_OnTarget/$sample.variants.chr$which_chr.INDEL.filter.i.c.vcf $TempReports/$sample.chr${which_chr}.indel $run_info INDEL	
		sam=$sample
	fi	
	
	if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
	then
		snv_var=$sam.chr${which_chr}.snv
		## add rsids
		$script_path/add.rsids_snvs.sh $TempReports $snv_var $which_chr $run_info
		## add allele frequency
		$script_path/add.frequencies.sh $TempReports $snv_var $which_chr $run_info
		## merge sift snpeff polyphen codon UCSC tracks and other tracks
		$script_path/merge.snv.sh $TempReports $sam $which_chr $sift $snpeff $poly $snv_var $run_info
	fi
	
	if [ $variant_type == "BOTH" -o $variant_type = "INDEL" ]
	then	
		indel_var=$sam.chr${which_chr}.indel
		## add rsIDs
		$script_path/add.rsids_indels.sh $TempReports $indel_var $which_chr $run_info
		$script_path/merge.indel.sh $TempReports $sam $which_chr $snpeff $indel_var $run_info
	fi
    echo `date`
fi	
    



            
