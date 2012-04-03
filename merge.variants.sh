#!/bin/sh
#	INFO
#	script merge all the variants per sample to be used in merged report

################################
#		$1		=		sift dir
#		$2		=		sseq dir
#		$3		=		chomosome index
#		$4		=		run information
#################################

if [ $# != 4 ];
then
	echo "Usage: <output dir><temp reports><chromossome><run_info>";
else
	set -x
	echo `date`
	output_dir=$1
	TempReports=$2
	which_chr=$3
	run_info=$4
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)	
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	SNV_caller=$( cat $run_info | grep -w '^SNV_CALLER' | cut -d '=' -f2 )
	variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
	
	if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
	then
		ls $output_dir/*.chr${which_chr}.raw.snvs.bed.i.ToMerge | sort > $TempReports/chr${which_chr}.snvs
		perl $script_path/snv.merger.pl -i $TempReports/chr${which_chr}.snvs -c quality -o $TempReports/list.chr${which_chr}.snvs
		rm $TempReports/chr${which_chr}.snvs
	fi
	if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
	then
		ls $output_dir/*.chr${which_chr}.raw.indels.bed.i.ToMerge | sort > $TempReports/chr${which_chr}.indels 
		perl $script_path/indel.merger.pl -i $TempReports/chr${which_chr}.indels -o $TempReports/list.chr${which_chr}.indels
		rm $TempReports/chr${which_chr}.indels
	fi
	echo `date`
fi