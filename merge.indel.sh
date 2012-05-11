#!/bin/sh
#	INFO
#	script to make per sample indel report

if [ $# != 6 ]; 
then
	echo "Usage <TempReports> <sample> <which_chr> <sseq> <indel_file><run info>";
else
	set -x
	echo `date`
	TempReports=$1 
	sample=$2 
	which_chr=$3
	snpeff=$4 
	indel_file=$5
	run_info=$6
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	GeneIdMap=$( cat $tool_info | grep -w '^GeneIdMap' | cut -d '=' -f2)
	GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
	cosmic=$( cat $tool_info | grep -w '^COSMIC_INDEL_REF' | cut -d '=' -f2) 
	## add cosmic data
	num=`cat $TempReports/$indel_file.rsIDs |wc -l `
	cat $TempReports/$indel_file.rsIDs | awk 'NR>1' | cut -f 1,2,3,4,5 > $TempReports/$indel_file.rsIDs.forfrequencies.temp
	perl $script_path/add.cosmic.pl $TempReports/$indel_file.rsIDs.forfrequencies.temp 0 $cosmic $GenomeBuild 1 $TempReports/$indel_file.cosmic.txt
	rm $TempReports/$indel_file.rsIDs.forfrequencies.temp
	perl $script_path/extract.allele_freq.pl -i $TempReports/$indel_file.rsIDs -f $TempReports/$indel_file.cosmic.txt -o $TempReports/$indel_file.rsIDs.frequencies -v INDEL
	num_a=`cat $TempReports/$indel_file.rsIDs.frequencies |wc -l `
	if [ $num == $num_a ]
	then
		rm $TempReports/$indel_file.rsIDs
		rm $TempReports/$indel_file.cosmic.txt
	else
		echo "ERROR: adding cosmic data for indel failed for $indel_file"
	fi	
	## add snpeff prediction
	perl $script_path/add_snpeff_indel.pl $TempReports/$indel_file.rsIDs.frequencies 
	perl $script_path/MergeIndelReport_SSeq.pl $TempReports/$indel_file.rsIDs.frequencies $sseq/$sample.chr${which_chr}.indels.sseq $chr $pos $ref $alt > $TempReports/$sample.chr${which_chr}.INDEL.report
	rm $TempReports/$indel_file.cosmic.txt
	report=$TempReports/$sample.chr${which_chr}.INDEL.report
	perl $script_path/add_entrezID.pl -i $report -m $GeneIdMap -o $report.entrezid
	mv $report.entrezid $report
	perl $script_path/to.exclude.redundant.columns.from.report.pl $report $report.formatted
	mv $report.formatted $report
	echo `date`
fi	