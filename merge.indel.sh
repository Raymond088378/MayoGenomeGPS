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
	    $script_path/errorlog.sh $TempReports/$indel_file.rsIDs.frequencies merge.indel.sh ERROR "failed to create"
		exit 1;
	fi	
	## add snpeff prediction
	perl $script_path/add_snpeff_indel.pl -i $TempReports/$indel_file.rsIDs.frequencies -s $snpeff/$sample.chr${which_chr}.indel.eff -o $TempReports/$sample.chr${which_chr}.INDEL.report 
	perl $script_path/add_snpeff_indel_filter.pl -i $TempReports/$indel_file.rsIDs.frequencies -s $snpeff/$sample.chr${which_chr}.indel.filtered.eff -o $TempReports/$sample.chr${which_chr}.filtered.INDEL.report 
	
	num=`cat $TempReports/$sample.chr${which_chr}.INDEL.report | awk '{print $1"_"$2"_"$3"_"$9"_"$10}' | sort | uniq | wc -l`
	num_b=`cat $TempReports/$sample.chr${which_chr}.filtered.INDEL.report  | wc -l `
	if [[ $num == $num_a  && $num_a == $num_b ]]
	then
	    rm $TempReports/$indel_file.rsIDs.frequencies
	else
	    $script_path/errorlog.sh $TempReports/$sample.chr${which_chr}.filtered.INDEL.report merge.indel.sh ERROR "failed to create" 
	fi	
	for report in $TempReports/$sample.chr${which_chr}.INDEL.report $TempReports/$sample.chr${which_chr}.filtered.INDEL.report
	do
		perl $script_path/add_entrezID.pl -i $report -m $GeneIdMap -o $report.entrezid
		mv $report.entrezid $report
		perl $script_path/to.exclude.redundant.columns.from.report.pl $report $report.formatted
		mv $report.formatted $report
	done
	perl $script_path/add.cols.pl $TempReports/$sample.chr${which_chr}.INDEL.report $run_info INDEL > $TempReports/$sample.chr${which_chr}.INDEL.xls
	perl $script_path/add.cols.pl $TempReports/$sample.chr${which_chr}.filtered.INDEL.report $run_info INDEL > $TempReports/$sample.chr${which_chr}.filtered.INDEL.xls
	rm $TempReports/$sample.chr${which_chr}.INDEL.report $TempReports/$sample.chr${which_chr}.filtered.INDEL.report
	echo `date`
fi	