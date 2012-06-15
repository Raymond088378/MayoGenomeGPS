#!/bin/sh

##	INFO
##	to add rsids to per sample report
	
###########################
#		$1		=		TempFolder
#		$2		=		snv input file
#		$3		=		indel input file
#		$3		=		chromomse index
#		$4		=		run info
###############################

if [ $# != 4 ];
then
    echo "Usage<TempReportDir> <snv file><chromosome> <run info>";
else	
    set -x
    echo `date`
    TempReports=$1
    snv=$2
    chr=$3
    run_info=$4
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    dbsnp_rsids_snv=$( cat $tool_info | grep -w '^dbSNP_SNV_rsIDs' | cut -d '=' -f2)
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
    dbsnp_rsids_disease=$( cat $tool_info | grep -w '^dbSNP_disease_rsIDs' | cut -d '=' -f2) 
    
    num=`cat $TempReports/$snv | wc -l` 
    cat $TempReports/$snv | awk 'NR>1' > $TempReports/$snv.forrsIDs
	len=`cat $TempReports/$snv.forrsIDs | wc -l`
	if [ $len -gt 1 ]
	then
		perl $script_path/add_dbsnp_snv.pl -i $TempReports/$snv.forrsIDs -b 1 -s $dbsnp_rsids_snv -c 1 -p 2 -o $TempReports/$snv.forrsIDs.added -r $chr -h 1 
		## add column to add flag for disease variant
		perl $script_path/add.dbsnp.disease.snv.pl -i $TempReports/$snv.forrsIDs.added -b 1 -s $dbsnp_rsids_disease -c 1 -p 2 -o $TempReports/$snv.forrsIDs.added.disease -r $chr
    else
		value=`echo $dbsnp_rsids_snv | perl -wlne 'print $1 if /.+dbSNP(\d+)/'`
		echo "dbsnp${value}\tdbsnp${value}Alleles" > $TempReports/$snv.forrsIDs.added
		cat $TempReports/$snv.forrsIDs | sed 's/[ \t]*$//' > $TempReports/$snv.forrsIDs.tmp
		mv $TempReports/$snv.forrsIDs.tmp $TempReports/$snv.forrsIDs
		paste $TempReports/$snv.forrsIDs $TempReports/$snv.forrsIDs.added > $TempReports/$snv.forrsIDs.added.tmp
		mv $TempReports/$snv.forrsIDs.added.tmp $TempReports/$snv.forrsIDs.added
		echo "DiseaseVariant" > $TempReports/$snv.forrsIDs.added.disease
		paste $TempReports/$snv.forrsIDs.added $TempReports/$snv.forrsIDs.added.disease > $TempReports/$snv.forrsIDs.added.disease.tmp
		mv $TempReports/$snv.forrsIDs.added.disease.tmp $TempReports/$snv.forrsIDs.added.disease 
	fi
	perl $script_path/extract.rsids.pl -i $TempReports/$snv -r $TempReports/$snv.forrsIDs.added.disease -o $TempReports/$snv.rsIDs -v SNV
    num_a=`cat $TempReports/$snv.rsIDs |wc -l `
    if [ $num == $num_a ]
    then
        rm $TempReports/$snv
        rm $TempReports/$snv.forrsIDs.added
        rm $TempReports/$snv.forrsIDs
        rm $TempReports/$snv.forrsIDs.added.disease
    fi
    echo `date`
fi	
	
	
