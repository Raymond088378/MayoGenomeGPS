#!/bin/sh

##	INFO
##	to add rsids to per sample report
	
###########################
#		$1		=		TempFolder
#		$2		=		sample	
#		$3		=		chromossome
#		$4		=		run info
###############################

if [ $# != 4 ];
then
    echo "Usage<TempReportDir><variant file with rsids ><chromosome><run info> ";
else	
    set -x
    echo `date`
    TempReports=$1
    var=$2
    chr=$3
    run_info=$4
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    bgi=$( cat $tool_info | grep -w '^BGI_REF' | cut -d '=' -f2) 
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
    cosmic=$( cat $tool_info | grep -w '^COSMIC_SNV_REF' | cut -d '=' -f2)
    esp=$( cat $tool_info | grep -w '^ESP' | cut -d '=' -f2)
    
    num=`cat $TempReports/$var.rsIDs | wc -l `
    if [ $num -eq 0 ]
	then
		$script_path/errorlog.sh $TempReports/$var.rsIDs add.frequencies.sh ERROR "not created"
		exit 1;
	fi	
	
	cat $TempReports/$var.rsIDs | awk 'NR>1' | cut -f 1-5 > $TempReports/$var.forFrequencies.temp
	len=`cat $TempReports/$var.forFrequencies.temp | wc -l`
	if [ $len -gt 1 ]
	then
		##1KGenome and Hapmap
		$script_path/add_allele_freq_3populations.sh $TempReports/$var.forFrequencies.temp $TempReports/$var.forFrequencies.allele.frequency $run_info $chr
		## BGI
		perl $script_path/add_bgi_freq.pl -i $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.txt -r $bgi -c $chr -o $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.txt
		## Cosmic
		perl $script_path/add.cosmic.pl $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.txt 1 $cosmic $GenomeBuild 1 $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.txt
		## ESP
		perl $script_path/add_esp.pl -i $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.txt -d $esp -o $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt
    else
		touch $TempReports/$var.forFrequencies.allele.frequency.CEU
		touch $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI
		touch $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.txt
		touch $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.txt
		touch $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.txt
		echo -e "HapMap_CEU_allele_freq\t1kgenome_CEU_allele_freq\tHapMap_YRI_allele_freq\t1kgenome_YRI_allele_freq\tHapMap_JPT+CHB_allele_freq\t1kgenome_JPT+CHB_allele_freq\tBGI200_Danish\tCOSMIC\tESP5400_EUR_maf\tESP5400_AFR_maf" > $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt
		paste $TempReports/$var.forFrequencies.temp $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt > $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt.tmp
		mv $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt.tmp $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt
	fi
	### arrange the columns
    perl $script_path/extract.allele_freq.pl -i $TempReports/$var.rsIDs -f $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt -o $TempReports/$var.rsIDs.allele_frequencies -v SNV
    num_a=`cat $TempReports/$var.rsIDs.allele_frequencies | wc -l`
    if [ $num == $num_a ]
    then
        rm $TempReports/$var.rsIDs
        rm $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt
        rm $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.Cosmic.txt
        rm $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.BGI.txt
        rm $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI.CHBJPT.txt
        rm $TempReports/$var.forFrequencies.allele.frequency.CEU.YRI
        rm $TempReports/$var.forFrequencies.allele.frequency.CEU
        rm $TempReports/$var.forFrequencies.temp
    else
        $script_path/errorlog.sh $TempReports/$var.rsIDs.allele_frequencies add.frequencies.sh ERROR "failed to create"
		exit 1;
    fi    
    echo `date`
fi	