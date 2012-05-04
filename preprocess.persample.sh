#!/bin/sh
##	INFO
##	script is a wrapper script to annotate variants per sample 

#########################################
#	$1		=		Sample Name
#	$2		=		SIFT directory	
#	$3		=		SSEQ directory	
#	$4		=		TempReports directory
#	$5		=		run information file to get filenames
#	$6		=		input directory run dir
#	$7		=		chromosome
###############################################	

if [ $# != 5 ];
then
    echo "Usage:<sample> <tempReport dir> <run info> <input variant folder><chromosome>";
else			
    set -x
    echo `date`
    sample=$1		
    TempReports=$2 	#Tempreport output folder
    run_info=$3
    input_dir=$4
    chr=$5
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2| tr "[a-z]" "[A-Z])
    
   
    if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
    then
        ## SNVs  parse the vcf file for on target
        var=$sample.chr$chr.SNV.filter.i.c.vcf
        perl $script_path/parse.vcf.SNV.pl -i $input_dir/$var -o $TempReports/${sample}.chr$chr.snv -s $sample -h 1
    fi
    if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
    then
        ## INDELs
        indel=$sample.chr$chr.INDEL.filter.i.c.vcf
        perl $script_path/parse.vcf.INDEL.pl -i $input_dir/$indel -o $TempReports/$sample.chr$chr.indel -s $sample -h 1
    fi
    echo `date`
fi	