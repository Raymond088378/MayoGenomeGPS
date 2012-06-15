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

if [ $# -le 4 ]
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
	if [ $6 ]
	then
		group=$6
    fi
	if [ $7 ]
	then
		prefix=$7
	fi
	
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
	
	if [[ $6 && $7 ]]
	then
		sam=$prefix.$group.$sample
		samp=$prefix.$group.$sample
	elif [ $6 ]
	then
		sam=$group.$sample
		samp=$group.$sample
	else
		sam=$sample
		samp=$sample
	fi	
	
   
    if [[ $variant_type == 'BOTH' || $variant_type == 'SNV' ]]
    then
        ## SNVs  parse the vcf file for on target
		var=$sam.variants.chr$chr.SNV.filter.i.c.vcf
		perl $script_path/parse.vcf.SNV.pl -i $input_dir/$var -o $TempReports/$samp.chr$chr.snv -s $sample -h 1
	fi
    if [[ $variant_type == 'BOTH' || $variant_type == 'INDEL' ]]
    then
        ## INDELs
        indel=$sam.variants.chr$chr.INDEL.filter.i.c.vcf
		perl $script_path/parse.vcf.INDEL.pl -i $input_dir/$indel -o $TempReports/$samp.chr$chr.indel -s $sample -h 1
	fi
    echo `date`
fi	