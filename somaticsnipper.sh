#!/bin/sh

if [ $# != 8 ]
then
    echo "Usage: <normal bam> <tumor bam > <output dir> <chromosome> <tumor sample name> <normal sample name> <output vcf file name> <run info>"
else
    set -x
    echo `date`
    normal_bam=$1
    tumor_bam=$2
    output=$3
    chr=$4
    tumor_sample=$5
    normal_sample=$6
    output_file=$7
    run_info=$8
	
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    somatic_sniper=$( cat $tool_info | grep -w '^SOMATIC_SNIPER' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    squal=$( cat $tool_info | grep -w 'SOMATIC_QUALITY' | cut -d '=' -f2)
    mqual=$( cat $tool_info | grep -w 'MAPPING_QUALITY' | cut -d '=' -f2)
    snv=$tumor_sample.chr$chr.snv.output
	
    $somatic_sniper/bam-somaticsniper -q $mqual -Q $squal -F vcf -f $ref $tumor_bam $normal_bam $output/$snv
    cat $output/$snv | sed -e "/NORMAL/s//$normal_sample/g" | sed -e "/TUMOR/s//$tumor_sample/g"  | awk '$0 ~ /^#/ || $5 !~ /,/' | $script_path/ssniper_vcf_add_AD.pl > $output/$output_file
    cat $output/$snv | sed -e "/NORMAL/s//$normal_sample/g" | sed -e "/TUMOR/s//$tumor_sample/g"  | awk '$0 ~ /^#/ || $5 ~ /,/' | $script_path/ssniper_vcf_add_AD.pl > $output/$output_file.multi.vcf    
	
    if [ ! -s $output/$output_file ]
    then
        echo "ERROR: $output/$output_file not found"
	exit 1;
    else
        rm $output/$snv   
    fi
    echo `date`
fi	