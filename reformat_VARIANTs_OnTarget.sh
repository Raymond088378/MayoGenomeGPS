#!/bin/sh

if [ $# != 4 ]
then
	echo "usage: <output><sample><run_info><marker>";
else	
    set -x
    echo `date`
    output=$1
    sample=$2
    run_info=$3
    marker=$4
    
    input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
    variant_type=`echo "$variant_type" | tr "[a-z]" "[A-Z]"`
    dbsnp_rsids_indel=$( cat $tool_info | grep -w '^dbSNP_INDEL_rsIDs' | cut -d '=' -f2)
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" " " )
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    vcftools=$( cat $tool_info | grep -w '^VCFTOOLS' | cut -d '=' -f2)
    perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
    
    if [ $marker == 2 ]
    then
        snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
        indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)
            
        if [ $input/$snv_file == $input/$indel_file ]
        then
            for chr in $chrs
            do
                cat $input/$snv_file | awk -v num=chr${chr} '$0 ~ /#/ || $1 == num' > $output/$sample/$sample.variants.chr$chr.filter.vcf
            done
        else
            ### concat both the vcf files to get one vcf file
            for chr in $chrs
            do
                cat $input/$snv_file | awk -v num=chr${chr} '$0 ~ /#/ || $1 == num' > $output/$sample/$sample.variants.chr$chr.filter.SNV.vcf
                cat $input/$indel_file | awk -v num=chr${chr} '$0 ~ /#/ || $1 == num' > $output/$sample/$sample.variants.chr$chr.filter.INDEL.vcf
                in="-V $output/$sample/$sample.variants.chr$chr.filter.SNV.vcf -V $output/$sample/$sample.variants.chr$chr.filter.INDEL.vcf"
                $script_path/combinevcf.sh $in $output/$sample/$sample.variants.chr$chr.filter.vcf $run_info yes
            done    
        fi
    else
        if [ $variant_type == "SNV" ]
        then
            snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
            for chr in $chrs
            do
                cat $input/$snv_file | awk -v num=chr${chr} '$0 ~ /#/ || $1 == num' > $output/$sample/$sample.variants.chr$chr.filter.vcf
            done
        elif [ $variant_type == "INDEL" ]
        then
            indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)
            for chr in $chrs
            do
                cat $input/$indel_file | awk -v num=chr${chr} '$0 ~ /#/ || $1 == num' >  $output/$sample/$sample.variants.chr$chr.filter.vcf
            done
        fi
    fi
    echo `date`
fi    
            