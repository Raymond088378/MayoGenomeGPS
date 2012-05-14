#!/bin/sh

#	INFO
#	reformat the inputs to chop into chromsome to make it faster fro vraiant module

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
    dbsnp_rsids_snv=$( cat $tool_info | grep -w '^dbSNP_SNV_rsIDs' | cut -d '=' -f2)
    SNV_caller=$( cat $run_info | grep -w '^SNV_CALLER' | cut -d '=' -f2)
    variant_type=`echo "$variant_type" | tr "[a-z]" "[A-Z]"`
    dbsnp_rsids_indel=$( cat $tool_info | grep -w '^dbSNP_INDEL_rsIDs' | cut -d '=' -f2)
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" " " )
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    distance=$( cat $tool_info | grep -w '^SNP_DISTANCE_INDEL' | cut -d '=' -f2 )
	
    if [ $marker -eq 2 ]
    then
        snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
        indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)
        ## format the vcf file to text delimited file
        ## determine if the file is vcf or text input
        for chr in $chrs
        do
            ## extract the file for the chromosome
            ##SNV
            cat $input/$snv_file | awk -v num=chr$chr '$0 ~ /^#/ || ( length($4) == 1 && length($5) == 1 ) || $1 == "num" ' > $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf
            ## INDEL
			cat $input/$indel_file | awk -v num=chr$chr '$0 ~ /^#/ || ( length($4) > 1 || length($5) > 1 ) || $1 == "num" ' > $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
            
			perl $script_path/markSnv_IndelnPos.pl -s $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf -i $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf \
				-n $distance -o $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
			cat $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CLOSE2INDEL="$NF,$9,$10;}' \
				> $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf  
			rm $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf		
        done
    else
        if [ $variant_type == "SNV" ]
        then
            snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
			for chr in $chrs	
            do
                cat $input/$snv_file | awk -v num=chr$chr '$0 ~ /^#/ || ( length($4) == 1 && length($5) == 1 ) || $1 == "num" ' | > $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf
                cat $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CLOSE2INDEL=0",$9,$10;}' \
					> $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
				mv $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf	
				
            done
        elif [ $variant_type == "INDEL" ]
        then
            indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)
            for chr in $chrs		
            do
                cat $input/$indel_file | awk -v num=chr$chr '$0 ~ /^#/ || ( length($4) > 1 || length($5) > 1 ) || $1 == "num" ' > $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
            done
        fi		
    fi
    echo `date`	
fi	
	
	
	
	
	
	
	
		
