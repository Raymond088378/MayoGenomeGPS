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
        type=`cat $input/$snv_file | head -1 | awk '{if ($0 ~ /^##/) print "vcf";else print "txt"}'`
        if [ $type == "txt" ]
        then
            perl $script_path/convert_txt_vcf.pl $input/$snv_file $sample > $output/$sample.SNV.vcf
            perl $script_path/convert_txt_vcf.pl $input/$indel_file $sample > $output/$sample.INDEL.vcf
        else
            ln -s $input/$snv_file $output/$sample.SNV.vcf
            ln -s $input/$indel_file $output/$sample.INDEL.vcf 
        fi
        for chr in $chrs
        do
            perl $script_path/vcf_to_variant_vcf.pl -i $output/$sample.SNV.vcf -v $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf -l $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -c chr$chr -s $sample
            rm $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
            perl $script_path/vcf_to_variant_vcf.pl -i $output/$sample.INDEL.vcf -v $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf.temp -l $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -c chr$chr -s $sample
            rm $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf.temp
            cat $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1",$9,$10;}' > $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
            perl $script_path/add_format_field_vcf.pl $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp INDEL > $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
            rm $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
            perl $script_path/markSnv_IndelnPos.pl -s $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf -i $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -n $distance -o $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
            cat $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1;CLOSE2INDEL="$NF,$9,$10;}' > $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf  
			perl $script_path/add_format_field_vcf.pl $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf SNV > $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf 
            mv $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf		
        done
        rm $output/$sample.SNV.vcf $output/$sample.INDEL.vcf 
    else
        if [ $variant_type == "SNV" ]
        then
            snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
            type=`cat $input/$snv_file | head -1 | awk '{if ($0 ~ /^##/) print "vcf";else print "txt"}'`
            if [ $type == "txt" ]
            then
                perl $script_path/convert_txt_vcf.pl $input/$snv_file $sample > $output/$sample.SNV.vcf
            else
                ln -s $input/$snv_file $output/$sample.SNV.vcf
            fi    
            for chr in $chrs	
            do
                perl $script_path/vcf_to_variant_vcf.pl -i $output/$sample.SNV.vcf -v $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf -l $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -c chr$chr -s $sample
                rm $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
                cat $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1;CLOSE2INDEL=0",$9,$10;}' > $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
				perl $script_path/add_format_field_vcf.pl $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf SNV > $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf
                rm $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf			
            done
            rm $output/$sample.SNV.vcf
        elif [ $variant_type == "INDEL" ]
        then
            indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)
            type=`cat $input/$indel_file | head -1 | awk '{if ($0 ~ /^##/) print "vcf";else print "txt"}'`
            if [ $type == "txt" ]
            then
                perl $script_path/convert_txt_vcf.pl $input/$indel_file $sample > $output/$sample.INDEL.vcf 
            else
                ln -s $input/$indel_file $output/$sample.INDEL.vcf 
            fi    
            for chr in $chrs		
            do
                perl $script_path/vcf_to_variant_vcf.pl -i $output/$sample.INDEL.vcf -v $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf -l $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -c chr$chr -s $sample
                rm $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf
                cat $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1",$9,$10;}' > $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
				perl $script_path/add_format_field_vcf.pl $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp INDEL > $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf  
                rm $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
            done
            rm $output/$sample.INDEL.vcf 
        fi		
    fi
    echo `date`	
fi	
	
	
	
	
	
	
	
		
