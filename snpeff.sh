#!/bin/sh

##	INFO
##	script used to annotate both SNVs and INDELs using snpeff jar script
###############################
#	$1		=		snpeff output directory	
#	$3		=		sample name
#	$2		=		directory for input fie
#	$4		=		run_innfo
################################# 

if [ $# != 4 ];
then
    echo "Usage:<snpeff dir> <input dir><sample><run_info> ";
else
    set -x
    echo `date`
    snpeff=$1
    input=$2
    sample=$3
    run_info=$4
    #SGE_TASK_ID=1
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    snpeff_path=$( cat $tool_info | grep -w '^SNPEFF' | cut -d '=' -f2)
    genome_version=$(cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
	
    if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
    then
        snv_file=$sample.variants.chr$chr.SNV.filter.i.c.vcf
        num_snvs=`cat $input/$snv_file | awk '$0 !~ /^#/' | wc -l`
        if [ $num_snvs -gt 1 ]
        then
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -onlyCoding true -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$snv_file > $snpeff/$sample.chr${chr}.snv.eff
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -onlyCoding true -o vcf -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$snv_file > $snpeff/$sample.chr${chr}.snv.eff.vcf
            perl $script_path/snpeff.pl $snpeff/$sample.chr${chr}.snv.eff > $snpeff/$sample.chr${chr}.snv.eff.fill
            mv $snpeff/$sample.chr${chr}.snv.eff.fill $snpeff/$sample.chr${chr}.snv.eff
        else
            echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tfunctionGVS\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList" > $snpeff/$sample.chr${chr}.snv.eff
            cp $input/$snv_file $snpeff/$sample.chr${chr}.snv.eff.vcf
        fi    
    fi

    if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]	
    then
	indel_file=$sample.variants.chr$chr.INDEL.filter.i.c.vcf
        num_indels=`cat $input/$indel_file | awk '$0 !~ /^#/' | wc -l`
        if [ $num_indels -gt 1 ]
        then
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -onlyCoding true -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$indel_file > $snpeff/$sample.chr${chr}.indel.eff
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -onlyCoding true -o vcf -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$indel_file > $snpeff/$sample.chr${chr}.indel.eff.vcf
            perl $script_path/snpeff.pl $snpeff/$sample.chr${chr}.indel.eff > $snpeff/$sample.chr${chr}.indel.eff.fill
            mv $snpeff/$sample.chr${chr}.indel.eff.fill $snpeff/$sample.chr${chr}.indel.eff
        else
            echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tfunctionGVS\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList" > $snpeff/$sample.chr${chr}.indel.eff
            cp $input/$indel_file $snpeff/$sample.chr${chr}.indel.eff.vcf
        fi
    fi    
    echo `date`
fi	
		
    
