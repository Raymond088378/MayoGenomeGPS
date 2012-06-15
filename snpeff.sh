#!/bin/sh

##	INFO
##	script used to annotate both SNVs and INDELs using snpeff jar script
###############################
#	$1		=		snpeff output directory	
#	$3		=		sample name
#	$2		=		directory for input fie
#	$4		=		run_innfo
################################# 

if [ $# -le 3 ];
then
    echo "Usage:<snpeff dir> <input dir><sample><run_info> ";
else
    set -x
    echo `date`
    snpeff=$1
    input=$2
    sample=$3
    run_info=$4
    if [ $5 ]
	then
		prefix=$5
	fi	
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    snpeff_path=$( cat $tool_info | grep -w '^SNPEFF' | cut -d '=' -f2)
    genome_version=$(cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    vcftools=$( cat $tool_info | grep -w '^VCFTOOLS' | cut -d '=' -f2)
    perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
    export PERL5LIB=$perllib:$PERL5LIB
	
	if [ $5 ]
	then
		sam=$prefix.$sample
	else
		sam=$sample
	fi	
	
    if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
    then
        snv_file=$sam.variants.chr$chr.SNV.filter.i.c.vcf
        num_snvs=`cat $input/$snv_file | awk '$0 !~ /^#/' | wc -l`
        if [ $num_snvs -ge 1 ]
        then
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -onlyCoding true -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$snv_file > $snpeff/$sam.chr${chr}.snv.eff
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -onlyCoding true -o vcf -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$snv_file > $snpeff/$sam.chr${chr}.snv.eff.vcf
			cat $snpeff/$sam.chr${chr}.snv.eff.vcf | awk '{if ($0 ~ /##SnpEffVersion/) print "##SnpEffVersion=\"2.0.5 (build 2012-01-19), by Pablo Cingolani\""; else print $0;}' > $snpeff/$sam.chr${chr}.snv.eff.vcf.tmp
			mv $snpeff/$sam.chr${chr}.snv.eff.vcf.tmp $snpeff/$sam.chr${chr}.snv.eff.vcf	
            perl $script_path/snpeff.pl $snpeff/$sam.chr${chr}.snv.eff > $snpeff/$sam.chr${chr}.snv.eff.fill
            mv $snpeff/$sam.chr${chr}.snv.eff.fill $snpeff/$sam.chr${chr}.snv.eff
			
            ### use GATK to filter the multiple transcript
            $java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
            -T VariantAnnotator \
            -et NO_ET \
            -K $gatk/Hossain.Asif_mayo.edu.key \
            -R $ref \
            -A SnpEff \
            --variant $input/$snv_file \
            --snpEffFile $snpeff/$sam.chr${chr}.snv.eff.vcf \
            -L $input/$snv_file \
            -o $snpeff/$snv_file.annotate.vcf
            
            $vcftools/bin/vcftools --vcf $snpeff/$snv_file.annotate.vcf \
            --get-INFO SNPEFF_AMINO_ACID_CHANGE \
            --get-INFO SNPEFF_EFFECT \
            --get-INFO SNPEFF_EXON_ID \
            --get-INFO SNPEFF_FUNCTIONAL_CLASS \
            --get-INFO SNPEFF_GENE_BIOTYPE \
            --get-INFO SNPEFF_GENE_NAME \
            --get-INFO SNPEFF_IMPACT \
            --get-INFO SNPEFF_TRANSCRIPT_ID \
            --out $snpeff/$snv_file.annotate
            
            perl $script_path/parse_snpeffect.pl $snpeff/$snv_file.annotate.INFO > $snpeff/$sam.chr${chr}.snv.filtered.eff
            rm $snpeff/$snv_file.annotate.INFO
            rm $snpeff/$snv_file.annotate.log
            rm $snpeff/$sam.chr${chr}.snv.eff.vcf
            rm $snpeff/$sam.chr${chr}.snv.eff.vcf.idx
            rm $snpeff/$snv_file.annotate.vcf.vcfidx
            rm $snpeff/$snv_file.annotate.vcf.idx
        else
            echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList" > $snpeff/$sam.chr${chr}.snv.eff
            cp $input/$snv_file $snpeff/$snv_file.annotate.vcf
			echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\tFunctionalClass\tFunctionalImpact\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList\n" > $snpeff/$sam.chr${chr}.snv.filtered.eff
        fi    
    fi

    if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]	
    then
		indel_file=$sam.variants.chr$chr.INDEL.filter.i.c.vcf
        num_indels=`cat $input/$indel_file | awk '$0 !~ /^#/' | wc -l`
        if [ $num_indels -ge 1 ]
        then
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -onlyCoding true -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$indel_file > $snpeff/$sam.chr${chr}.indel.eff
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -onlyCoding true -o vcf -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$indel_file > $snpeff/$sam.chr${chr}.indel.eff.vcf
            perl $script_path/snpeff.pl $snpeff/$sam.chr${chr}.indel.eff > $snpeff/$sam.chr${chr}.indel.eff.fill
            mv $snpeff/$sam.chr${chr}.indel.eff.fill $snpeff/$sam.chr${chr}.indel.eff
			
            cat $snpeff/$sam.chr${chr}.indel.eff.vcf | awk '{if ($0 ~ /##SnpEffVersion/) print "##SnpEffVersion=\"2.0.5 (build 2012-01-19), by Pablo Cingolani\""; else print $0;}' > $snpeff/$sam.chr${chr}.indel.eff.vcf.tmp
            mv $snpeff/$sam.chr${chr}.indel.eff.vcf.tmp $snpeff/$sam.chr${chr}.indel.eff.vcf
            
            ### use GATK to filter the multiple transcript
            $java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
            -T VariantAnnotator \
            -et NO_ET \
            -K $gatk/Hossain.Asif_mayo.edu.key \
            -R $ref \
            -A SnpEff \
            --variant $input/$indel_file \
            --snpEffFile $snpeff/$sam.chr${chr}.indel.eff.vcf \
            -L $input/$indel_file \
            -o $snpeff/$indel_file.annotate.vcf
            
            $vcftools/bin/vcftools --vcf $snpeff/$indel_file.annotate.vcf \
            --get-INFO SNPEFF_AMINO_ACID_CHANGE \
			--get-INFO SNPEFF_EFFECT \
            --get-INFO SNPEFF_EXON_ID \
            --get-INFO SNPEFF_FUNCTIONAL_CLASS \
            --get-INFO SNPEFF_GENE_BIOTYPE \
            --get-INFO SNPEFF_GENE_NAME \
            --get-INFO SNPEFF_IMPACT \
            --get-INFO SNPEFF_TRANSCRIPT_ID \
            --out $snpeff/$indel_file.annotate
            
            perl $script_path/parse_snpeffect.pl $snpeff/$indel_file.annotate.INFO > $snpeff/$sam.chr${chr}.indel.filtered.eff
            rm $snpeff/$indel_file.annotate.INFO
            rm $snpeff/$indel_file.annotate.log
            rm $snpeff/$sam.chr${chr}.indel.eff.vcf
            rm $snpeff/$sam.chr${chr}.indel.eff.vcf.idx
            rm $snpeff/$indel_file.annotate.vcf.vcfidx
            rm $snpeff/$indel_file.annotate.vcf.idx
		else
			echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList" > $snpeff/$sam.chr${chr}.indel.eff
			cp $input/$indel_file $snpeff/$indel_file.annotate.vcf
			echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\tFunctionalClass\tFunctionalImpact\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList\n" >  $snpeff/$sam.chr${chr}.indel.filtered.eff
        fi
    fi    
    echo `date`
fi	
		
    
