#!/bin/bash

##	INFO
##	script used to annotate both SNVs and INDELs using snpeff jar script
###############################
#	$1		=		snpeff output directory	
#	$3		=		sample name
#	$2		=		directory for input fie
#	$4		=		run_innfo
################################# 

if [ $# -le 4 ];
then
    echo -e "script to run snpeff on a vcf file in a folder\nUsage:<snpeff dir> <input dir><sample><run_info><somatic/germline><SGE_TASK_ID(optional)>";
else
    set -x
    echo `date`
    snpeff=$1
    input=$2
    sample=$3
    run_info=$4
    type=`echo $5 | tr "[A-Z]" "[a-z]"`	
	if [ $type == "somatic" ]
	then
		prefix="TUMOR"
		sam=$prefix.$sample
	else
		sam=$sample
	fi	
	if [ $6 ]
    then
        SGE_TASK_ID=$6
    fi	
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
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

    if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
    then
        snv_file=$sam.variants.chr$chr.SNV.filter.i.c.vcf
        if [ ! -s $input/$snv_file ]
        then
            $script_path/errorlog.sh $input/$snv_file snpeff.sh ERROR "not found"
			exit 1;
        fi
        num_snvs=`cat $input/$snv_file | awk '$0 !~ /^#/' | wc -l`
		$script_path/filesize.sh snpeff $sam $input $snv_file $JOB_ID $run_info
		
        if [ $num_snvs -ge 1 ]
        then
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -o txt -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$snv_file | $script_path/snpeff.pl > $snpeff/$sam.chr${chr}.snv.eff
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -o vcf -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$snv_file | awk '{if ($0 ~ /##SnpEffVersion/) print "##SnpEffVersion=\"3.0c (build 2012-07-30), by Pablo Cingolani\""; else print $0;}' > $snpeff/$sam.chr${chr}.snv.eff.vcf
			
            ### use GATK to filter the multiple transcript
            $java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar -T VariantAnnotator \
            -et NO_ET -K $gatk/Hossain.Asif_mayo.edu.key \
            -R $ref -A SnpEff --variant $input/$snv_file \
            --snpEffFile $snpeff/$sam.chr${chr}.snv.eff.vcf \
            -L $input/$snv_file -o $snpeff/$snv_file.annotate.vcf
            
            $vcftools/bin/vcftools --vcf $snpeff/$snv_file.annotate.vcf \
            --get-INFO SNPEFF_AMINO_ACID_CHANGE --get-INFO SNPEFF_EFFECT \
            --get-INFO SNPEFF_EXON_ID --get-INFO SNPEFF_FUNCTIONAL_CLASS --get-INFO SNPEFF_GENE_BIOTYPE \
            --get-INFO SNPEFF_GENE_NAME --get-INFO SNPEFF_IMPACT --get-INFO SNPEFF_TRANSCRIPT_ID \
            --out $snpeff/$snv_file.annotate
            
            $script_path/parse_snpeffect.pl $snpeff/$snv_file.annotate.INFO > $snpeff/$sam.chr${chr}.snv.filtered.eff
            rm $snpeff/$snv_file.annotate.INFO $snpeff/$snv_file.annotate.log $snpeff/$sam.chr${chr}.snv.eff.vcf
            rm $snpeff/$sam.chr${chr}.snv.eff.vcf.idx $snpeff/$snv_file.annotate.vcf.vcfidx $snpeff/$snv_file.annotate.vcf.idx $snpeff/$snv_file.annotate.vcf
        else
            echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList" > $snpeff/$sam.chr${chr}.snv.eff
			echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\tFunctionalClass\tFunctionalImpact\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList\n" > $snpeff/$sam.chr${chr}.snv.filtered.eff
        fi  
		if [[ ! -s $snpeff/$sam.chr${chr}.snv.eff || ! -s $snpeff/$sam.chr${chr}.snv.filtered.eff ]]
		then
			$script_path/errorlog.sh "$snpeff/$sam.chr${chr}.snv.eff $snpeff/$sam.chr${chr}.snv.filtered.eff" snpeff.sh ERROR "failed to create" 
			exit 1;
		else
			$script_path/filesize.sh snpeff.out $sam $snpeff $sam.chr${chr}.snv.eff $JOB_ID $run_info
		fi
	fi

    if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]	
    then
	indel_file=$sam.variants.chr$chr.INDEL.filter.i.c.vcf
        if [ ! -s $input/$indel_file ]
        then
            $script_path/errorlog.sh $input/$indel_file snpeff.sh ERROR "not found"
			exit 1;
        fi
        num_indels=`cat $input/$indel_file | awk '$0 !~ /^#/' | wc -l`
        if [ $num_indels -ge 1 ]
        then
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -o txt -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$indel_file  | $script_path/snpeff.pl > $snpeff/$sam.chr${chr}.indel.eff
            $java/java -Xmx2g -Xms512m -jar $snpeff_path/snpEff.jar eff -o vcf -chr chr -noStats -noLog \
                -c $snpeff_path/snpEff.config $genome_version $input/$indel_file | awk '{if ($0 ~ /##SnpEffVersion/) print "##SnpEffVersion=\"3.0c (build 2012-07-30), by Pablo Cingolani\""; else print $0;}' > $snpeff/$sam.chr${chr}.indel.eff.vcf
            
            ### use GATK to filter the multiple transcript
            $java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar -T VariantAnnotator \
            -et NO_ET -K $gatk/Hossain.Asif_mayo.edu.key \
            -R $ref -A SnpEff --variant $input/$indel_file \
            --snpEffFile $snpeff/$sam.chr${chr}.indel.eff.vcf -L $input/$indel_file \
            -o $snpeff/$indel_file.annotate.vcf
            
            $vcftools/bin/vcftools --vcf $snpeff/$indel_file.annotate.vcf \
            --get-INFO SNPEFF_AMINO_ACID_CHANGE --get-INFO SNPEFF_EFFECT \
            --get-INFO SNPEFF_EXON_ID --get-INFO SNPEFF_FUNCTIONAL_CLASS \
            --get-INFO SNPEFF_GENE_BIOTYPE --get-INFO SNPEFF_GENE_NAME \
            --get-INFO SNPEFF_IMPACT --get-INFO SNPEFF_TRANSCRIPT_ID \
            --out $snpeff/$indel_file.annotate
            
            $script_path/parse_snpeffect.pl $snpeff/$indel_file.annotate.INFO > $snpeff/$sam.chr${chr}.indel.filtered.eff
            rm $snpeff/$indel_file.annotate.INFO $snpeff/$indel_file.annotate.log $snpeff/$sam.chr${chr}.indel.eff.vcf
            rm $snpeff/$sam.chr${chr}.indel.eff.vcf.idx $snpeff/$indel_file.annotate.vcf.vcfidx $snpeff/$indel_file.annotate.vcf.idx $snpeff/$indel_file.annotate.vcf
		else
			echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList" > $snpeff/$sam.chr${chr}.indel.eff
			echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\tFunctionalClass\tFunctionalImpact\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList\n" >  $snpeff/$sam.chr${chr}.indel.filtered.eff
        fi
		if [[ ! -s $snpeff/$sam.chr${chr}.indel.eff || ! -s $snpeff/$sam.chr${chr}.indel.filtered.eff ]]
		then
			$script_path/errorlog.sh "$snpeff/$sam.chr${chr}.indel.eff $snpeff/$sam.chr${chr}.indel.filtered.eff" snpeff.sh ERROR "failed to create" 
			exit 1;
		else
			$script_path/filesize.sh snpeff.out $sam $snpeff $sam.chr${chr}.indel.eff $JOB_ID $run_info
		fi
    fi    
    echo `date`
fi	
		
    
