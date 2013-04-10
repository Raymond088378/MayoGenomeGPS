#!/bin/sh
java=$1
snpeff=$2
GenomeBuild=$3
output=$4
ff=$5
gatk=$6
ref=$7
vcftools=$8
script_path=$9
sample=${10}
	
	
	if [ `cat $output/$ff.SNV.vcf | awk '$0 !~ /^#/' | head -100 | wc -l`  -gt 0 ]
	then
		$java/java -Xmx10g -Xms512m -Djava.io.tmpdir=$output -jar $snpeff/snpEff.jar eff -o vcf -chr chr -noStats -noLog -c $snpeff/snpEff.config $GenomeBuild $output/$ff.SNV.vcf | awk '{if ($0 ~ /##SnpEffVersion/) print "##SnpEffVersion=\"3.0c (build 2012-07-30), by Pablo Cingolani\""; else print $0;}' >  $output/$ff.SNV.vcf.eff.vcf
		cat $output/$ff.SNV.vcf.eff.vcf | $script_path/snpeff.pl > $output/$sample.SNV.eff
		
		$java/java -Xmx10g -Xms512m -Djava.io.tmpdir=$output -jar $gatk/GenomeAnalysisTK.jar \
		-T VariantAnnotator \
		-et NO_ET \
		-K $gatk/Hossain.Asif_mayo.edu.key \
		-R $ref \
		-A SnpEff \
		--variant $output/$ff.SNV.vcf \
		--snpEffFile $output/$ff.SNV.vcf.eff.vcf \
		-L $output/$ff.SNV.vcf \
		-o $output/$ff.SNV.vcf.annot.vcf >> $output/log 2>&1
	
		$vcftools/bin/vcftools --vcf $output/$ff.SNV.vcf.annot.vcf \
		--get-INFO SNPEFF_AMINO_ACID_CHANGE \
		--get-INFO SNPEFF_EFFECT \
		--get-INFO SNPEFF_EXON_ID \
		--get-INFO SNPEFF_FUNCTIONAL_CLASS \
		--get-INFO SNPEFF_GENE_BIOTYPE \
		--get-INFO SNPEFF_GENE_NAME \
		--get-INFO SNPEFF_IMPACT \
		--get-INFO SNPEFF_TRANSCRIPT_ID \
		--out $output/$ff.SNV.vcf.annot >> $output/log 2>&1
            
		perl $script_path/parse_snpeffect.pl $output/$ff.SNV.vcf.annot.INFO > $output/$sample.SNV.filtered.eff
		rm $output/$ff.SNV.vcf.annot.INFO $output/$ff.SNV.vcf.annot.log
		rm $output/$ff.SNV.vcf.annot.vcf $output/$ff.SNV.vcf.annot.vcf.idx $output/$ff.SNV.vcf.annot.vcf.vcfidx
		rm  $output/$ff.SNV.vcf.idx $output/$ff.SNV.vcf.eff.vcf $output/$ff.SNV.vcf.eff.vcf.idx 
	else
		echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList" > $output/$sample.SNV.eff
		echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\tFunctionalClass\tFunctionalImpact\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList\n" > $output/$sample.SNV.filtered.eff
	fi
	echo " Filtering SNPEFF output from multiple transcript to most impacting transcript "
	## INDELs 
	if [ `cat $output/$ff.INDEL.vcf| awk '$0 !~ /^#/' | head -100 | wc -l` -ge 1 ]
	then
		$java/java -Xmx10g -Xms512m -Djava.io.tmpdir=$output -jar $snpeff/snpEff.jar eff -o vcf -chr chr -noStats -noLog -c $snpeff/snpEff.config $GenomeBuild $output/$ff.INDEL.vcf | awk '{if ($0 ~ /##SnpEffVersion/) print "##SnpEffVersion=\"3.0c (build 2012-07-30), by Pablo Cingolani\""; else print $0;}' > $output/$ff.INDEL.vcf.eff.vcf

		cat $output/$ff.INDEL.vcf.eff.vcf | $script_path/snpeff.pl > $output/$sample.INDEL.eff	
		
		### use GATK to filter the multiple transcript
		$java/java -Xmx10g -Xms512m -Djava.io.tmpdir=$output  -jar $gatk/GenomeAnalysisTK.jar \
		-T VariantAnnotator \
		-et NO_ET \
		-K $gatk/Hossain.Asif_mayo.edu.key \
		-R $ref \
		-A SnpEff \
		--variant $output/$ff.INDEL.vcf \
		--snpEffFile $output/$ff.INDEL.vcf.eff.vcf \
		-L $output/$ff.INDEL.vcf \
		-o $output/$ff.INDEL.vcf.annot.vcf >> $output/log 2>&1

		$vcftools/bin/vcftools --vcf $output/$ff.INDEL.vcf.annot.vcf \
		--get-INFO SNPEFF_AMINO_ACID_CHANGE \
		--get-INFO SNPEFF_EFFECT \
		--get-INFO SNPEFF_EXON_ID \
		--get-INFO SNPEFF_FUNCTIONAL_CLASS \
		--get-INFO SNPEFF_GENE_BIOTYPE \
		--get-INFO SNPEFF_GENE_NAME \
		--get-INFO SNPEFF_IMPACT \
		--get-INFO SNPEFF_TRANSCRIPT_ID \
		--out $output/$ff.INDEL.vcf.annot >> $output/log 2>&1

		perl $script_path/parse_snpeffect.pl $output/$ff.INDEL.vcf.annot.INFO > $output/$sample.INDEL.filtered.eff
		rm $output/$ff.INDEL.vcf.annot.INFO $output/$ff.INDEL.vcf.annot.log
		rm $output/$ff.INDEL.vcf.annot.vcf $output/$ff.INDEL.vcf.annot.vcf.idx $output/$ff.INDEL.vcf.annot.vcf.vcfidx
		rm  $output/$ff.INDEL.vcf.idx $output/$ff.INDEL.vcf.eff.vcf $output/$ff.INDEL.vcf.eff.vcf.idx 
	else
		echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList" > $output/$sample.INDEL.eff
	    echo -e "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\tFunctionalClass\tFunctionalImpact\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList\n" >   $output/$sample.INDEL.filtered.eff
   fi
	
	echo " SNPEFF is done "
