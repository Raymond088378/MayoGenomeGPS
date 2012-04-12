#!/bin/sh

if [ $# != 5 ]
then
	echo -e "Usage: to annotate vcf files \n </path/to/tumor bam> </path/to/normal bam> <input vcf file><chromosome><run info file>"
else
	set -x
	echo `date`
	
	tumor_bam=$1
	normal_bam=$2
	vcf=$3
	chr=$4
	run_info=$5
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)


	# ## annotate SNVs
	$java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
	-R $ref \
	-et NO_ET \
	-T VariantAnnotator \
	-I $tumor_bam \
	-I $normal_bam \
	-V $vcf \
	--dbsnp $dbSNP \
	-L chr$chr \
	-A QualByDepth -A HaplotypeScore -A DepthOfCoverage -A MappingQualityZero -A RMSMappingQuality -A FisherStrand \
	--out $vcf.temp	

	if [ -s $vcf.temp ]
	then
		mv $vcf.temp $vcf
	else		
		echo "ERROR : variants.sh SNV VariantAnnotator failed, file: $vcf.temp"
		exit 1
	fi

	if [ -s $vcf.temp.idx ]
	then
		mv $vcf.temp.idx $vcf.idx
	else	
		echo "ERROR: variants.sh SNV VariantAnnotator did not generated index, file: $vcf.idx" 
		exit 1
	fi
	echo `date`
fi	