#!/bin/bash

if [ $# != 4 ]
then
	echo "Usage : </path/to/input folder> <sample name> </path/to/output folder> </path/to/run info file>"
else
	set -x
	echo `date`
	input_dir=$1
	sample=$2
	output_dir=$3
	run_info=$4
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	chrs=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2| tr ":" " ")	
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	vcftools=$( cat $tool_info | grep -w '^VCFTOOLS' | cut -d '=' -f2)
	perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
	tabix=$( cat $tool_info | grep -w '^TABIX' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	export PERL5LIB=$PERL5LIB:$perllib
	PATH=$tabix/:$PATH
        
	inputargs=""
	indexes=""
	for chr in $chrs
	do
		inputfile=$input_dir/$sample.variants.chr${chr}.raw.all.vcf.gz
		$tabix/bgzip -d -c $inputfile >  $output_dir/$sample.variants.chr${chr}.raw.all.vcf
		rm $inputfile
		perl $script_path/vcf_to_variant_vcf.pl -i $output_dir/$sample.variants.chr${chr}.raw.all.vcf -v $output_dir/$sample.variants.chr${chr}.raw.all.snv.vcf -t snv
		rm $output_dir/$sample.variants.chr${chr}.raw.all.vcf 
		inputfile=$output_dir/$sample.variants.chr${chr}.raw.all.snv.vcf
		$tabix/bgzip $inputfile
		$tabix/tabix -p vcf $inputfile.gz
		inputfile_i=$inputfile.gz.tbi
		indexes=$indexes" $inputfile_i"
		inputfile=$output_dir/$sample.variants.chr${chr}.raw.all.snv.vcf.gz
		inputargs=$inputargs" $inputfile"
	done
	
	$vcftools/bin/vcf-concat $inputargs | $script_path/filter.vcf.mrna.pl | $tabix/bgzip > $output_dir/$sample.vcf.gz
	rm $inputargs
	rm $indexes
	### filter the variants using filters # of alternate reads >=2 and ratio of alt/total >= 0.1
    $tabix/tabix -p vcf $output_dir/$sample.vcf.gz
	$vcftools/bin/vcf-query $output_dir/$sample.vcf.gz -f '%CHROM\t%POS\t%INFO/CAPTURE\t%INFO/ED\t%REF\t%ALT\t[%GT\t%AD\t%DP\t%GQ\t%C2I]\t%INFO/RATIO\n' | sed -e '/\//s///g' | tr "," "\t" > $output_dir/$sample.snvs
	rm $output_dir/$sample.vcf.gz $output_dir/$sample.vcf.gz.tbi
    echo `date`    
fi        