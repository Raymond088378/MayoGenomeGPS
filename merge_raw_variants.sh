#!/bin/sh

if [ $# -le 1 ]
then
	echo "Usage : </path/to/output folder> </path/to/run info file>"
else
	set -x
	echo `date`
	output_dir=$1
	run_info=$2
        if [ $3 ]
        then
            SGE_TASK_ID=$3
        fi    
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	samples=$(cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2| tr ":" " ")
	chrs=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2| tr ":" " ")	
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	vcftools=$( cat $tool_info | grep -w '^VCFTOOLS' | cut -d '=' -f2)
	perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
	tabix=$( cat $tool_info | grep -w '^TABIX' | cut -d '=' -f2)
	chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
	export PERL5LIB=$perllib
	PATH=$tabix/:$PATH
	var_dir=$output_dir/variants
	inputargs=""
	input_index=""
	for sample in $samples
	do
		input=$var_dir/$sample
		inputfile=$input/$sample.variants.chr$chr.raw.all.vcf.gz
		input_indexfile=$input/$sample.variants.chr$chr.raw.all.vcf.gz.tbi
			
		if [ ! -s $input_indexfile ]
		then
			$tabix/tabix -p vcf $input/$sample.variants.chr$chr.raw.all.vcf.gz
		fi	
		if [ -s $inputfile ]
		then
			inputargs="$inputfile "$inputargs
			input_index=" $input_indexfile"$input_index 
		else
			$script_path/errorlog.sh $inputfile merge_raw_variants.sh WARNING "does not exist"
		fi			
	done
	
	$vcftools/bin/vcf-merge -s $inputargs > $var_dir/raw.chr$chr.vcf
	
	### raw SNV
	perl $script_path/vcf_to_variant_vcf.pl -i $var_dir/raw.chr$chr.vcf -v $var_dir/raw.chr$chr.SNV.vcf -l $var_dir/raw.chr$chr.INDEL.vcf
	$tabix/bgzip $var_dir/raw.chr$chr.SNV.vcf  
	$tabix/tabix -p vcf $var_dir/raw.chr$chr.SNV.vcf.gz	
	### raw INDEL
	$tabix/bgzip $var_dir/raw.chr$chr.INDEL.vcf
	$tabix/tabix -p vcf $var_dir/raw.chr$chr.INDEL.vcf.gz  
	if [[ -s $var_dir/raw.chr$chr.INDEL.vcf.gz && -s $var_dir/raw.chr$chr.SNV.vcf.gz ]]
    then
        rm $var_dir/raw.chr$chr.vcf 
		rm $inputargs
		rm $input_index
	fi
	echo `date`
fi	