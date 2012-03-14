#!/bin/sh

if [ $# != 2 ]
then
	echo "Usage : </path/to/output folder> </path/to/run info file>"
else
	set -x
	echo `date`
	output_dir=$1
	run_info=$2
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	samples=$(cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2| tr ":" " ")
	chrs=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2| tr ":" " ")	
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	
	var_dir=$output_dir/variants
	inputargs=""
	input_index=""
	for sample in $samples
	do
		input=$var_dir/$sample
		for chr in $chrs
		do
			inputfile=$input/$sample.variants.chr$chr.raw.all.vcf
			input_indexfile=$input/$sample.variants.chr$chr.raw.all.vcf.idx
			if [ -s $inputfile ]
			then
				inputargs="-V $inputfile "$inputargs
				input_index=" $input_indexfile"$input_index 
			else
				echo "WARNING : $inputfile not there"
			fi			
		done
	done

	### merge all the variants using GATK
	$java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
	-R $ref \
	-et NO_ET \
	-T CombineVariants \
	$inputargs \
	-o $var_dir/raw.SNV.vcf
	
	gzip $var_dir/raw.SNV.vcf
	if [ -s $var_dir/raw.SNV.vcf.gz ]
	then
		file=`echo $inputargs | sed -e '/-V/s///g'`
		rm $file
		rm $input_index
	fi	
	echo `date`
fi	