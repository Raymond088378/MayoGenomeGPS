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
	
	if [ $marker -eq 2 ]
	then
		snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
		indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)
		## format the vcf file to text delimited file
		perl $script_path/parse.vcf.INDEL.pl -i $input/$indel_file -o $output/$sample.indels -s $sample
		perl $script_path/parse.vcf.SNV.pl -i $input/$snv_file -o $output/$sample.snvs -s $sample
		for chr in $chrs
		do
			## extract the file for the chromosome
			##SNV
			cat $output/$sample.snvs | grep -w chr${chr} | awk '{print $0"\t1"}' | sort -T $output -n -k 2,12n > $output/$sample.chr${chr}.raw.snvs.bed.i.ToMerge
			`dos2unix $output/$sample.chr${chr}.raw.snvs`
			## INDEL
			cat $output/$sample.indels | grep -w chr${chr} | awk '{print $0"\t1"}' | sort -T $output -n -k 2,12n > $output/$sample.chr${chr}.raw.indels.bed.i.ToMerge
			`dos2unix $output/$sample.chr${chr}.raw.indels`
			
			perl $script_path/markSnv_IndelnPos.pl -s $output/$sample.chr${chr}.raw.snvs.bed.i.ToMerge -i $output/$sample.chr${chr}.raw.indels.bed.i.ToMerge -n 10 -p 2 -o $output/$sample.chr${chr}.raw.snvs.bed.i.ToMerge.pos
			mv $output/$sample.chr${chr}.raw.snvs.bed.i.ToMerge.pos $output/$sample.chr${chr}.raw.snvs.bed.i.ToMerge
		done
		rm $output/$sample.indels
		rm $output/$sample.snvs
	else
		if [ $variant_type == "SNV" ]
		then
			snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
			perl $script_path/parse.vcf.SNV.pl -i $input/$snv_file -o $output/$sample.snvs -s $sample
			for chr in $chrs	
			do
				cat $output/$sample.snvs | grep -w chr${chr} | awk '{print $0"\t1\t0"}' | sort -T $output -n -k 2,12n > $output/$sample.chr${chr}.raw.snvs.bed.i.ToMerge
				`dos2unix $output/$sample.chr${chr}.raw.snvs.bed.i.ToMerge`	
			done
			rm $output/$sample.snvs
		elif [ $variant_type == "INDEL" ]
		then
			indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)
			perl $script_path/parse.vcf.INDEL.pl -i $input/$indel_file -o $output/$sample.indels -s $sample
			for chr in $chrs		
			do
				cat $output/$sample.indels | grep -w chr${chr} | awk '{print $0"\t1"}' | sort -T $output -n -k 2,12n > $output/$sample.chr${chr}.raw.indels.bed.i.ToMerge
				`dos2unix $output/$sample.chr${chr}.raw.indels`
			done
			rm $output/$sample.indels
		fi		
	fi
	echo `date`	
fi	
	
	
	
	
	
	
	
		
