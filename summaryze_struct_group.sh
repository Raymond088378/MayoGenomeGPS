#!/bin/sh
if [ $# != 3 ]
then	
	echo "Usage: <group name> <base dir> <path to run_info file> ";
	exit 1
else
	set -x
	echo `date`	
	group=$1
	basedir=$2
	run_info=$3

	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" )
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2 | tr ":" "\n" )
	samples=$( cat $sample_info | grep -w "^$group" | cut -d '=' -f2)
	output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	ref_genome=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	min_read=$( cat $tool_info | grep -w '^STRUCT_MIN_SUPPORT' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
	ref_genome=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	let num_tumor=`echo $samples|tr " " "\n"|wc -l`-1
	normal_sample=`echo $samples| tr " " "\n" | head -n 1 `
	tumor_list=`echo $samples | tr " " "\n" | tail -$num_tumor`

	#Summaryzing CNVs
	output=$basedir/Reports_per_Sample
	mkdir -p $output/SV
	
	for sample in $tumor_list
	do
		inputargs=""
		inputargs_filter=""
		input=""
	
		for chr in $chrs
		do
			inputfile=$basedir/cnv/$group/$sample.$chr.del.vcf
			input=$basedir/cnv/$group/$sample.$chr.del.bed
			
			## deltion files
			if [ ! -f $inputfile ]
			then
				touch $inputfile.fix.log
				$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
				$script_path/wait.sh $inputfile.fix.log 
				inputargs="-V $inputfile "$inputargs  
				cat $input >> $basedir/cnv/$group/$sample.cnv.bed
				rm $input
			else
				inputargs="-V $inputfile "$inputargs  
				cat $input >> $basedir/cnv/$group/$sample.cnv.bed
				rm $input
			fi	
			inputfile=$basedir/cnv/$group/$sample.$chr.dup.vcf 
			input=$basedir/cnv/$group/$sample.$chr.dup.bed
			if [ ! -f $inputfile ]
			then	
				touch $inputfile.fix.log
				$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
				$script_path/wait.sh $inputfile.fix.log 
				inputargs="-V $inputfile "$inputargs  
				cat $input >> $basedir/cnv/$group/$sample.cnv.bed
				rm $input
			else
				inputargs="-V $inputfile "$inputargs  
				cat $input >> $basedir/cnv/$group/$sample.cnv.bed
				rm $input
			fi
			inputfile=$basedir/cnv/$group/$sample.$chr.filter.del.vcf
			input=$basedir/cnv/$group/$sample.$chr.filter.del.bed
			if [ ! -f $inputfile ]
			then	
				touch $inputfile.fix.log
				$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
				$script_path/wait.sh $inputfile.fix.log 
				cat $input >> $basedir/cnv/$group/$sample.cnv.filter.bed
				rm $input
				inputargs_filter="-V $inputfile "$inputargs_filter  
			else
				cat $input >> $basedir/cnv/$group/$sample.cnv.filter.bed
				rm $input
				inputargs_filter="-V $inputfile "$inputargs_filter  
			fi
			inputfile=$basedir/cnv/$group/$sample.$chr.filter.dup.vcf
			input=$basedir/cnv/$group/$sample.$chr.filter.dup.bed
			if [ ! -f $inputfile ]
			then
				touch $inputfile.fix.log
				$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
				$script_path/wait.sh $inputfile.fix.log 
				cat $input >> $basedir/cnv/$group/$sample.cnv.filter.bed
				rm $input
				inputargs_filter="-V $inputfile "$inputargs_filter 
			else
				cat $input >> $basedir/cnv/$group/$sample.cnv.filter.bed
				rm $input
				inputargs_filter="-V $inputfile "$inputargs_filter 
			fi
		done
		
		$script_path/combinevcf.sh "$inputargs" $output/SV/$group.$sample.cnv.vcf $run_info yes
		$script_path/combinevcf.sh "$inputargs_filter" $output/SV/$group.$sample.cnv.filter.vcf $run_info yes
	done

	#Summaryzing Breakdancer
	
	for sample in $samples
	do
		inputargs=""
		for chr in $chrs
		do
			inputfile=$basedir/struct/break/$group/$sample/$sample.$chr.break
			input_vcf=$basedir/struct/break/$group/$sample/$sample.$chr.break.vcf
			if [ ! -f $inputfile ]
			then
				touch $inputfile.fix.log
				$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
				$script_path/wait.sh $inputfile.fix.log 
				cat $inputfile >> $basedir/struct/$group.$sample.break  
				cat $input_vcf | awk '$0 ~ /^#/' > $basedir/struct/$group.$sample.break.header
	            cat $input_vcf | awk '$0 !~ /^#/' >> $output/SV/$group.$sample.break.vcf
	            rm $inputfile $input_vcf 
			fi
			cat $inputfile >> $basedir/struct/$group.$sample.break  
			cat $input_vcf | awk '$0 ~ /^#/' > $basedir/struct/$group.$sample.break.header
            cat $input_vcf | awk '$0 !~ /^#/' >> $output/SV/$group.$sample.break.vcf
            rm $inputfile $input_vcf 
		done
        input_vcf=$basedir/struct/break/$group/$sample/$sample.inter.break.vcf 
        inputfile=$basedir/struct/break/$group/$sample/$sample.inter.break
        cat $inputfile >> $basedir/struct/$group.$sample.break
        rm $inputfile
        cat $input_vcf | awk '$0 !~ /^#/' >> $output/SV/$group.$sample.break.vcf
        rm $input_vcf 
        cat $basedir/struct/$group.$sample.break.header $output/SV/$group.$sample.break.vcf > $output/SV/$group.$sample.break.vcf.temp
        mv $output/SV/$group.$sample.break.vcf.temp $output/SV/$group.$sample.break.vcf
        rm $basedir/struct/$group.$sample.break.header
	done
	
	### somatic variants subtract normal and tumor to return only tumor variants
	let num_tumor=`echo $samples|tr " " "\n"|wc -l`-1
	normal=`echo $samples| tr " " "\n" | head -n 1 `
	tumor_list=`echo $samples | tr " " "\n" | tail -$num_tumor`
	for tumor in $tumor_list
	do
		perl $script_path/subtract_break.pl $basedir/struct/$group.$tumor.break $basedir/struct/$group.$normal.break > $basedir/struct/$group.$tumor.somatic.break
		perl $script_path/Breakdancer2VCF.pl -i $basedir/struct/$group.$tumor.somatic.break -f $ref_genome -o $output/SV/$group.$tumor.somatic.break.vcf -s $tumor -t $samtools
		perl $script_path/vcfsort.pl ${ref_genome}.fai $output/SV/$group.$tumor.somatic.break.vcf > $output/SV/$group.$tumor.somatic.break.vcf.sort
		mv $output/SV/$group.$tumor.somatic.break.vcf.sort $output/SV/$group.$tumor.somatic.break.vcf
		if [ ! -s $output/SV/$group.$tumor.somatic.break.vcf.fail ]
		then
			rm $output/SV/$group.$tumor.somatic.break.vcf.fail
		fi
	done	
	
	#Summaryzing Crest
	for tumor in $tumor_list
	do
		for chr in $chrs
		do
			inputfile=$basedir/struct/crest/$group/$tumor.$chr.predSV.txt
			inputvcf=$basedir/struct/crest/$group/$tumor.$chr.raw.vcf
			if [ ! -s $inputvcf ]
			then
				touch $inputvcf.fix.log
				$script_path/email.sh $inputvcf "not exist" $JOB_NAME $JOB_ID $run_info
				$script_path/wait.sh $inputvcf.fix.log 
			fi			
			inputfile_filter=$basedir/struct/crest/$group/$tumor.$chr.filter.predSV.txt
			inputvcf_filter=$basedir/struct/crest/$group/$tumor.$chr.filter.vcf
			if [ ! -s $inputvcf_filter ]
			then
				touch $inputvcf_filter.fix.log
				$script_path/email.sh $inputvcf_filter "not exist" $JOB_NAME $JOB_ID $run_info
				$script_path/wait.sh $inputvcf_filter.fix.log
			fi	
			cat $inputfile >> $basedir/struct/$group.$tumor.somatic.raw.crest
			rm $inputfile
			cat $inputfile_filter >> $basedir/struct/$group.$tumor.somatic.filter.crest
			rm $inputfile_filter
			cat $inputvcf | awk '$0 ~ /^#/' > $basedir/struct/$group.$tumor.somatic.header
			cat $inputvcf | awk '$0 !~ /^#/' >> $output/SV/$group.$tumor.somatic.raw.crest.vcf
			cat $inputvcf_filter | awk '$0 !~ /^#/' >> $output/SV/$group.$tumor.somatic.filter.crest.vcf
			rm $inputvcf $inputvcf_filter
		done
		cat $basedir/struct/$group.$tumor.somatic.header $output/SV/$group.$tumor.somatic.raw.crest.vcf > $output/SV/$group.$tumor.somatic.raw.crest.vcf.temp
		mv $output/SV/$group.$tumor.somatic.raw.crest.vcf.temp $output/SV/$group.$tumor.somatic.raw.crest.vcf
		cat $basedir/struct/$group.$tumor.somatic.header $output/SV/$group.$tumor.somatic.filter.crest.vcf > $output/SV/$group.$tumor.somatic.filter.crest.vcf.temp
		mv $output/SV/$group.$tumor.somatic.filter.crest.vcf.temp $output/SV/$group.$tumor.somatic.filter.crest.vcf
		rm $basedir/struct/$group.$tumor.somatic.header
	done	
	echo `date`
fi	




