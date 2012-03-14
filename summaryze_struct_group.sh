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
	runinfo=$3

	chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" )
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2 | tr ":" "\n" )
	samples=$( cat $sample_info | grep -w "^$group" | cut -d '=' -f2)
	output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	min_read=$( cat $tool_info | grep -w '^STRUCT_MIN_SUPPORT' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
	
	let num_tumor=`echo $samples|tr "\t" "\n"|wc -l`-1
	normal_sample=`echo $samples| tr "\t" "\n" | head -n 1 `
	tumor_list=`echo $sample | tr "\t" "\n" | tail -$num_tumor`

	#Summaryzing CNVs
	output=$basedir/Reports_per_Sample
	mkdir $output/SV
	
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
				echo "ERROR:summaryze_struct. CNV: File $inputfile does not exist "
				exit 1
			fi
			cat $inputfile >> $basedir/cnv/$group/$sample.del.bed  

			inputargs="-V $inputfile "$inputargs  
			cat $input >> $basedir/cnv/$group/$sample.cnv.bed
			rm $input
			
			inputfile=$basedir/cnv/$group/$sample.$chr.dup.vcf 
			input=$basedir/cnv/$group/$sample.$chr.dup.bed
			if [ ! -f $inputfile ]
			then	
				echo "ERROR :summaryze_struct. CNV: File $inputfile does not exist "
				exit 1
			fi
		
			inputargs="-V $inputfile "$inputargs  
			cat $input >> $basedir/cnv/$group/$sample.cnv.bed
			rm $input
			inputfile=$basedir/cnv/$group/$sample.$chr.filter.del.vcf
			input=$basedir/cnv/$group/$sample.$chr.filter.del.bed
			if [ ! -f $inputfile ]
			then	
				echo "ERROR :summaryze_struct. CNV: File $inputfile does not exist "
				exit 1
			fi
			cat $input >> $basedir/cnv/$group/$sample.cnv.filter.bed
			rm $input
			inputargs_filter="-V $inputfile "$inputargs_filter  
			inputfile=$basedir/cnv/$group/$sample.$chr.filter.dup.vcf
			input=$basedir/cnv/$group/$sample.$chr.filter.dup.bed
			if [ ! -f $inputfile ]
			then
				echo "ERROR :summaryze_struct. CNV: File $inputfile does not exist "
				exit 1
			fi
			cat $input >> $basedir/cnv/$group/$sample.cnv.filter.bed
			rm $input
			inputargs_filter="-V $inputfile "$inputargs_filter 
		done
		
		$java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
		-R $ref \
		-et NO_ET \
		-T CombineVariants \
		$inputargs \
		-o $output/SV/$group.$sample.cnv.vcf
		
		$java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
		-R $ref \
		-et NO_ET \
		-T CombineVariants \
		$inputargs_filter \
		-o $output/SV/$group.$sample.cnv.filter.vcf
		
		if [ -s $output/SV/$group.$sample.cnv.vcf ]
		then
			file=`echo $inputargs | sed -e '/-V/s///g'`
			rm $file
		fi

		if [ -s $output/SV/$group.$sample.cnv.filter.vcf ]
		then
			file=`echo $inputargs_filter | sed -e '/-V/s///g'`
			rm $file
		fi    
		rm $basedir/cnv/$group/$sample.*.idx
	done

	#Summaryzing Breakdancer
	
	for sample in $samples
	do
		inputargs=""
		rm $basedir/struct/$group.$sample.break
		for chr in $chrs
		do
			inputfile=$basedir/struct/break/$group/$sample/$sample.$chr.break
			input_vcf=$basedir/struct/break/$group/$sample/$sample.$chr.break.vcf
			if [ ! -f $inputfile ]
			then
				echo "ERROR:summaryze_struct. Breakdancer: File $inputfile does not exist "
				exit 1
			fi
			cat $inputfile >> $basedir/struct/$group.$sample.break  
			inputargs="-V $input_vcf "$inputargs
		done
		### merge all the vcfs
		$java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
		-R $ref \
		-et NO_ET \
		-T CombineVariants \
		$inputargs \
		-o $output/SV/$group.$sample.break.vcf	
		if [ -s $output/SV/$group.$sample.break.vcf	]
		then
			file=`echo $inputargs | sed -e '/-V/s///g'`
			rm $file
		fi
		rm $basedir/struct/break/$group/$sample/$sample.*.idx	
	done

	#Summaryzing Crest
	for sample in $tumor_list
	do
		rm $basedir/struct/$group.$sample.predSV.txt
		inputargs=""
		inputargs_filter=""
		for chr in $chrs
		do
			inputfile=$basedir/struct/crest/group/$sample.$chr.predSV.txt
			inputfile_filter=$basedir/struct/crest/group/$sample.$chr.filter.predSV.txt
			input_vcf=$basedir/struct/crest/group/$sample.$chr.raw.vcf
			input__filter_vcf=$basedir/struct/crest/group/$sample.$chr.filter.vcf
			if [ ! -f $inputfile ]
			then
				echo "ERROR:summaryze_struct. Crest: File $inputfile does not exist "
				exit 1
			fi
			cat $inputfile >> $basedir/struct/$group.raw.crest
			rm $inputfile
			cat $inputfile_filter >> $basedir/struct/$group.filter.crest
			rm $inputfile_filter
			inputargs="-V $input_vcf "$inputargs
			inputargs_filter="-V $input_vcf "$inputargs_filter
		done
		$java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
		-R $ref \
		-et NO_ET \
		-T CombineVariants \
		$inputargs \
		-o $output/SV/$group.$sample.crest.vcf
		
		$java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
		-R $ref \
		-et NO_ET \
		-T CombineVariants \
		$inputargs_filter \
		-o $output/SV/$group.$sample.crest.filter.vcf	
		
		if [ -s $output/SV/$group.$sample.crest.vcf ]
		then
			file=`echo $inputargs | sed -e '/-V/s///g'`
			rm $file
		fi
		
		if [ -s $output/SV/$group.$sample.crest.filter.vcf ]
		then
			file=`echo $inputargs_filter | sed -e '/-V/s///g'`
			rm $file
		fi	
		rm $basedir/struct/crest/group/$sample.*.idx
	done
	
	#Output a bedpe 
	cat $basedir/struct/$group.$sample.predSV.txt |\
	    awk 'BEGIN{OFS="\t"} {score=sprintf ("%.2f",(($4/$10)+($8/$11))/2); print $1,$2,$2+1,$5,$6,$6+1,"sample",score,$3,$7}' |\
	    sed -e "s/sample/$group.$sample/g" > $basedir/struct/$group.$sample.crest.bedpe
	echo `date`
fi	




