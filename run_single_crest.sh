#!/bin/sh

########################################################
###### 	SV CALLER FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			run_cnvnator.sh
######		Date:				09/26/2011
######		Summary:			Calls CNVnator
######		Input 
######		$1	=	samplename
######		$2	=	input_bam
######		$3	=	/path/to/output directory
######		$4	=	/path/to/run_info.txt
######		Output files:	BAM files. 
########################################################

if [ $# != 4 ]
then
	echo "\nUsage: samplename </path/to/realign directory/></path/to/output directory> </path/to/run_info.txt>";
else
	set -x
	echo `date`
	sample=$1
	input=$2
	output_dir=$3
	run_info=$4
	
	
########################################################	
######		Reading run_info.txt and assigning to variables
#SGE_TASK_ID=2
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
	queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	crest=$( cat $tool_info | grep -w '^CREST' | cut -d '=' -f2 )
	perllib=$( cat $tool_info | grep -w '^PERLLIB' | cut -d '=' -f2 )
	cap3=$( cat $tool_info | grep -w '^CAP3' | cut -d '=' -f2 )
	blat=$( cat $tool_info | grep -w '^BLAT' | cut -d '=' -f2 )
	blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
	blat_ref=$( cat $tool_info | grep -w '^BLAT_REF' | cut -d '=' -f2 )
	blat_server=$( cat $tool_info | grep -w '^BLAT_SERVER' | cut -d '=' -f2 )
	chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
	ref_genome=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
	output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)

	min_read=$( cat $tool_info | grep -w '^STRUCT_MIN_SUPPORT' | cut -d '=' -f2)
	min_id=$( cat $tool_info | grep -w '^STRUCT_MIN_IDENTITY' | cut -d '=' -f2)
	blacklist_sv=$( cat $tool_info | grep -w '^BLACKLIST_SV' | cut -d '=' -f2 )
	bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    PATH=$bedtools/:$PATH
#filter_sv_crest ()
#{
    # Filter SV's from file ($1), that have at leat $min_read supporting reads with
    # $min_id identity and that do not intersect $blacklist_sv and puts the result in
    # outfile ($2)
#    file=$1
#    outfile=$2
    
#    awk "((\$10>=$min_read)&&(\$11>=$min_read)&&(\$14>=$min_id)&&(\$16>=$min_id))" $file |\
#	awk '{ print $1"\t"$2"\t"$2+1"\t"$5"\t"$6"\t"$6+1}' |\ 
#	$bedtools/pairToBed -a stdin -b $blacklist_sv -type neither |\ 
#	$script_path/report_original.pl $file > $outfile
#}

########################################################	
######		
	input_bam=$input/chr${chr}.cleaned.bam
	
	PERL5LIB=$perllib
	PATH=$PATH:$blat:$crest
	mkdir -p $output_dir/$sample

	SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $input_bam -s $samtools`
	if [ $SORT_FLAG == 0 ]
	then
		echo "ERROR : run_crest_multi $file should be sorted"
		exit 1;
	fi


    if [ ! -s $input_bam.bai ]
	then
	    $samtools/samtools index $input_bam 
	fi

	mkdir -p $output_dir/$sample/log

	range=20000
	let blat_port+=$RANDOM%range
	status=`$blat/gfServer status localhost $blat_port | wc -l`;
    
    ln -s $input_bam $output_dir/$sample/$sample.chr${chr}.cleaned.bam
    ln -s $input_bam.bai $output_dir/$sample/$sample.chr${chr}.cleaned.bam.bai
    input_bam=$output_dir/$sample/$sample.chr${chr}.cleaned.bam
    
	if [ "$status" -le 1 ]
	then
		$blat/gfServer start localhost $blat_port -log=$output_dir/$sample/log/blat.$sample.$chr.txt $blat_ref  &
		sleep 2m
	fi

	$crest/extractSClip.pl -i $input_bam -r chr$chr --ref_genome $ref_genome -o $output_dir/$sample -p $sample

	status=`$blat/gfServer status localhost $blat_port | wc -l`;

	if [ "$status" -eq 0 ]
	then
        echo "WARNING : Crest_single: Blat server not available at port: $blat_port. Sample:$sample Chr:$chr "
        blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
        range=20000
        let blat_port+=$RANDOM%range
        status=`$blat/gfServer status localhost $blat_port | wc -l`;
        if [ "$status" -le 1 ]
        then
            rm $output_dir/$sample/log/blat.$sample.$chr.txt
			$blat/gfServer start localhost $blat_port -log=$output_dir/$sample/log/blat.$sample.$chr.txt $blat_ref  &
            sleep 2m
			echo "WARNING : Crest_single: Blat server available at port: $blat_port. Sample:$sample Chr:$chr "
        fi
    fi
	
	cd $crest    
	if [ -f $output_dir/$sample/$sample.chr$chr.cover ]
	then
		$crest/CREST.pl -f $output_dir/$sample/$sample.chr$chr.cover -d $input_bam \
		--ref_genome $ref_genome -t $blat_ref \
		--blatport $blat_port -blatserver localhost \
		--cap3 $cap3/cap3 \
		--min_sclip_len 12 \
		-o $output_dir/$sample -p $sample.$chr
		
        rm $input_bam
        rm $input_bam.bai
        perl $script_path/CREST2VCF.pl -i $output_dir/$sample/$sample.$chr.predSV.txt -f $ref_genome -o $output_dir/$sample/$sample.$chr.raw.vcf -s $sample -t $samtools
		if [ ! -s $output_dir/$sample/$sample.$chr.raw.vcf.fail ]
        then
            rm $output_dir/$sample/$sample.$chr.raw.vcf.fail
        fi  
		perl $script_path/vcfsort.pl ${ref_genome}.fai $output_dir/$sample/$sample.$chr.raw.vcf > $output_dir/$sample/$sample.$chr.raw.vcf.sort
		mv $output_dir/$sample/$sample.$chr.raw.vcf.sort $output_dir/$sample/$sample.$chr.raw.vcf
		
		awk "((\$10>=$min_read)&&(\$11>=$min_read)&&(\$14>=$min_id)&&(\$16>=$min_id))" $output_dir/$sample/$sample.$chr.predSV.txt | awk '{print $1"\t"$2"\t"$2+1"\t"$5"\t"$6"\t"$6+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $output_dir/$sample/$sample.$chr.predSV.txt > $output_dir/$sample/$sample.$chr.filter.predSV.txt
	
		### convert the output to VCF format
		perl $script_path/CREST2VCF.pl -i $output_dir/$sample/$sample.$chr.filter.predSV.txt -f $ref_genome -o $output_dir/$sample/$sample.$chr.filter.vcf -s $sample -t $samtools
		if [ ! -s $output_dir/$sample/$sample.$chr.filter.vcf.fail ]
        then
            rm $output_dir/$sample/$sample.$chr.filter.vcf.fail
        fi  
		perl $script_path/vcfsort.pl ${ref_genome}.fai $output_dir/$sample/$sample.$chr.filter.vcf > $output_dir/$sample/$sample.$chr.filter.vcf.sort
		mv $output_dir/$sample/$sample.$chr.filter.vcf.sort $output_dir/$sample/$sample.$chr.filter.vcf
		
		if [ -f $output_dir/$sample/$sample.$chr.predSV.txt ]
		then
			rm $output_dir/$sample/$sample.chr$chr.cover $output_dir/$sample/$sample.chr$chr.sclip.txt 
		else
			echo "ERROR: $output_dir/$sample/$sample.$chr.predSV.txt not created"
		fi			
	else
		echo "ERROR : Crest_single: extractSClip.pl, $output_dir/$sample/$sample.chr$chr.cover not created "
		touch $output_dir/$sample/$sample.$chr.predSV.txt
		perl $script_path/CREST2VCF.pl -i $output_dir/$sample/$sample.$chr.predSV.txt -f $ref_genome -o $output_dir/$sample/$sample.$chr.raw.vcf -s $sample -t $samtools
		rm $output_dir/$sample/$sample.$chr.raw.vcf.fail
		touch  $output_dir/$sample/$sample.$chr.filter.predSV.txt
		perl $script_path/CREST2VCF.pl -i $output_dir/$sample/$sample.$chr.filter.predSV.txt -f $ref_genome -o $output_dir/$sample/$sample.$chr.filter.vcf -s $sample -t $samtools
		rm $output_dir/$sample/$sample.$chr.filter.vcf.fail
	fi
    echo `date`
fi
