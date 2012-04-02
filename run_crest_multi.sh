#!/bin/sh

########################################################
###### 	SV CALLER FOR TUMOR/NORMAL PAIR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			run_crest_multi.sh
######		Date:				09/26/2011
######		Summary:			Calls Crest
######		Input 
######		$1	=	group name
######		$2	=	bam list (first is normal) : separated
######		$3	=	names of the samples : separated
######		$4	=	/path/to/output directory
######		$5	=	/path/to/run_info.txt
########################################################

if [ $# != 4 ]
then
	echo -e "Usage: Script to run crest on a paired sample \n <group name> </path/to/input directory> </path/to/output directory> </path/to/run_info.txt>";
else
	set -x
	echo `date`
	group=$1
	input=$2
	output_dir=$3
	run_info=$4
	
	

	########################################################	
	######		Reading run_info.txt and assigning to variables
#SGE_TASK_ID=1
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	crest=$( cat $tool_info | grep -w '^CREST' | cut -d '=' -f2 )
	picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
	perllib=$( cat $tool_info | grep -w '^PERLLIB' | cut -d '=' -f2 )
	blat=$( cat $tool_info | grep -w '^BLAT' | cut -d '=' -f2 )
	cap3=$( cat $tool_info | grep -w '^CAP3' | cut -d '=' -f2 )
	blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
	blat_ref=$( cat $tool_info | grep -w '^BLAT_REF' | cut -d '=' -f2 )
	
	blat_server=$( cat $tool_info | grep -w '^BLAT_SERVER' | cut -d '=' -f2 )
	chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
	ref_genome=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
	samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2)
	#bam=$input/chr$chr.cleaned.bam

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
#	cut -f 1,2,5,6 |\ 
#        awk '{ print $1"\t"$2"\t"($2+1)"\t"$3"\t"$4"\t"($4+1)}' |\ 
#	$bedtools/pairToBed -a stdin -b $blacklist_sv -type neither |\ 
#	$script_path/report_original.pl $file > $outfile
#}
	
	########################################################	
	######		

	echo `date`
	export PERL5LIB=$perllib
	PATH=$PATH:$blat:$crest:$perllib
	mkdir -p $output_dir/$group

	mkdir -p $output_dir/$group/log

	range=20000
	let blat_port+=$RANDOM%range

	status=`$blat/gfServer status localhost $blat_port | wc -l`;

	if [ "$status" -eq 0 ]
	then
		$blat/gfServer start localhost $blat_port -log=$output_dir/$group/log/blat.$group.$chr.txt $blat_ref  &
		sleep 2m
	fi

#i=1
#	for sample in $samples
#	do
#		ln -s $input/$group.$sample.chr$chr.bam $output_dir/$group/$sample.chr$chr.bam
		# check if BAM is sorted
#		file=$output_dir/$group/$sample.chr$chr.bam
#		SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $file -s $samtools`
#		if [ $SORT_FLAG == 0 ]
#		then
#			echo "ERROR : run_crest_multi $file should be sorted"
#			exit 1;
#		fi
		# check if BAM has an index
#		if [ ! -s $file.bai ]
#		then
#			$samtools/samtools index $file
#		fi

#		$crest/extractSClip.pl -i $file -r chr$chr --ref_genome $ref_genome -o $output_dir/$group -p $sample
#		let i=i+1
#	done

	let num_tumor=`echo $samples|tr " " "\n"|wc -l`-1
	normal_sample=`echo $samples| tr " " "\n" | head -n 1 `
	tumor_list=`echo $samples | tr " " "\n" | tail -$num_tumor`

	i=2
	for file in $tumor_list
	do
		sample=`echo $samples| tr " " "\n"| head -n $i | tail -1`
		if [ ! -f $output_dir/$group/$sample.chr$chr.cover ]
		then
			echo "ERROR : run_crest_multi, file $output_dir/$group/$sample.chr$chr.cover was not created, error in extracSClip.pl"
		fi

		if [ ! -f $output_dir/$group/$sample.chr$chr.sclip.txt ]
		then
			echo "ERROR: run_crest_multi, file $output_dir/$group/$sample.chr$chr.sclip.txt does not exist, error in extracSClip.pl "
		fi

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
			fi
		fi

		$crest/CREST.pl -f $output_dir/$group/$sample.chr$chr.cover \
        -d $output_dir/$group/${file}.chr$chr.bam -g $output_dir/$group/${normal_sample}.chr$chr.bam \
		--ref_genome $ref_genome -t $blat_ref \
		--blatport $blat_port -blatserver localhost \
		--min_sclip_len 12 \
		--cap3 $cap3/cap3 --min_sclip_reads 14 \
		-r chr$chr -o $output_dir/$group -p $sample.$chr

        if [ -f $output_dir/$group/${file}.$chr.predSV.txt ] 
        then
            rm $output_dir/$group/${file}.chr$chr.bam $output_dir/$group/${file}.chr$chr.bam.bai
        fi
        perl $script_path/CREST2VCF.pl -i $output_dir/$group/${file}.$chr.predSV.txt -f $ref_genome -o $output_dir/$group/$file.$chr.raw.vcf -s $file -t $samtools
        perl $script_path/vcfsort.pl ${ref_genome}.fai $output_dir/$group/$file.$chr.raw.vcf > $output_dir/$group/$file.$chr.raw.vcf.sort
        mv $output_dir/$group/$file.$chr.raw.vcf.sort $output_dir/$group/$file.$chr.raw.vcf
		if [ ! -s $output_dir/$group/$file.$chr.raw.vcf.fail ]
        then
            rm $output_dir/$group/$file.$chr.raw.vcf.fail
        fi  

		#awk "((\$10>=$min_read)&&(\$11>=$min_read)&&(\$14>=$min_id)&&(\$16>=$min_id))" $output_dir/$group/$sample.$chr.predSV.txt >> $output_dir/$group/$sample.$chr.filter.predSV.txt
		
		if [ -s $output_dir/$group/${file}.$chr.predSV.txt ]
        then
            awk "((\$10>=$min_read)&&(\$11>=$min_read)&&(\$14>=$min_id)&&(\$16>=$min_id))" $output_dir/$group/${file}.$chr.predSV.txt | awk '{print $1"\t"$2"\t"$2+1"\t"$5"\t"$6"\t"$6+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $output_dir/$sample/$sample.$chr.predSV.txt > $output_dir/$group/${file}.$chr.filter.predSV.txt
		else
            touch $output_dir/$group/${file}.$chr.filter.predSV.txt
        fi

		#filter_sv_crest $output_dir/$sample/$sample.$chr.predSV.txt $output_dir/$sample/$sample.$chr.filter.predSV.txt

		perl $script_path/CREST2VCF.pl -i $output_dir/$group/${file}.$chr.filter.predSV.txt -f $ref_genome -o $output_dir/$group/${file}.$chr.filter.vcf -s $file -t $samtools
		if [ ! -s $output_dir/$group/${file}.$chr.filter.vcf.fail ]
        then
            rm $output_dir/$group/${file}.$chr.filter.vcf.fail
        fi  
		perl $script_path/vcfsort.pl ${ref_genome}.fai $output_dir/$group/${file}.$chr.filter.vcf > $output_dir/$group/${file}.$chr.filter.vcf.sort
		mv $output_dir/$group/${file}.$chr.filter.vcf.sort $output_dir/$group/${file}.$chr.filter.vcf
		
		### vcf converter for CREST output to VCF
		if [ ! -f $output_dir/$group/${file}.$chr.filter.predSV.txt ]
		then
			echo "ERROR: run_crest_multi, file $output_dir/$group/${file}.$chr.filter.predSV.txt was not created, error in CREST.pl" 
			exit 1
		else
			rm $output_dir/$group/$file.chr$chr.sclip.txt $output_dir/$group/$file.chr$chr.cover
		fi
		let i=i+1
	done
    rm $output_dir/$group/${normal_sample}.chr$chr.bam $output_dir/$group/${normal_sample}.chr$chr.bam.bai
    rm $output_dir/$group/${normal_sample}.chr$chr.sclip.txt $output_dir/$group/${normal_sample}.chr$chr.cover
fi
