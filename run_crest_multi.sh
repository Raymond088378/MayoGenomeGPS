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
	SGE_TASK_ID=2
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
	queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
	lqueue=$( cat $run_info | grep -w '^LQUEUE' | cut -d '=' -f2)
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
	min_read=$( cat $tool_info | grep -w '^STRUCT_MIN_SUPPORT' | cut -d '=' -f2)
	min_id=$( cat $tool_info | grep -w '^STRUCT_MIN_INDENTITY' | cut -d '=' -f2)
	samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2)
	bam=$input/chr$chr.cleaned.bam
	
	########################################################	
	######		

	echo `date`
	PERL5LIB=$perllib
	PATH=$PATH:$blat:$crest
	mkdir $output_dir/$group

	mkdir -p $output_dir/$group/log

	range=20000
	let blat_port+=$RANDOM%range

	status=`$blat/gfServer status localhost $blat_port | wc -l`;

	if [ "$status" -eq 0 ]
	then
		$blat/gfServer start localhost $blat_port -log=$output_dir/$group/log/blat.$group.$chr.txt $blat_ref  &
		sleep 2m
	fi

	i=1
	for sample in $samples
	do
		sam=`echo $samples | tr " " "\n" | grep -v "$sample" | tr "\n" " " `
		gr=""
		for s in $sam
		do
			a="ID:$s|";
			gr="$gr $a"
		done
		gr=`echo $gr |  sed "s/|$//"`
		$samtools/samtools view -b -r $sample $bam > $output_dir/$group/$sample.chr$chr.bam
        $samtools/samtools view -H $output_dir/$group/$sample.chr$chr.bam | grep -E -v $gr | $samtools/samtools reheader - $output_dir/$group/$sample.chr$chr.bam > $output_dir/$group/$sample.chr$chr.re.bam
        mv $output_dir/$group/$sample.chr$chr.re.bam $output_dir/$group/$sample.chr$chr.bam
		# check if BAM is sorted
		file=$output_dir/$group/$sample.chr$chr.bam
		SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $file -s $samtools`
		if [ $SORT_FLAG == 0 ]
		then
			echo "ERROR : run_crest_multi $file should be sorted"
			exit 1;
		fi

		# check if BAM has an index
		if [ ! -s $file.bai ]
		then
			$samtools/samtools index $file $file.bai
		fi

		$crest/extractSClip.pl -i $file -r chr$chr --ref_genome $ref_genome -o $output_dir/$group -p $sample
		let i=i+1
	done

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
			exit 1
		fi

		if [ ! -f $output_dir/$group/$sample.chr$chr.sclip.txt ]
		then
			echo "ERROR: run_crest_multi, file $output_dir/$group/$sample.chr$chr.sclip.txt does not exist, error in extracSClip.pl "
			exit 1
		fi

		status=`$blat/gfServer status localhost $blat_port | wc -l`;

		if [ "$status" -eq 0 ]
		then
			echo "ERROR: run_crest_multi: Blat server not available at port: $blat_port. Group:$group Chr:$chr "
			exit 1
		fi

		$crest/CREST.pl -f $output_dir/$group/$sample.chr$chr.cover \
		-d $file -g ${normal_sample}.chr$chr.bam \
		--ref_genome $ref_genome -t $blat_ref \
		--blatport $blat_port -blatserver localhost \
		--min_sclip_len 12 \
		--cap3 $cap3/cap3 --min_sclip_reads 14 \
		-r chr$chr -o $output_dir/$group -p $sample.$chr

		#Filter output files with at least $min_reads reads of support
        perl $script_path/CREST2VCF.pl -i $output_dir/$group/$sample.$chr.predSV.txt -f $ref_genome -o $output_dir/$group/$sample.$chr.raw.vcf -s $sample -t $samtools
        perl $script_path/vcfsort.pl ${ref_genome}.fai $output_dir/$group/$sample.$chr.raw.vcf > $output_dir/$group/$sample.$chr.raw.vcf.sort
        mv $output_dir/$group/$sample.$chr.raw.vcf.sort $output_dir/$group/$sample.$chr.raw.vcf
		if [ ! -s $output_dir/$group/$sample.$chr.raw.vcf.fail ]
        then
            rm $output_dir/$group/$sample.$chr.raw.vcf.fail
        fi  

		awk "((\$10>=$min_read)&&(\$11>=$min_read)&&(\$14>=$min_id)&&(\$16>=$min_id))" $output_dir/$group/$sample.$chr.predSV.txt >> $output_dir/$group/$sample.$chr.filter.predSV.txt
		perl $script_path/CREST2VCF.pl -i $output_dir/$group/$sample.$chr.filter.predSV.txt -f $ref_genome -o $output_dir/$group/$sample.$chr.filter.vcf -s $sample -t $samtools
		if [ ! -s $output_dir/$group/$sample.$chr.filter.vcf.fail ]
        then
            rm $output_dir/$group/$sample.$chr.filter.vcf.fail
        fi  
		perl $script_path/vcfsort.pl ${ref_genome}.fai $output_dir/$group/$sample.$chr.filter.vcf > $output_dir/$group/$sample.$chr.filter.vcf.sort
		mv $output_dir/$group/$sample.$chr.filter.vcf.sort $output_dir/$group/$sample.$chr.filter.vcf
		
		### vcf converter for CREST output to VCF
		if [ ! -f $output_dir/$group/$sample.$chr.predSV.txt ]
		then
			echo "ERROR: run_crest_multi, file $output_dir/$group/$sample.$chr.predSV.txt was not created, error in CREST.pl" 
			exit 1
		else
			rm $output_dir/$group/$sample.chr$chr.sclip.txt $output_dir/$group/$sample.chr$chr.cover
		fi
		let i=i+1
	done
fi
