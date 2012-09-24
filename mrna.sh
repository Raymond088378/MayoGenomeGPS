#!/bin/bash

########################################################
###### 	MASTER SCRIPT FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			whole_genome_pipeline.sh
######		Date:				06/05/2012
######		Summary:			Master script encompassing subscripts for alignment, remove duplicates, realignment, 
######                          	recalibration, fastqc, variant calling and final filtering of variants.
######		Input files:		$1	=	/path/to/run_info.txt
######		Output files:		variant VCF files. Look at subscripts for granular description of output files.
######		TWIKI:				http://bioinformatics.mayo.edu/BMI/bin/view/Main/BioinformaticsCore/Analytics/WholeGenomeWo
########################################################

if [ $# != 1 ]
then	
	echo "Usage: <Please specify path to run_info.txt file> ";
else
	set -x
	echo `date`
	run_info=$1
	dos2unix $run_info
	if [ `perl -v 1 | tr " " "\n" | grep '^v'` != "v5.10.0" ]
	then
		echo -e "\nperl path is not correct in your enviornment"
		echo "Perl path should point to /usr/local/biotools/perl/5.10.0/bin/perl if the user use the command which perl, user can change this using mayobiotools"
		exit 1;
	fi
	dir_info=`dirname $run_info`
	if [ "$dir_info" = "." ]
	then
		echo "ERROR : run_info=$run_info should be specified as a complete path\n";
		exit 1;
	fi
	if [ ! -s $run_info ]
	then
		echo "ERROR : run_info=$run_info does not exist \n";
		exit 1;
	fi
	## removing trailing and leading spaces from run ifno file
	cat $run_info|sed 's/[ \t]*$//' |sed 's/^[ \t]*//;s/[ \t]*$//'> $run_info.tmp
	mv $run_info.tmp $run_info
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	dos2unix $sample_info
	dos2unix $tool_info
	## removing trailing and leading spaces
	cat $sample_info|sed 's/[ \t]*$//' |sed 's/^[ \t]*//;s/[ \t]*$//'> $sample_info.tmp
	mv $sample_info.tmp $sample_info
	cat $tool_info |sed 's/[ \t]*$//' |sed 's/^[ \t]*//;s/[ \t]*$//' > $tool_info.tmp
	mv $tool_info.tmp $tool_info
	
	input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
	output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2)
	groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	type=$( cat $run_info | grep -w '^TOOL' | cut -d '=' -f2|tr "[a-z]" "[A-Z]")
	version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
	queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	all_sites=$( cat $tool_info | grep -w '^EMIT_ALL_SITES' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	aligner=$( cat $run_info | grep -w '^ALIGNER' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]")
	upload_tb=$( cat $run_info | grep -w '^UPLOAD_TABLEBROWSER' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	numchrs=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
	paired=$( cat $run_info | grep -w '^PAIRED' | cut -d '=' -f2)
	threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2)
	variant_type=$(cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")   
	somatic_caller=$(cat $run_info | grep -w '^SOMATIC_CALLER' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]") 
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
	if [[ $somatic_caller == "JOINTSNVMIX" || $somatic_caller == "BEAUTY_EXOME" ]]
	then
		python_path=`which python`
		if [ $python_path != "/usr/local/biotools/python/2.7/bin/python" ]
		then
			echo -e "\n python path is not correct in your enviorment"
			echo " Python path should point to /usr/local/biotools/python/2.7/bin/python if the user use the command which python, user can change this using mayobiotools"
			exit 1;
		fi    
	fi
	#################################################
	### validate the config file
	perl $script_path/check_config.pl $run_info > $run_info.configuration_errors.log
	if [ `cat $run_info.configuration_errors.log | wc -l` -gt 0 ]
	then
		echo "Configuration files are malformed: look at the erros in $run_info.configuration_errors.log "
		exit 1;
	else
		rm $run_info.configuration_errors.log
	fi	
	### create folders
	$script_path/create_folder.sh $run_info
	output_dir=$output/$PI/$tool/$run_num
	if [ -f $output_dir/folder_exist.log ]
	then
		echo "ERROR: folder already exist"
		exit 1;
	fi	
	## copy cofig files
	$script_path/copy_config.sh $output_dir $run_info
	job_ids_dir=$output_dir/job_ids
	output_align=$output_dir/alignment
	if [ $analysis != "alignment" ]
	then
		output_OnTarget=$output_dir/OnTarget
		output_annot=$output_dir/annotation
		TempReports=$output_dir/TempReports
		sift=$output_annot/SIFT
		snpeff=$output_annot/SNPEFF
		polyphen=$output_annot/POLYPHEN
	fi
	##########################################################
	echo -e "${tool} analysis for ${run_num} for ${PI} " >> $output_dir/log.txt
	START=`date`
	echo -e "Analysis started at:" >> $output_dir/log.txt
	echo -e "${START}" >>  $output_dir/log.txt
	echo -e "TOOL INFO file used : $tool_info" >>  $output_dir/log.txt
	echo -e "SAMPLE INFO file used : $sample_info" >>  $output_dir/log.txt
	echo -e "RUN INFO  file used : $run_info" >>  $output_dir/log.txt
	
	###########################################################
	#### sge paramtersff
	TO=`id |awk -F '(' '{print $2}' | cut -f1 -d ')'`
	args="-V -wd $output_dir/logs -q $queue -m a -M $TO -l h_stack=10M"
	#############################################################
	### get the identification number for this run.
	
	if [ $multi_sample != "YES" ]
	then
		echo "Single sample"
		numsamples=$(cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
		for sample in `echo $samples | tr ":" "\n"`
		do            
			sleep 5
			if [ $analysis != "annotation" ]
			then
				align_dir=$output_dir/alignment/$sample
				bamfile=$sample.sorted.bam
			fi
			if [ $analysis == "mayo" -o $analysis == "external" -o $analysis == "alignment" ]
			then
				mkdir -p $align_dir
				if [ $paired == 1 ]
				then
					let numfiles=(`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" |wc -l`)/2
				else
					let numfiles=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" |wc -l`
				fi	

				if [ $analysis == "mayo" ]
				then
					for ((i=1; i <=$numfiles; i++));
					do
						$script_path/dashboard.sh $sample $run_info Beginning started $i
					done
				fi	
				if [ $aligner == "novoalign" ]
				then
					echo "novoalign is used as aligner"
					sleep 5
					qsub $args -N $type.$version.align_novo.$sample.$run_num -l h_vmem=4G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_novo.sh $sample $output_dir $run_info
				elif [ $aligner == "bwa" ]
				then
					echo "bwa is used as aligner"
					sleep 5
					qsub $args -N $type.$version.align_read_bwa.R1.$sample.$run_num -l h_vmem=1G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_read_bwa.sh $sample $output_dir 1 $run_info
					if [ $paired == 1 ]
					then
						sleep 5
						qsub $args -N $type.$version.align_read_bwa.R2.$sample.$run_num -l h_vmem=1G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_read_bwa.sh $sample $output_dir 2 $run_info
						hold="-hold_jid $type.$version.align_read_bwa.R1.$sample.$run_num,$type.$version.align_read_bwa.R2.$sample.$run_num"
					else
						hold="-hold_jid $type.$version.align_read_bwa.R1.$sample.$run_num"
					fi	
					sleep 5
					qsub $args -N $type.$version.align_bwa.$sample.$run_num -l h_vmem=3G -pe threaded $threads $hold -t 1-$numfiles:1 $script_path/align_bwa.sh $sample $output_dir $run_info
				else
					echo "Doesn't support the aligner"
				fi	
				if [ $aligner == "bwa" ]
				then
					hold="-hold_jid $type.$version.align_bwa.$sample.$run_num"
				elif [ $aligner == "novoalign" ]
				then
					hold="-hold_jid $type.$version.align_novo.$sample.$run_num"
				fi    
				sleep 5
				qsub $args -N $type.$version.processBAM.$sample.$run_num -pe threaded $threads -l h_vmem=4G $hold $script_path/processBAM.sh $align_dir $sample $run_info 	
				if [ $analysis != "alignment" ]
				then
					sleep 5
					qsub $args -N $type.$version.extract_reads_bam.$sample.$run_num -l h_vmem=8G -hold_jid $type.$version.processBAM.$sample.$run_num $script_path/extract_reads_bam.sh $align_dir $bamfile $run_info $output_dir/IGV_BAM
				fi
			elif [ $analysis == "realignment" -o $analysis == "realign-mayo" ]
			then
				mkdir -p $align_dir
				infile=`cat $sample_info | grep -w ^BAM:${sample} | cut -d '=' -f2`
				num_bams=`echo $infile | tr " " "\n" | wc -l`
				for ((i=1; i <=$num_bams; i++));
				do
					bam=`echo $infile | awk -v num=$i '{print $num}'`
					$samtools/samtools view -H $input/$bam 2> $align_dir/$sample.$i.sorted.bam.log
					if [ `cat $align_dir/$sample.$i.sorted.bam.log | wc -l` -gt 0 ]
					then
						echo "$input/$bam : bam file is truncated or corrupted" 	
						exit 1;
					else
						rm $align_dir/$sample.$i.sorted.bam.log
					fi	
					ln -s $input/$bam $align_dir/$sample.$i.sorted.bam
					$script_path/dashboard.sh $sample $run_info Beginning started $i
				done  
				sleep 5
				qsub $args -N $type.$version.processBAM.$sample.$run_num -pe threaded $threads -l h_vmem=4G $script_path/processBAM.sh $align_dir $sample $run_info
				sleep 5
				qsub $args -N $type.$version.extract_reads_bam.$sample.$run_num -l h_vmem=8G -hold_jid $type.$version.processBAM.$sample.$run_num $script_path/extract_reads_bam.sh $align_dir $bamfile $run_info $output_dir/IGV_BAM
			fi    
			if [[ $analysis == "mayo" || $analysis == "external" || $analysis == "realignment" || $analysis == "variant" || $analysis == "realign-mayo" ]]
			then
				realign_dir=$output_dir/realign/$sample
				variant_dir=$output_dir/variants/$sample
				mkdir -p $realign_dir $variant_dir
				if [ $analysis == "variant" ]
				then
					infile=`cat $sample_info | grep -w ^BAM:${sample} | cut -d '=' -f2`
					num_bams=`echo $infile | tr " " "\n" | wc -l`
					for ((i=1; i <=$num_bams; i++));
					do
						bam=`echo $infile | awk -v num=$i '{print $num}'`
						$samtools/samtools view -H $input/$bam 2> $align_dir/$sample.$i.sorted.bam.fix.log
						if [ `cat $align_dir/$sample.$i.sorted.bam.fix.log | wc -l` -gt 0 ]
						then
							echo "$input/$bam : bam file is truncated or corrupted" 	
							exit 1;
						else
							rm $align_dir/$sample.$i.sorted.bam.fix.log
						fi
						ln -s $input/$bam $realign_dir/$sample.$i.sorted.bam						
					done
					sleep 5
					qsub $args -N $type.$version.reformat_BAM.$sample.$run_num -l h_vmem=8G $script_path/reformat_BAM.sh $realign_dir $sample $run_info	
					sleep 5
					qsub $args -N $type.$version.extract_reads_bam.$sample.$run_num -l h_vmem=8G -hold_jid $type.$version.reformat_BAM.$sample.$run_num $script_path/extract_reads_bam.sh $realign_dir $bamfile $run_info $output_dir/IGV_BAM
					sleep 5
					qsub $args -N $type.$version.split_bam_chr.$sample.$run_num -hold_jid $type.$version.reformat_BAM.$sample.$run_num -l h_vmem=2G -t 1-$numchrs:1 $script_path/split_bam_chr.sh $realign_dir $sample $run_info
					variant_id="$type.$version.split_bam_chr.$sample.$run_num"
				else
					sleep 5
					qsub $args -N $type.$version.realign_recal.$sample.$run_num -hold_jid $type.$version.processBAM.$sample.$run_num -l h_vmem=8G -t 1-$numchrs:1 $script_path/realign_recal.sh $align_dir $bamfile $sample $realign_dir $run_info 1	
					variant_id="$type.$version.realign_recal.$sample.$run_num"
				fi
				sleep 5
				qsub $args -N $type.$version.igv_bam.$sample.$run_num -l h_vmem=2G -hold_jid $variant_id,$type.$version.extract_reads_bam.$sample.$run_num $script_path/igv_bam.sh $output_dir/realign $output_dir/IGV_BAM $sample $output_dir/alignment $run_info
				sleep 5
				qsub $args -N $type.$version.variants.$sample.$run_num -hold_jid $variant_id -pe threaded $threads -l h_vmem=6G -t 1-$numchrs:1 $script_path/variants.sh $realign_dir $sample $variant_dir 1 $run_info
				sleep 5
				qsub $args -N $type.$version.concat_mrna_variants.$sample.$run_num -hold_jid $type.$version.variants.$sample.$run_num -l h_vmem=4G $script_path/concat_mrna_variants.sh $variant_dir $sample $output_dir/variants/ $run_info
			fi        
		done	
	fi
	echo `date`
fi
