#!/bin/bash

########################################################
###### 	MASTER SCRIPT FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			whole_genome_pipeline.sh
######		Date:				09/17/2012
######		Summary:			Master script encompassing subscripts for alignment, remove duplicates, realignment, 
######                          	recalibration, fastqc, variant calling and final filtering of variants.
######		Input files:		$1	=	/path/to/run_info.txt
######		Output files:		variant VCF files. Look at subscripts for granular description of output files.
######		TWIKI:				http://bioinformatics.mayo.edu/BMI/bin/view/Main/BioinformaticsCore/Analytics/WholeGenomeWo
########################################################

if [ $# != 1 ]
then	
	echo -e "Wrapper script to submit all the jobs for dna-seq workflow\nUsage: ./whole_genome.sh <Please specify full /path/to/run_info file>";
else
	set -x
	echo `date`
	run_info=$1
	dos2unix $run_info
	if [ `perl -v 2>&1 | tr " " "\n" | grep '^v'` != "v5.10.0" ]
	then
		echo -e "\nperl path is not correct in your enviornment"
		echo "Perl path should point to /usr/local/biotools/perl/5.10.0/bin/perl if the user use the command which perl and user can change this using mayobiotools"
		exit 1;
	fi
	dir_info=`dirname $run_info`
	if [ "$dir_info" = "." ]
	then
		echo -e "ERROR : run_info=$run_info should be specified as a complete path\n";
		exit 1;
	fi
	if [ ! -s $run_info ]
	then
		echo -e "ERROR : run_info=$run_info does not exist\n";
		exit 1;
	fi
	
	## removing trailing and leading spaces from run ifno file
	cat $run_info | sed 's/^[ \t]*//;s/[ \t]*$//' > $run_info.tmp
	mv $run_info.tmp $run_info
	
	#### check for unique identification number
	
	identify=$( cat $run_info | grep -w '^IDENTIFICATION_NUMBER' | cut -d '=' -f2)
	if echo $identify | egrep -q '^[0-9]+$'
	then
		echo -e "HURRAY !!! you are good to run the workflow"
	else
		echo -e "ERROR : unique identification for the workflow was not generated, please run the unique_id.sh script to generate the same before running the workflow script."
		exit 1;
	fi
	
	temp_id=$( cat $run_info | grep -w '^TEMPORARY_ID' | cut -d '=' -f2)
	if echo $temp_id | egrep -q '^[0-9]+$'
	then
		echo -e "ERROR : unique identification for the workflow is not new, please run the unique_id.sh script to generate the same before running the workflow script."
		exit 1;	
	else
		echo -e "HURRAY !!! you are good to run the workflow"
		echo -e "TEMPORARY_ID=$identify" >> $run_info
	fi
	
	### check for tool info file
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	if [ ! -s $tool_info ]
	then
		echo "ERROR : tool_info=$tool_info does not exist \n";
		exit 1;
	else
		dos2unix $tool_info
		cat $tool_info | sed 's/^[ \t]*//;s/[ \t]*$//' > $tool_info.tmp
		mv $tool_info.tmp $tool_info
	fi
	
	### check for sample info file
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	if [ ! -s $sample_info ]
	then
		echo "ERROR : sample_info=$sample_info does not exist \n";
		exit 1;
	else
		dos2unix $sample_info
		cat $sample_info | sed 's/^[ \t]*//;s/[ \t]*$//' > $sample_info.tmp
		mv $sample_info.tmp $sample_info
	fi
	
	### check for sample info file
	memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	if [ ! -s $memory_info ]
	then
		echo "ERROR : memory_info=$memory_info does not exist \n";
		exit 1;
	else
		dos2unix $memory_info
		cat $memory_info | sed 's/^[ \t]*//;s/[ \t]*$//' > $memory_info.tmp
		mv $memory_info.tmp $memory_info
	fi
	
	#### extract paths
	input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
	output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2)
	groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	type=$( cat $run_info | grep -w '^TOOL' | cut -d '=' -f2|tr "[a-z]" "[A-Z]")
	version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
	queue=$( cat $tool_info | grep -w '^QUEUE' | sed -e '/QUEUE=/s///g')
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	all_sites=$( cat $tool_info | grep -w '^EMIT_ALL_SITES' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	aligner=$( cat $run_info | grep -w '^ALIGNER' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]")
	upload_tb=$( cat $tool_info | grep -w '^UPLOAD_TABLEBROWSER' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	numchrs=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
	paired=$( cat $run_info | grep -w '^PAIRED' | cut -d '=' -f2)
	threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2)
	variant_type=$(cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")   
	somatic_caller=$(cat $run_info | grep -w '^SOMATIC_CALLER' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]") 
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
	limit=$( cat $tool_info | grep -w '^JOB_LIMIT' | cut -d '=' -f2 )
	info=$(cat $run_info | grep -w '^SAMPLEINFORMATION' | cut -d '=' -f2 )
	workflow=$( cat $run_info | grep '^TOOL=' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
	if [[ $somatic_caller == "JOINTSNVMIX" || $somatic_caller == "BEAUTY_EXOME" ]]
	then
		python_path=`which python`
		if [ `python -V  2>&1 | tr " " "\n" | grep -v Python` !=  "2.7" ]
		then
			echo -e "\n python path is not correct in your enviorment"
			echo " Python path should point to /usr/local/biotools/python/2.7/bin/python if the user use the command which python, user can change this using mayobiotools"
			exit 1;
		fi    
	fi
	
	#################################################
	### validate the config file
	$script_path/check_config.pl $run_info > $run_info.configuration_errors.log
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
	config=$output_dir/config
        
	if [ -f $output_dir/folder_exist.log ]
	then
		echo "ERROR: folder already exist"
		exit 1;
	fi	
	
	## copy cofig files
	$script_path/copy_config.sh $output_dir $run_info
	run_info=$output_dir/run_info.txt
	
	### modify the run info file to use configurations in the output folder
	add=`date +%D`
	cat $run_info | grep -w -v -E '^TOOL_INFO|^SAMPLE_INFO|^MEMORY_INFO' > $run_info.tmp
	echo -e "TOOL_INFO=$config/tool_info.txt\nSAMPLE_INFO=$config/sample_info.txt\nMEMORY_INFO=$config/memory_info.txt\nDATE=$add" | cat $run_info.tmp - > $run_info
	rm $run_info.tmp
	mv $run_info $config/run_info.txt
	run_info=$config/run_info.txt
	## tool info file
	tool_info=$output_dir/tool_info.txt
	mv $tool_info $config/tool_info.txt
	tool_info=$config/tool_info.txt
	## sample info file
	sample_info=$output_dir/sample_info.txt
	mv $sample_info $config/sample_info.txt
	sample_info=$config/sample_info.txt
	### memory info 
	memory_info=$output_dir/memory_info.txt
	mv $memory_info $config/memory_info.txt
	memory_info=$config/memory_info.txt
	output_align=$output_dir/alignment
	if [ $analysis != "alignment" ]
	then
		output_realign=$output_dir/realign/
		output_variant=$output_dir/variants
		output_OnTarget=$output_dir/OnTarget
		output_annot=$output_dir/annotation
		TempReports=$output_dir/TempReports
		sift=$output_annot/SIFT
		snpeff=$output_annot/SNPEFF
		polyphen=$output_annot/POLYPHEN
		igv=$output_dir/IGV_BAM
		RSample=$output_dir/Reports_per_Sample/
		annot=$output_dir/Reports_per_Sample/ANNOT
		sv=$output_dir/Reports_per_Sample/SV
		numbers=$output_dir/numbers
		struct=$output_dir/struct/
		cnv=$output_dir/cnv/
		circos=$output_dir/circos/
	fi
	
	##########################################################
	echo -e "${tool} analysis for ${run_num} for ${PI}\n" >> $output_dir/log.txt
	START=`date`
	echo -e "Analysis started at:" >> $output_dir/log.txt
	echo -e "${START}\n" >>  $output_dir/log.txt
	echo -e "Configuration files used in the analysis:\n" >> $output_dir/log.txt
	echo -e "TOOL INFO file used : $tool_info" >>  $output_dir/log.txt
	echo -e "SAMPLE INFO file used : $sample_info" >>  $output_dir/log.txt
	echo -e "RUN INFO  file used : $run_info" >>  $output_dir/log.txt
	echo -e "MEMORY INFO  file used : $memory_info" >>  $output_dir/log.txt
	###########################################################
	#### sge paramters
	TO=$USER
	args="-V -wd $output_dir/logs -q $queue -m a -M $TO -l h_stack=10M"
	echo -e "\nRCF arguments used : $args\n" >> $output_dir/log.txt
	echo -e "Started the ${tool} analysis for ${run_num} for ${PI}\n\n${info}\n\nCourtesy: $workflow $version" | mailx -v -s "Analysis Started" -c Kahl.Jane@mayo.edu "$TO"
        #############################################################
	
	if [ $multi_sample != "YES" ]
	then
		echo "Single sample"
		numsamples=$(cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
		for sample in `echo $samples | tr ":" "\n"`
		do            
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
				if [ $numfiles -eq 0 ]
				then
					echo "sample info file is not properly configured"
					exit 1;
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
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^align_novo' | cut -d '=' -f2)
					qsub_args="-N $type.$version.align_novo.$sample.$run_num.$identify -pe threaded $threads -t 1-$numfiles:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/align_novo.sh $sample $output_dir $run_info
				elif [ $aligner == "bwa" ]
				then
					echo "bwa is used as aligner"
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^align_read_bwa' | cut -d '=' -f2)
					qsub_args="-N $type.$version.align_read_bwa.R1.$sample.$run_num.$identify -pe threaded $threads -t 1-$numfiles:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/align_read_bwa.sh $sample $output_dir 1 $run_info
					if [ $paired == 1 ]
					then
						$script_path/check_qstat.sh $limit
						mem=$( cat $memory_info | grep -w '^align_read_bwa' | cut -d '=' -f2)
						qsub_args="-N $type.$version.align_read_bwa.R2.$sample.$run_num.$identify -pe threaded $threads -t 1-$numfiles:1 -l h_vmem=$mem"
						qsub $args $qsub_args $script_path/align_read_bwa.sh $sample $output_dir 2 $run_info
						hold="-hold_jid $type.$version.align_read_bwa.R1.$sample.$run_num.$identify,$type.$version.align_read_bwa.R2.$sample.$run_num.$identify"
					else
						hold="-hold_jid $type.$version.align_read_bwa.R1.$sample.$run_num.$identify"
					fi	
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^align_bwa' | cut -d '=' -f2)
					qsub_args="-N $type.$version.align_bwa.$sample.$run_num.$identify -t 1-$numfiles:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/align_bwa.sh $sample $output_dir $run_info
				else
					echo "Doesn't support the aligner"
				fi	
				if [ $aligner == "bwa" ]
				then
					hold="-hold_jid $type.$version.align_bwa.$sample.$run_num.$identify"
				elif [ $aligner == "novoalign" ]
				then
					hold="-hold_jid $type.$version.align_novo.$sample.$run_num.$identify"
				fi    
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^processBAM' | cut -d '=' -f2)
				qsub_args="-N $type.$version.processBAM.$sample.$run_num.$identify $hold -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/processBAM.sh $align_dir $sample $run_info 	
				if [ $analysis != "alignment" ]
				then
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^extract_reads_bam' | cut -d '=' -f2)
					qsub_args="-N $type.$version.extract_reads_bam.$sample.$run_num.$identify -hold_jid $type.$version.processBAM.$sample.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/extract_reads_bam.sh $align_dir $bamfile $run_info $igv
				fi
			elif [[ $analysis == "realignment" || $analysis == "realign-mayo" ]]
			then
				mkdir -p $align_dir
				infile=`cat $sample_info | grep -w ^BAM:${sample} | cut -d '=' -f2`
				num_bams=`echo $infile | tr " " "\n" | wc -l`
				if [ $num_bams -eq 0 ]
				then
					echo "sample info file is not properly configured"
					exit 1;
				fi
				if [ $analysis == "realign-mayo" ]
				then
					$script_path/dashboard.sh $sample $run_info Beginning started 
				fi	
				for ((i=1; i <=$num_bams; i++));
				do	
					bam=`echo $infile | awk -v num=$i '{print $num}'`
					$samtools/samtools view -H $input/$bam 1>$align_dir/$sample.$i.sorted.header 2> $align_dir/$sample.$i.sorted.bam.log
					if [[ `cat $align_dir/$sample.$i.sorted.bam.log | wc -l` -gt 0 || `cat $align_dir/$sample.$i.sorted.header | wc -l` -le 0 ]]
					then
						$script_path/errorlog.sh $input/$bam whole_genome.sh ERROR "truncated or corrupted"
						exit 1;
					else
						rm $align_dir/$sample.$i.sorted.bam.log
					fi	
					rm $align_dir/$sample.$i.sorted.header
					ln -s $input/$bam $align_dir/$sample.$i.sorted.bam
					$script_path/dashboard.sh $sample $run_info Beginning started $i
				done  
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^processBAM' | cut -d '=' -f2)
				qsub_args="-N $type.$version.processBAM.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/processBAM.sh $align_dir $sample $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^extract_reads_bam' | cut -d '=' -f2)
				qsub_args="-N $type.$version.extract_reads_bam.$sample.$run_num.$identify -hold_jid $type.$version.processBAM.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/extract_reads_bam.sh $align_dir $bamfile $run_info $igv
			fi    
			if [[ $analysis == "mayo" || $analysis == "external" || $analysis == "realignment" || $analysis == "variant" || $analysis == "realign-mayo" ]]
			then
				realign_dir=$output_realign/$sample
				variant_dir=$output_variant/$sample
				mkdir -p $realign_dir $variant_dir
				if [ $analysis == "variant" ]
				then
					infile=`cat $sample_info | grep -w ^BAM:${sample} | cut -d '=' -f2`
					num_bams=`echo $infile | tr " " "\n" | wc -l`
					if [ $num_bams -eq 0 ]
					then
						echo "sample info file is not properly configured"
						exit 1;
					fi	
					for ((i=1; i <=$num_bams; i++));
					do
						bam=`echo $infile | awk -v num=$i '{print $num}'`
						$samtools/samtools view -H $input/$bam 1>$realign_dir/$sample.$i.sorted.bam.header 2> $realign_dir/$sample.$i.sorted.bam.fix.log
						if [[ `cat $realign_dir/$sample.$i.sorted.bam.fix.log | wc -l` -gt 0 || `cat $realign_dir/$sample.$i.sorted.bam.header | wc -l` -le 0 ]]
						then
							$script_path/errorlog.sh $input/$bam whole_genome.sh ERROR "truncated or corrupted"
							exit 1;
						else
							rm $realign_dir/$sample.$i.sorted.bam.fix.log
						fi
						rm $realign_dir/$sample.$i.sorted.bam.header
						ln -s $input/$bam $realign_dir/$sample.$i.sorted.bam						
					done
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^reformat_BAM' | cut -d '=' -f2)
					qsub_args="-N $type.$version.reformat_BAM.$sample.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/reformat_BAM.sh $realign_dir $sample $run_info	
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^extract_reads_bam' | cut -d '=' -f2)
					qsub_args="-N $type.$version.extract_reads_bam.$sample.$run_num.$identify -hold_jid $type.$version.reformat_BAM.$sample.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/extract_reads_bam.sh $realign_dir $bamfile $run_info $igv
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^split_bam_chr' | cut -d '=' -f2)
					qsub_args="-N $type.$version.split_bam_chr.$sample.$run_num.$identify -hold_jid $type.$version.reformat_BAM.$sample.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/split_bam_chr.sh $realign_dir $sample $run_info
					variant_id="$type.$version.split_bam_chr.$sample.$run_num.$identify"
				else
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^realign_recal' | cut -d '=' -f2)
					qsub_args="-N $type.$version.realign_recal.$sample.$run_num.$identify -hold_jid $type.$version.processBAM.$sample.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/realign_recal.sh $align_dir $bamfile $sample $realign_dir $run_info 1	
					variant_id="$type.$version.realign_recal.$sample.$run_num.$identify"
				fi
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^igv_bam' | cut -d '=' -f2)
				qsub_args="-N $type.$version.igv_bam.$sample.$run_num.$identify -hold_jid $variant_id,$type.$version.extract_reads_bam.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/igv_bam.sh $output_realign $igv $sample $output_align $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^variants' | cut -d '=' -f2)
				qsub_args="-N $type.$version.variants.$sample.$run_num.$identify -hold_jid $variant_id -pe threaded $threads -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/variants.sh $realign_dir $sample $variant_dir 1 $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^merge_variant_single' | cut -d '=' -f2)
				qsub_args="-N $type.$version.merge_variant_single.$sample.$run_num.$identify -pe threaded $threads -hold_jid $type.$version.variants.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/merge_variant_single.sh $output_variant $sample $RSample $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^OnTarget_BAM' | cut -d '=' -f2)
				qsub_args="-N $type.$version.OnTarget_BAM.$sample.$run_num.$identify -hold_jid $variant_id -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args  $script_path/OnTarget_BAM.sh $realign_dir $output_OnTarget $sample $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^OnTarget_PILEUP' | cut -d '=' -f2)
				qsub_args="-N $type.$version.OnTarget_PILEUP.$sample.$run_num.$identify -hold_jid $variant_id -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/OnTarget_PILEUP.sh $realign_dir $output_OnTarget $sample $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^getCoverage' | cut -d '=' -f2)
				qsub_args="-N $type.$version.getCoverage.$sample.$run_num.$identify -hold_jid $type.$version.OnTarget_PILEUP.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/getCoverage.sh $output_OnTarget $numbers $sample $run_info    
			fi
			if [ $analysis == "ontarget" ]
			then
				if [ $variant_type == "BOTH" ]
				then
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^reformat_VARIANTs_OnTarget' | cut -d '=' -f2)
					qsub_args="-N $type.$version.reformat_VARIANTs_OnTarget.$sample.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/reformat_VARIANTs_OnTarget.sh $output_variant $RSample $sample $run_info 2
				elif [ $variant_type == "SNV" -o $variant_type == "INDEL" ]
				then
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^reformat_VARIANTs_OnTarget' | cut -d '=' -f2)
					qsub_args="-N $type.$version.reformat_VARIANTs_OnTarget.$sample.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/reformat_VARIANTs_OnTarget.sh $output_variant $RSample $sample $run_info 1
				fi
				hold_args="-hold_jid $type.$version.reformat_VARIANTs_OnTarget.$sample.$run_num.$identify"				
			elif [[ $analysis != "alignment" && $analysis != "annotation" ]]
			then
				hold_args="-hold_jid $type.$version.merge_variant_single.$sample.$run_num.$identify"	
			fi
			if [[ $analysis != "alignment" && $analysis != "annotation" ]]
			then
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^OnTarget_variant' | cut -d '=' -f2)
				qsub_args="-N $type.$version.OnTarget_variant.$sample.$run_num.$identify -t 1-$numchrs:1 $hold_args -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/OnTarget_variant.sh $output_variant $output_OnTarget $sample $run_info
			fi
			if [ $analysis == "annotation" ]
			then
				if [ $variant_type == "BOTH" ]
				then
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^reformat_VARIANTs' | cut -d '=' -f2)
					qsub_args="-N $type.$version.reformat_VARIANTs.$sample.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/reformat_VARIANTs.sh $output_OnTarget $sample $run_info 2
				elif [ $variant_type == "SNV" -o $variant_type == "INDEL" ]
				then
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^reformat_VARIANTs' | cut -d '=' -f2)
					qsub_args="-N $type.$version.reformat_VARIANTs.$sample.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/reformat_VARIANTs.sh $output_OnTarget $sample $run_info 1
				fi
				hold_args="-hold_jid $type.$version.reformat_VARIANTs.$sample.$run_num.$identify"
			elif [ $analysis != "alignment" ]
			then
				hold_args="-hold_jid $type.$version.OnTarget_variant.$sample.$run_num.$identify"
			fi
			if [ $analysis != "alignment" ]
			then
				if [ $variant_type == "SNV" -o $variant_type == "BOTH" ]
				then
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^sift' | cut -d '=' -f2)
					qsub_args="-N $type.$version.sift.$sample.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/sift.sh $sift $output_OnTarget $sample $run_info germline
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^polyphen' | cut -d '=' -f2)
					qsub_args="-N $type.$version.polyphen.$sample.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/polyphen.sh $polyphen $output_OnTarget $sample $run_info germline	    	
				fi
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^snpeff' | cut -d '=' -f2)
				qsub_args="-N $type.$version.snpeff.$sample.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/snpeff.sh $snpeff $output_OnTarget $sample $run_info germline		
				if [ $variant_type == "SNV" -o $variant_type == "BOTH" ]
				then
					hold="-hold_jid $type.$version.sift.$sample.$run_num.$identify,$type.$version.polyphen.$sample.$run_num.$identify,$type.$version.snpeff.$sample.$run_num.$identify"
				else
					hold="-hold_jid $type.$version.snpeff.$sample.$run_num.$identify"
				fi	
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^reports' | cut -d '=' -f2)
				qsub_args="-N $type.$version.reports.$sample.$run_num.$identify $hold -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/reports.sh $run_info $sample $TempReports $output_OnTarget $sift $snpeff $polyphen $output_dir germline
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^sample_report' | cut -d '=' -f2)
				qsub_args="-N $type.$version.sample_report.$sample.$run_num.$identify -hold_jid $type.$version.reports.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/sample_report.sh $output_dir $TempReports $sample $run_info germline
				if [[ $tool == "whole_genome"  && $analysis != "annotation" ]]
				then
					crest=$output_dir/struct/crest
					break=$output_dir/struct/break
					cnv=$output_dir/cnv/$sample
					mkdir -p $break $crest $cnv
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^run_single_crest' | cut -d '=' -f2)
					if [ $analysis == "variant" ]
					then
						qsub_args="-N $type.$version.run_single_crest.$sample.$run_num.$identify -hold_jid $type.$version.reformat_BAM.$sample.$run_num.$identify -t 1-$numchrs:1 -pe threaded 2 -l h_vmem=$mem"
						qsub $args $qsub_args $script_path/run_single_crest.sh $sample $realign_dir $bamfile $crest $run_info
					else
						qsub_args="-N $type.$version.run_single_crest.$sample.$run_num.$identify -hold_jid $type.$version.processBAM.$sample.$run_num.$identify -t 1-$numchrs:1 -pe threaded 2 -l h_vmem=$mem"
						qsub $args $qsub_args $script_path/run_single_crest.sh $sample $align_dir $bamfile $crest $run_info
					fi
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^run_cnvnator' | cut -d '=' -f2)
					qsub_args="-N $type.$version.run_cnvnator.$sample.$run_num.$identify -hold_jid $variant_id -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/run_cnvnator.sh $sample $realign_dir $cnv $run_info
					let nump=$numchrs+1;
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^run_breakdancer' | cut -d '=' -f2)
					qsub_args="-N $type.$version.run_breakdancer.$sample.$run_num.$identify -hold_jid $variant_id -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/run_breakdancer.sh $sample $output_realign $break $run_info
					$script_path/check_qstat.sh $limit
					qsub_args="-N $type.$version.run_breakdancer_in.$sample.$run_num.$identify -hold_jid $type.$version.igv_bam.$sample.$run_num.$identify -t $nump-$nump:$nump -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/run_breakdancer.sh $sample $igv $break $run_info
					### merge the structural variants
					hold="-hold_jid $type.$version.run_single_crest.$sample.$run_num.$identify,$type.$version.run_cnvnator.$sample.$run_num.$identify,"
					hold=$hold"$type.$version.run_breakdancer.$sample.$run_num.$identify,$type.$version.run_breakdancer_in.$sample.$run_num.$identify"
					mkdir -p $sv
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^summaryze_struct_single' | cut -d '=' -f2)
					qsub_args="-N $type.$version.summaryze_struct_single.$sample.$run_num.$identify -l h_vmem=$mem $hold"
					qsub $args $qsub_args $script_path/summaryze_struct_single.sh $sample $output_dir $run_info
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^plot_circos_cnv_sv' | cut -d '=' -f2)
					qsub_args="-N $type.$version.plot_circos_cnv_sv.$sample.$run_num.$identify -hold_jid $type.$version.summaryze_struct_single.$sample.$run_num.$identify -l h_vmem=$mem"
					break_file=$break/$sample/$sample.break
					crest_file=$crest/$sample/$sample.filter.crest
					cnv_file=$cnv/$sample.cnv.filter.bed
					qsub $args $qsub_args $script_path/plot_circos_cnv_sv.sh $break_file $crest_file $cnv_file $sample $circos $run_info	
				fi
				if [[ $tool == "whole_genome" && $analysis != "alignment" && $analysis != "annotation" && $analysis != "ontarget" ]]
				then
					mkdir -p $output_dir/Reports_per_Sample/ANNOT
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^annotation_CNV' | cut -d '=' -f2)
					qsub_args="-N $type.$version.annotation_CNV.$sample.$run_num.$identify -hold_jid $type.$version.plot_circos_cnv_sv.$sample.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/annotation_CNV.sh $sv $run_info $annot $sample
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^annotation_SV' | cut -d '=' -f2)
					qsub_args="-N $type.$version.annotation_SV.sh.$sample.$run_num.$identify -hold_jid $type.$version.plot_circos_cnv_sv.$sample.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/annotation_SV.sh $output_dir $run_info $annot $sample
				fi	
			fi
			if [[ $analysis != "annotation" && $analysis != "alignment" ]]
			then
				if [ $tool == "whole_genome" ]
				then
					hold_args="-hold_jid $type.$version.plot_circos_cnv_sv.$sample.$run_num.$identify,$type.$version.sample_report.$sample.$run_num.$identify,"
					hold_args=$hold_args"$type.$version.annotation_CNV.$sample.$run_num.$identify,$type.$version.annotation_SV.sh.$sample.$run_num.$identify"
				else
					hold_args="-hold_jid $type.$version.sample_report.$sample.$run_num.$identify"
				fi
			elif [ $analysis == "annotation" -o $analysis == "ontarget" ]
			then
				hold_args="-hold_jid $type.$version.sample_report.$sample.$run_num.$identify"
			elif [ $analysis == "alignment" ]
			then
				hold_args="-hold_jid $type.$version.processBAM.$sample.$run_num.$identify"
			fi	
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^sample_numbers' | cut -d '=' -f2)
			qsub_args="-N $type.$version.sample_numbers.$sample.$run_num.$identify $hold_args -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/sample_numbers.sh $output_dir $sample $run_info $numbers
			if [ $analysis != "alignment" ]
			then
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^gene_summary' | cut -d '=' -f2)
				qsub_args="-N $type.$version.gene_summary.$sample.$run_num.$identify $hold_args -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/gene_summary.sh $output_dir $sample $run_info $RSample		
			fi
		done
		### concat raw varaints
		if [[  $tool == "exome"  && $all_sites == "YES" ]]
		then
			id=""
			for s in `echo $samples | tr ":" "\n"`
			do
				id=$id"$type.$version.variants.$s.$run_num.$identify,"
			done
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^merge_raw_variants' | cut -d '=' -f2)
			qsub_args="-N $type.$version.merge_raw_variants.$run_num.$identify -t 1-$numchrs:1 -hold_jid $id -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/merge_raw_variants.sh $output_dir $run_info
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^concat_raw_variants' | cut -d '=' -f2)
			qsub_args="-N $type.$version.concat_raw_variants.$run_num.$identify -hold_jid $type.$version.merge_raw_variants.$run_num.$identify -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/concat_raw_variants.sh $output_dir $run_info	    
		fi	
		if [ $analysis != "alignment" ]
		then
			id=""
			for s in `echo $samples | tr ":" "\n"`
			do
				id=$id"$type.$version.sample_report.$s.$run_num.$identify,"
			done 
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^merge_sample' | cut -d '=' -f2)
			qsub_args="-N $type.$version.merge_sample.$run_num.$identify -hold_jid $id -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/merge_sample.sh $output_dir $run_info    
		fi
		id_igv=""
		id_numbers=""
		id_gene_summary=""
		id_coverage=""
		id_reads=
		for s1 in `echo $samples | tr ":" "\n"`
		do
			id_igv=$id_igv"$type.$version.igv_bam.$s1.$run_num.$identify,"
			id_numbers=$id_numbers"$type.$version.sample_numbers.$s1.$run_num.$identify,"
			id_gene_summary=$id_gene_summary"$type.$version.gene_summary.$s1.$run_num.$identify,"
			id_coverage=$id_coverage"$type.$version.getCoverage.$s1.$run_num.$identify,"
			id_reads=$id_reads"$type.$version.extract_reads_bam.$s1.$run_num.$identify,"
		done 
		if [ $analysis == "alignment" ]
		then
			hold="-hold_jid $id_numbers,$id_gene_summary"
		elif [ $analysis == "annotation" -o $analysis == "ontarget" ]
		then
			hold="-hold_jid $type.$version.merge_sample.$run_num.$identify,$id_numbers,$id_gene_summary"
		elif [[ $analysis == "mayo" || $analysis == "external" || $analysis == "variant" || $analysis == "realign-mayo" || $analysis == "realignment" ]]
		then
			if [ $tool == "whole_genome" ]
			then
				hold="-hold_jid $id_coverage,$id_igv,$id_numbers,$id_gene_summary,$type.$version.merge_sample.$run_num.$identify,$type.$version.annotation.CNV.sh.$run_num.$identify,$type.$version.annotation.SV.sh.$run_num.$identify,$id_reads"
			elif [[  $tool == "exome"  && $all_sites == "YES" ]]
			then
				hold="-hold_jid $id_coverage,$id_igv,$id_numbers,$id_gene_summary,$type.$version.merge_sample.$run_num.$identify,$type.$version.concat_raw_variants.$run_num.$identify,$id_reads"
			elif  [[  $tool == "exome"  && $all_sites == "NO" ]]
			then
				hold="-hold_jid $id_coverage,$id_igv,$id_numbers,$id_gene_summary,$type.$version.merge_sample.$run_num.$identify,$id_reads"
			fi
		fi
		## generate html page for all the module
		$script_path/check_qstat.sh $limit
		mem=$( cat $memory_info | grep -w '^generate_html' | cut -d '=' -f2)
		qsub_args="-N $type.$version.generate_html.$run_num.$identify $hold -l h_vmem=$mem"
		qsub $args $qsub_args $script_path/generate_html.sh $output_dir $run_info
	else
		echo "Multi-sample"
		numgroups=$(cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
		for sample in `echo $samples | tr ":" "\n"`
		do
			bamfile=$sample.sorted.bam
			align_dir=$output_dir/alignment/$sample;
			mkdir -p $align_dir
			if [[ $analysis == "mayo" || $analysis == "external" ]]
			then
				if [ $paired == 1 ]
				then
					let numfiles=(`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" |wc -l`)/2
				else
					let numfiles=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" |wc -l`
				fi
				if [ $numfiles -eq 0 ]
				then
					echo "sample info file is not properly configured"
					exit 1;
				fi	
				if [ $aligner == "novoalign" ]
				then
					echo "novoalign is used as aligner"
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^align_novo' | cut -d '=' -f2)
					qsub_args="-N $type.$version.align_novo.$sample.$run_num.$identify -pe threaded $threads -t 1-$numfiles:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/align_novo.sh $sample $output_dir $run_info
					hold="$type.$version.align_novo.$sample.$run_num.$identify"
				elif [ $aligner == "bwa" ]
				then
					echo "bwa is used as aligner"
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^align_read_bwa' | cut -d '=' -f2)
					qsub_args="-N $type.$version.align_read_bwa.R1.$sample.$run_num.$identify -pe threaded $threads -t 1-$numfiles:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/align_read_bwa.sh $sample $output_dir 1 $run_info
					if [ $paired == 1 ]
					then
						$script_path/check_qstat.sh $limit
						mem=$( cat $memory_info | grep -w '^align_read_bwa' | cut -d '=' -f2)
						qsub_args="-N $type.$version.align_read_bwa.R2.$sample.$run_num.$identify -pe threaded $threads -t 1-$numfiles:1 -l h_vmem=$mem"
						qsub $args $qsub_args $script_path/align_read_bwa.sh $sample $output_dir 2 $run_info	
						hold="-hold_jid $type.$version.align_read_bwa.R2.$sample.$run_num.$identify,$type.$version.align_read_bwa.R1.$sample.$run_num.$identify"
					else
						hold="-hold_jid $type.$version.align_read_bwa.R1.$sample.$run_num.$identify"
					fi	
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^align_bwa' | cut -d '=' -f2)
					qsub_args="-N $type.$version.align_bwa.$sample.$run_num.$identify $hold -t 1-$numfiles:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/align_bwa.sh $sample $output_dir $run_info
					hold="$type.$version.align_bwa.$sample.$run_num.$identify"
				fi	    
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^processBAM' | cut -d '=' -f2)
				qsub_args="-N $type.$version.processBAM.$sample.$run_num.$identify -hold_jid $hold -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/processBAM.sh $align_dir $sample $run_info 
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^extract_reads_bam' | cut -d '=' -f2)
				qsub_args="-N $type.$version.extract_reads_bam.$sample.$run_num.$identify -hold_jid $type.$version.processBAM.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/extract_reads_bam.sh $align_dir $bamfile $run_info $igv		
			elif [[ $analysis == "realignment" || $analysis == "realign-mayo" ]]
			then
				infile=`cat $sample_info | grep -w ^BAM:${sample} | cut -d '=' -f2`
				num_bams=`echo $infile | tr " " "\n" | wc -l`
				if [ $num_bams -eq 0 ]
				then
					echo "sample info file is not properly configured"
					exit 1;
				fi	
				for ((i=1; i <=$num_bams; i++));
				do
					bam=`echo $infile | awk -v num=$i '{print $num}'`
					$samtools/samtools view -H $input/$bam 1>$align_dir/$sample.$i.sorted.header 2> $align_dir/$sample.$i.sorted.bam.log
					if [[ `cat $align_dir/$sample.$i.sorted.bam.log | wc -l` -gt 0 || `cat $align_dir/$sample.$i.sorted.header | wc -l` -le 0 ]]
					then
						echo "$input/$bam : bam file is truncated or corrupted" 	
						exit 1;
					else
						rm $align_dir/$sample.$i.sorted.bam.log
					fi
					rm $align_dir/$sample.$i.sorted.header
					ln -s $input/$bam $align_dir/$sample.$i.sorted.bam
				done
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^processBAM' | cut -d '=' -f2)
				qsub_args="-N $type.$version.processBAM.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/processBAM.sh $align_dir $sample $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^extract_reads_bam' | cut -d '=' -f2)
				qsub_args="-N $type.$version.extract_reads_bam.$sample.$run_num.$identify -hold_jid $type.$version.processBAM.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/extract_reads_bam.sh $align_dir $bamfile $run_info $igv
			fi
		done	
		if [ $analysis != "alignment" ]
		then
			for group in `echo $groups | tr ":" "\n"`
			do
				samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" "\n")
				bam_samples=""
				input_dirs=""
				names_samples=""
				for sample in $samples
				do
					bamfile=$sample.sorted.bam
					align_dir=$output_dir/alignment/$sample;	
					### setting the bams and its path for multiple sample anaylysis
					names_samples=$names_samples"$sample:"
					bam_samples=$bam_samples"$sample.sorted.bam:"
					input_dirs=$input_dirs"$output_dir/alignment/$sample:"
				done
				realign_dir=$output_dir/realign/$group
				variant_dir=$output_dir/variants/$group
				mkdir -p $realign_dir $variant_dir
				if [ $analysis == "variant" ]
				then
					vvid=""
                    infile=`cat $sample_info | grep -w ^BAM:${group} | cut -d '=' -f2`
					num_bams=`echo $infile | tr " " "\n" | wc -l`
					if [ $num_bams -eq 0 ]
					then
						echo "sample info file is not properly configured"
						exit 1;
					fi	
					for ((i=1; i <=$num_bams; i++));
					do
						bam=`echo $infile | awk -v num=$i '{print $num}'`
						$samtools/samtools view -H $input/$bam 1>$realign_dir/$group.$i.sorted.header 2> $realign_dir/$group.$i.sorted.bam.log
						if [[ `cat $realign_dir/$group.$i.sorted.bam.log | wc -l` -gt 0 || `cat $realign_dir/$group.$i.sorted.header | wc -l` -le 0 ]]
						then
							echo "$input/$bam : bam file is truncated or corrupted" 	
							exit 1;
						else
							rm $realign_dir/$group.$i.sorted.bam.log
						fi	
						rm $realign_dir/$group.$i.sorted.header
						ln -s $input/$bam $realign_dir/$group.$i.sorted.bam
					done
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^reformat_pairBAM' | cut -d '=' -f2)
					qsub_args="-N $type.$version.reformat_pairBAM.$group.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/reformat_pairBAM.sh $realign_dir $group $run_info
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^extract_reads_bam' | cut -d '=' -f2)
					qsub_args="-N $type.$version.extract_reads_bam.$group.$run_num.$identify -hold_jid $type.$version.reformat_pairBAM.$group.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/extract_reads_bam.sh $realign_dir $group.sorted.bam $run_info $igv $group
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^split_bam_chr' | cut -d '=' -f2)
					qsub_args="-N $type.$version.split_bam_chr.$group.$run_num.$identify -hold_jid $type.$version.reformat_pairBAM.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/split_bam_chr.sh $realign_dir $group $run_info
					variant_id="$type.$version.split_bam_chr.$group.$run_num.$identify"
					vvid="$type.$version.split_bam_chr.$group.$run_num.$identify"
				else        
					vvid=""
					for sample in $samples
					do
						vvid=$vvid"$type.$version.processBAM.$sample.$run_num.$identify,$type.$version.extract_reads_bam.$sample.$run_num.$identify"
					done    
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^realign_recal' | cut -d '=' -f2)
					qsub_args="-N $type.$version.realign_recal.$group.$run_num.$identify -hold_jid $vvid -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/realign_recal.sh $input_dirs $bam_samples $names_samples $realign_dir $run_info 1
					variant_id="$type.$version.realign_recal.$group.$run_num.$identify"
				fi
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^split_sample_pair' | cut -d '=' -f2)
				qsub_args="-N $type.$version.split_sample_pair.$group.$run_num.$identify -hold_jid $variant_id -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/split_sample_pair.sh $output_realign $igv $group $output_align $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^igv_bam' | cut -d '=' -f2)
				qsub_args="-N $type.$version.igv_bam.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/igv_bam.sh $output_realign $igv $group $output_align $run_info 
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^variants' | cut -d '=' -f2)
				qsub_args="-N $type.$version.variants.$group.$run_num.$identify -hold_jid $variant_id -pe threaded $threads -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/variants.sh $realign_dir $names_samples $variant_dir 1 $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^merge_variant_group' | cut -d '=' -f2)
				qsub_args="-N $type.$version.merge_variant_group.$group.$run_num.$identify -pe threaded $threads -hold_jid $type.$version.variants.$group.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/merge_variant_group.sh $output_variant $group $RSample $run_info 
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^OnTarget_BAM' | cut -d '=' -f2)
				qsub_args="-N $type.$version.OnTarget_BAM.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/OnTarget_BAM.sh $igv $output_OnTarget $group $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^OnTarget_PILEUP' | cut -d '=' -f2)
				qsub_args="-N $type.$version.OnTarget_PILEUP.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/OnTarget_PILEUP.sh $realign_dir $output_OnTarget $group $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^getCoverage' | cut -d '=' -f2)
				qsub_args="-N $type.$version.getCoverage.$group.$run_num.$identify -hold_jid $type.$version.OnTarget_PILEUP.$group.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/getCoverage.sh $output_OnTarget $numbers $group $run_info
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^OnTarget_variant' | cut -d '=' -f2)
				qsub_args="-N $type.$version.OnTarget_variant.$group.$run_num.$identify -hold_jid $type.$version.merge_variant_group.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/OnTarget_variant.sh $output_variant $output_OnTarget $group $run_info
				hold_args="-hold_jid $type.$version.OnTarget_variant.$group.$run_num.$identify"
				## SIFT
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^sift' | cut -d '=' -f2)
				qsub_args="-N $type.$version.sift.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/sift.sh $sift $output_OnTarget $group $run_info germline
				$script_path/check_qstat.sh $limit
				qsub_args="-N $type.$version.sift.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/sift.sh $sift $output_OnTarget $group $run_info somatic
				## SNPEFF
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^snpeff' | cut -d '=' -f2)
				qsub_args="-N $type.$version.snpeff.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/snpeff.sh $snpeff $output_OnTarget $group $run_info germline
				$script_path/check_qstat.sh $limit
				qsub_args="-N $type.$version.snpeff.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/snpeff.sh $snpeff $output_OnTarget $group $run_info somatic
				##POLYPHEN
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^polyphen' | cut -d '=' -f2)
				qsub_args="-N $type.$version.polyphen.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/polyphen.sh $polyphen $output_OnTarget $group $run_info germline  
				$script_path/check_qstat.sh $limit
				qsub_args="-N $type.$version.polyphen.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/polyphen.sh $polyphen $output_OnTarget $group $run_info somatic
				hold="$type.$version.sift.$group.$run_num.$identify,$type.$version.snpeff.$group.$run_num.$identify,$type.$version.polyphen.$group.$run_num.$identify"
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^reports' | cut -d '=' -f2)
				qsub_args="-N $type.$version.reports.$group.$run_num.$identify -hold_jid $hold -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/reports.sh $run_info $group $TempReports $output_OnTarget $sift $snpeff $polyphen $output_dir germline
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^sample_report' | cut -d '=' -f2)
				qsub_args="-N $type.$version.sample_reports.$sample.$run_num.$identify -hold_jid $type.$version.reports.$group.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/sample_report.sh $output_dir $TempReports $group $run_info germline
				hold="$type.$version.sift.${group}.$run_num.$identify,$type.$version.snpeff.$group.$run_num.$identify,$type.$version.polyphen.$group.$run_num.$identify"
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^reports' | cut -d '=' -f2)
				qsub_args="-N $type.$version.reports.$group.$run_num.$identify -hold_jid $hold -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/reports.sh $run_info $group $TempReports $output_OnTarget $sift $snpeff $polyphen $output_dir somatic
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^sample_report' | cut -d '=' -f2)
				qsub_args="-N $type.$version.sample_report.$group.$run_num.$identify -hold_jid $type.$version.reports.$group.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/sample_report.sh $output_dir $TempReports $group $run_info somatic
				if [ $tool == "whole_genome" ]
				then
					crest=$output_dir/struct/crest
					break=$output_dir/struct/break
					mkdir -p $break
					mkdir -p $crest
					id=""
					for sam in `cat $sample_info| grep -w "^$group" | cut -d '=' -f2`
					do
						$script_path/check_qstat.sh $limit
						mem=$( cat $memory_info | grep -w '^run_crest_multi_cover' | cut -d '=' -f2)
						if [ $analysis == "variant" ]
						then
							qsub_args="-N $type.$version.run_crest_multi_cover.$group.$sam.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
							qsub $args $qsub_args $script_path/run_crest_multi_cover.sh $sam $group $igv $crest $run_info
						else
							qsub_args="-N $type.$version.run_crest_multi_cover.$group.$sam.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
							qsub $args $qsub_args $script_path/run_crest_multi_cover.sh $sam $group $output_align/$sam/ $crest $run_info
						fi	
						id=$id"$type.$version.run_crest_multi_cover.$group.$sam.$run_num.$identify,"
					done
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^run_crest_multi' | cut -d '=' -f2)
					qsub_args="-N $type.$version.run_crest_multi.$group.$run_num.$identify -hold_jid $id -t 1-$numchrs:1 -pe threaded 2 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/run_crest_multi.sh $group $igv $crest $run_info
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^run_segseq' | cut -d '=' -f2)
					qsub_args="-N $type.$version.run_segseq.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l matlab_lic=1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/run_segseq.sh $group $igv $cnv $run_info    
					let nump=$numchrs+1;    
					mkdir -p $break/$group
					id=""
					for sam in `cat $sample_info| grep -w "^$group" | cut -d '=' -f2`
					do
						$script_path/check_qstat.sh $limit
						mem=$( cat $memory_info | grep -w '^run_breakdancer' | cut -d '=' -f2)
						qsub_args="-N $type.$version.run_breakdancer.$group.$sam.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
						qsub $args $qsub_args $script_path/run_breakdancer.sh $sam $igv $break/$group $run_info $group
						$script_path/check_qstat.sh $limit
						qsub_args="-N $type.$version.run_breakdancer_in.$group.$sam.$run_num.$identify -hold_jid $type.$version.igv_bam.$group.$run_num.$identify -t $nump-$nump:$nump -l h_vmem=$mem"
						qsub $args $qsub_args $script_path/run_breakdancer.sh $sam $igv $break/$group $run_info $group
						id=$id"$type.$version.run_breakdancer.$group.$sam.$run_num.$identify,$type.$version.run_breakdancer_in.$group.$sam.$run_num.$identify,"
					done
					hhold="$id,$type.$version.run_segseq.$group.$run_num.$identify,$type.$version.run_crest_multi.$group.$run_num.$identify"
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^summaryze_struct_group' | cut -d '=' -f2)
					qsub_args="-N $type.$version.summaryze_struct_group.$group.$run_num.$identify -hold_jid $hhold -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/summaryze_struct_group.sh $group $output_dir $run_info
					mkdir -p $output_dir/circos;
					samm=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" "\n")
					sampleArray=()
					i=1
					for sample in $samm
					do
						sampleArray[$i]=$sample
						let i=i+1
					done
					for i in $(seq 2 ${#sampleArray[@]})
					do  
						tumor=${sampleArray[$i]}
						$script_path/check_qstat.sh $limit
						mem=$( cat $memory_info | grep -w '^plot_circos_cnv_sv' | cut -d '=' -f2)
						qsub_args="-N $type.$version.plot_circos_cnv_sv.$group.$tumor.$i.$run_num.$identify -hold_jid $type.$version.summaryze_struct_group.$group.$run_num.$identify -l h_vmem=$mem"
						break_file=$struct/$group.$tumor.somatic.break
						crest_file=$struct/$group.$tumor.somatic.filter.crest
						cnv_file=$cnv/$group/$tumor.cnv.filter.bed
						qsub $args $qsub_args $script_path/plot_circos_cnv_sv.sh $break_file $crest_file $cnv_file $group.$tumor $circos $run_info
					done
				fi
			done
			mkdir -p $annot
			if [ $tool == "whole_genome" ]
			then
				id=""
				for group in `echo $groups | tr ":" "\n"`
				do
					samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" "\n")
					sampleArray=()
					i=1
					for sample in $samples
					do
						sampleArray[$i]=$sample
						let i=i+1
					done
					for i in $(seq 2 ${#sampleArray[@]})
					do
						tumor=${sampleArray[$i]}
						id=$id"$type.$version.plot_circos_cnv_sv.$group.$tumor.$i.$run_num.$identify,"
					done
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^annotation_CNV' | cut -d '=' -f2)
					qsub_args="-N $type.$version.annotation_CNV.$group.$run_num.$identify -hold_jid $id -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/annotation_CNV.sh $sv $run_info $annot $group
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^annotation_SV' | cut -d '=' -f2)
					qsub_args="-N $type.$version.annotation_SV.$group.$run_num.$identify -hold_jid $id -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/annotation_SV.sh $output_dir $run_info $annot $group
				done
			fi
			### generate reports for all the samples
			id=""
			for group in `echo $groups | tr ":" "\n"`
			do
				id=$id"$type.$version.sample_report.$group.$run_num.$identify,"	
			done
			for group in `echo $groups | tr ":" "\n"`
			do
				samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" " ")
				for sample in $samples
				do
					id=$id"$type.$version.extract_reads_bam.$sample.$run_num.$identify,"
				done
			done

			if [ $tool == "whole_genome" ]
			then
				for group in `echo $groups | tr ":" "\n"`
				do
					id=$id"$type.$version.summaryze_struct_group.$group.$run_num.$identify,$type.$version.annotation_CNV.$group.$run_num.$identify,$type.$version.annotation_SV.$group.$run_num.$identify,"
				done
			fi    
			for group in `echo $groups | tr ":" "\n"`
			do
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^sample_numbers' | cut -d '=' -f2)
				qsub_args="-N $type.$version.sample_numbers.$group.$run_num.$identify -hold_jid $id -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/sample_numbers.sh $output_dir $group $run_info $numbers
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^gene_summary' | cut -d '=' -f2)
				qsub_args="-N $type.$version.gene_summary.$group.$run_num.$identify -hold_jid $id -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/gene_summary.sh $output_dir $group $run_info $RSample
			done

			if [ $tool == "exome" ]
			then
				id=""
				for group in `echo $groups | tr ":" "\n"`
				do
					id=$id"$type.$version.getCoverage.$group.$run_num.$identify,$type.$version.sample_numbers.$group.$run_num.$identify,$type.$version.gene_summary.$group.$run_num.$identify,$type.$version.igv_bam.$group.$run_num.$identify,$type.$version.sample_report.$group.$run_num.$identify,$type.$version.annotation_CNV.$run_num.$identify,$type.$version.annotation_SV.$run_num.$identify,"
				done          
			elif [ $tool == "whole_genome" ]
			then
				id=""
				for group in `echo $groups | tr ":" "\n"`
				do
					id=$id"$type.$version.getCoverage.$group.$run_num.$identify,$type.$version.sample_numbers.$group.$run_num.$identify,$type.$version.gene_summary.$group.$run_num.$identify,$type.$version.igv_bam.$group.$run_num.$identify,$type.$version.sample_report.$group.$run_num.$identify,"
				done    
			fi
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^merge_sample' | cut -d '=' -f2)
			qsub_args="-N $type.$version.merge_sample.$run_num.$identify -hold_jid $id -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/merge_sample.sh $output_dir $run_info	
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^generate_html' | cut -d '=' -f2)
			qsub_args="-N $type.$version.generate_html.$run_num.$identify -hold_jid $type.$version.merge_sample.$run_num.$identify -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/generate_html.sh $output_dir $run_info 	
		fi
	fi
	echo `date`
fi


#### end of the script
