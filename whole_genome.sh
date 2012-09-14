#!/bin/bash

########################################################
###### 	MASTER SCRIPT FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			whole_genome_pipeline.sh
######		Date:				09/07/2012
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
	if [ `perl -v 2>&1 | tr " " "\n" | grep '^v'` != "v5.10.0" ]
	then
		echo -e "\nperl path is not correct in your enviornment"
		echo "Perl path should point to /usr/local/biotools/perl/5.10.0/bin/perl if the user use the command which perl and user can change this using mayobiotools"
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
	cat $run_info | sed 's/^[ \t]*//;s/[ \t]*$//' > $run_info.tmp
	mv $run_info.tmp $run_info
	#### check for unique identification number
	identify=$( cat $run_info | grep -w '^IDENTIFICATION_NUMBER' | cut -d '=' -f2)
	if echo $identify | egrep -q '^[0-9]+$'
	then
		echo "HURRAY !!! you are good to run the workflow"
	else
		echo "ERROR : unique identification for the workflow was not generated, please run the unique_id.sh script to generate the same before running the workflow script."
		exit 1;
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
	#### extract paths
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
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
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
	limit=$( cat $tool_info | grep -w '^JOB_LIMIT' | cut -d '=' -f2 )
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
	if [ -f $output_dir/folder_exist.log ]
	then
		echo "ERROR: folder already exist"
		exit 1;
	fi	
	## copy cofig files
	$script_path/copy_config.sh $output_dir $run_info
	run_info=$output_dir/run_info.txt
	cat $run_info | grep -w -v -E '^TOOL_INFO|^SAMPLE_INFO' > $run_info.tmp
	echo -e "TOOL_INFO=$output_dir/tool_info.txt\nSAMPLE_INFO=$output_dir/sample_info.txt" | cat $run_info.tmp - > $run_info
	rm $run_info.tmp
	tool_info=$output_dir/tool_info.txt
	sample_info=$output_dir/sample_info.txt
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
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.align_novo.$sample.$run_num -l h_vmem=4G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_novo.sh $sample $output_dir $run_info
				elif [ $aligner == "bwa" ]
				then
					echo "bwa is used as aligner"
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.align_read_bwa.R1.$sample.$run_num -l h_vmem=1G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_read_bwa.sh $sample $output_dir 1 $run_info
					if [ $paired == 1 ]
					then
						$script_path/check_qstat.sh $limit
						qsub $args -N $type.$version.align_read_bwa.R2.$sample.$run_num -l h_vmem=1G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_read_bwa.sh $sample $output_dir 2 $run_info
						hold="-hold_jid $type.$version.align_read_bwa.R1.$sample.$run_num,$type.$version.align_read_bwa.R2.$sample.$run_num"
					else
						hold="-hold_jid $type.$version.align_read_bwa.R1.$sample.$run_num"
					fi	
					$script_path/check_qstat.sh $limit
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
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.processBAM.$sample.$run_num -pe threaded $threads -l h_vmem=4G $hold $script_path/processBAM.sh $align_dir $sample $run_info 	
				if [ $analysis != "alignment" ]
				then
					$script_path/check_qstat.sh $limit
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
					$samtools/samtools view -H $input/$bam 1>$input/$bam.header 2> $align_dir/$sample.$i.sorted.bam.log
					if [ `cat $align_dir/$sample.$i.sorted.bam.log | wc -l` -gt 0 ]
					then
						echo "$input/$bam : bam file is truncated or corrupted" 	
						exit 1;
					else
						rm $align_dir/$sample.$i.sorted.bam.log
					fi	
					rm $input/$bam.header
					ln -s $input/$bam $align_dir/$sample.$i.sorted.bam
					$script_path/dashboard.sh $sample $run_info Beginning started $i
				done  
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.processBAM.$sample.$run_num -pe threaded $threads -l h_vmem=4G $script_path/processBAM.sh $align_dir $sample $run_info
				$script_path/check_qstat.sh $limit
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
						$samtools/samtools view -H $input/$bam 1>$input/$bam.header 2> $align_dir/$sample.$i.sorted.bam.fix.log
						if [ `cat $align_dir/$sample.$i.sorted.bam.fix.log | wc -l` -gt 0 ]
						then
							echo "$input/$bam : bam file is truncated or corrupted" 	
							exit 1;
						else
							rm $align_dir/$sample.$i.sorted.bam.fix.log
						fi
						rm $input/$bam.header
						ln -s $input/$bam $realign_dir/$sample.$i.sorted.bam						
					done
					
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.reformat_BAM.$sample.$run_num -l h_vmem=8G $script_path/reformat_BAM.sh $realign_dir $sample $run_info	
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.extract_reads_bam.$sample.$run_num -l h_vmem=8G -hold_jid $type.$version.reformat_BAM.$sample.$run_num $script_path/extract_reads_bam.sh $realign_dir $bamfile $run_info $output_dir/IGV_BAM
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.split_bam_chr.$sample.$run_num -hold_jid $type.$version.reformat_BAM.$sample.$run_num -l h_vmem=2G -t 1-$numchrs:1 $script_path/split_bam_chr.sh $realign_dir $sample $run_info
					variant_id="$type.$version.split_bam_chr.$sample.$run_num"
				else
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.realign_recal.$sample.$run_num -hold_jid $type.$version.processBAM.$sample.$run_num -l h_vmem=8G -t 1-$numchrs:1 $script_path/realign_recal.sh $align_dir $bamfile $sample $realign_dir $run_info 1	
					variant_id="$type.$version.realign_recal.$sample.$run_num"
				fi
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.igv_bam.$sample.$run_num -l h_vmem=2G -hold_jid $variant_id,$type.$version.extract_reads_bam.$sample.$run_num $script_path/igv_bam.sh $output_dir/realign $output_dir/IGV_BAM $sample $output_dir/alignment $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.variants.$sample.$run_num -hold_jid $variant_id -pe threaded $threads -l h_vmem=4G -t 1-$numchrs:1 $script_path/variants.sh $realign_dir $sample $variant_dir 1 $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.merge_variant_single.$sample.$run_num -l h_vmem=4G -pe threaded $threads -hold_jid $type.$version.variants.$sample.$run_num $script_path/merge_variant_single.sh $output_dir/variants $sample $output_dir/Reports_per_Sample/ $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.OnTarget_BAM.$sample.$run_num -hold_jid $variant_id -l h_vmem=3G -t 1-$numchrs:1 $script_path/OnTarget_BAM.sh $realign_dir $output_dir/OnTarget $sample $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.OnTarget_PILEUP.$sample.$run_num -hold_jid $variant_id -l h_vmem=6G -t 1-$numchrs:1 $script_path/OnTarget_PILEUP.sh $realign_dir $output_dir/OnTarget $sample $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.getCoverage.$sample.$run_num -hold_jid $type.$version.OnTarget_PILEUP.$sample.$run_num -l h_vmem=2G $script_path/getCoverage.sh $output_dir/OnTarget $output_dir/numbers $sample $run_info    
			fi
			if [ $analysis == "ontarget" ]
			then
				if [ $variant_type == "BOTH" ]
				then
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.reformat_VARIANTs_OnTarget.$sample.$run_num -l h_vmem=4G $script_path/reformat_VARIANTs_OnTarget.sh $output_dir/variants $output_dir/Reports_per_Sample $sample $run_info 2
				elif [ $variant_type == "SNV" -o $variant_type == "INDEL" ]
				then
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.reformat_VARIANTs_OnTarget.$sample.$run_num -l h_vmem=4G $script_path/reformat_VARIANTs_OnTarget.sh $output_dir/variants $output_dir/Reports_per_Sample $sample $run_info 1
				fi
				hold_args="-hold_jid $type.$version.reformat_VARIANTs_OnTarget.$sample.$run_num"				
			elif [[ $analysis != "alignment" && $analysis != "annotation" ]]
			then
				hold_args="-hold_jid $type.$version.merge_variant_single.$sample.$run_num"	
			fi
			if [[ $analysis != "alignment" && $analysis != "annotation" ]]
			then
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.OnTarget_variant.$sample.$run_num -t 1-$numchrs:1 $hold_args -l h_vmem=2G $script_path/OnTarget_variant.sh $output_dir/variants $output_dir/OnTarget $sample $run_info
			fi
			if [ $analysis == "annotation" ]
			then
				if [ $variant_type == "BOTH" ]
				then
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.reformat_VARIANTs.$sample.$run_num -l h_vmem=2G $script_path/reformat_VARIANTs.sh $output_OnTarget $sample $run_info 2
				elif [ $variant_type == "SNV" -o $variant_type == "INDEL" ]
				then
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.reformat_VARIANTs.$sample.$run_num -l h_vmem=2G $script_path/reformat_VARIANTs.sh $output_OnTarget $sample $run_info 1
				fi
				hold_args="-hold_jid $type.$version.reformat_VARIANTs.$sample.$run_num"
			elif [ $analysis != "alignment" ]
			then
				hold_args="-hold_jid $type.$version.OnTarget_variant.$sample.$run_num"
			fi
			if [ $analysis != "alignment" ]
			then
				if [ $variant_type == "SNV" -o $variant_type == "BOTH" ]
				then
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.sift.$sample.$run_num $hold_args -t 1-$numchrs:1 -l h_vmem=4G $script_path/sift.sh $sift $output_OnTarget $sample $run_info
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.polyphen.$sample.$run_num $hold_args -pe threaded $threads -t 1-$numchrs:1 -l h_vmem=4G $script_path/polyphen.sh $polyphen $output_OnTarget $sample $run_info	    	
				fi
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.snpeff.$sample.$run_num $hold_args -t 1-$numchrs:1 -l h_vmem=4G $script_path/snpeff.sh $snpeff $output_OnTarget $sample $run_info		
				if [ $variant_type == "SNV" -o $variant_type == "BOTH" ]
				then
					hold="-hold_jid $type.$version.sift.$sample.$run_num,$type.$version.polyphen.$sample.$run_num,$type.$version.snpeff.$sample.$run_num"
				else
					hold="-hold_jid $type.$version.snpeff.$sample.$run_num"
				fi	
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.sample_reports.$sample.$run_num $hold -t 1-$numchrs:1 -l h_vmem=4G $script_path/sample_reports.sh $run_info $sample $TempReports $output_OnTarget $sift $snpeff $polyphen $output_dir
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.sample_report.$sample.$run_num -hold_jid $type.$version.sample_reports.$sample.$run_num -l h_vmem=2G $script_path/sample_report.sh $output_dir $TempReports $sample $run_info
				if [[ $tool == "whole_genome"  && $analysis != "annotation" ]]
				then
					crest=$output_dir/struct/crest
					break=$output_dir/struct/break
					cnv=$output_dir/cnv/$sample
					mkdir -p $break $crest $cnv
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.run_single_crest.sh.$sample.$run_num -hold_jid $variant_id -t 1-$numchrs:1 -l h_vmem=6G $script_path/run_single_crest.sh $sample $realign_dir $crest $run_info
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.run_cnvnator.$sample.$run_num -hold_jid $variant_id -l h_vmem=3G -t 1-$numchrs:1 $script_path/run_cnvnator.sh $sample $realign_dir $cnv $run_info
					let nump=$numchrs+1;
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.run_breakdancer.$sample.$run_num -hold_jid $variant_id -l h_vmem=4G -t 1-$numchrs:1 $script_path/run_breakdancer.sh $sample $output_dir/realign $break $run_info
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.run_breakdancer_in.$sample.$run_num -hold_jid $type.$version.igv_bam.$sample.$run_num -l h_vmem=4G -t $nump-$nump:$nump $script_path/run_breakdancer.sh $sample $output_dir/IGV_BAM $break $run_info
					### merge the structural variants
					hold="-hold_jid $type.$version.run_single_crest.sh.$sample.$run_num,$type.$version.run_cnvnator.$sample.$run_num,$type.$version.run_breakdancer.$sample.$run_num,$type.$version.run_breakdancer_in.$sample.$run_num"
					mkdir -p $output_dir/Reports_per_Sample/SV
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.summaryze_struct_single.$sample.$run_num -l h_vmem=5G $hold $script_path/summaryze_struct_single.sh $sample $output_dir $run_info
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.plot_circos_cnv_sv.$sample.$run_num -hold_jid $type.$version.summaryze_struct_single.$sample.$run_num -l h_vmem=2G $script_path/plot_circos_cnv_sv.sh $break/$sample/$sample.break $crest/$sample/$sample.filter.crest $cnv/$sample.cnv.filter.bed $sample $output_dir/circos $run_info	
				fi
				if [[ $tool == "whole_genome" && $analysis != "alignment" && $analysis != "annotation" && $analysis != "ontarget" ]]
				then
					mkdir -p $output_dir/Reports_per_Sample/ANNOT
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.annotation_CNV.$sample.$run_num -l h_vmem=2G -hold_jid $type.$version.plot_circos_cnv_sv.$sample.$run_num $script_path/annotation_CNV.sh $output_dir/Reports_per_Sample/SV/ $run_info $output_dir/Reports_per_Sample/ANNOT/ $sample
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.annotation_SV.sh.$sample.$run_num -l h_vmem=2G -hold_jid $type.$version.plot_circos_cnv_sv.$sample.$run_num $script_path/annotation_SV.sh $output_dir $run_info $output_dir/Reports_per_Sample/ANNOT/ $sample
				fi	
			fi
			if [[ $analysis != "annotation" && $analysis != "alignment" ]]
			then
				if [ $tool == "whole_genome" ]
				then
					hold_args="-hold_jid $type.$version.plot_circos_cnv_sv.$sample.$run_num,$type.$version.sample_report.$sample.$run_num,$type.$version.annotation_CNV.$sample.$run_num,$type.$version.annotation_SV.sh.$sample.$run_num"
				else
					hold_args="-hold_jid $type.$version.sample_report.$sample.$run_num"
				fi
			elif [ $analysis == "annotation" -o $analysis == "ontarget" ]
			then
				hold_args="-hold_jid $type.$version.sample_report.$sample.$run_num"
			elif [ $analysis == "alignment" ]
			then
				hold_args="-hold_jid $type.$version.processBAM.$sample.$run_num"
			fi	
			$script_path/check_qstat.sh $limit
			qsub $args -N $type.$version.sample_numbers.$sample.$run_num $hold_args -l h_vmem=2G $script_path/sample_numbers.sh $output_dir $sample $run_info $output_dir/numbers
			if [ $analysis != "alignment" ]
			then
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.gene_summary.$sample.$run_num $hold_args -l h_vmem=2G $script_path/gene_summary.sh $output_dir $sample $run_info $output_dir/Reports_per_Sample		
			fi
		done
		### concat raw varaints
		if [[  $tool == "exome"  && $all_sites == "YES" ]]
		then
			id=""
			for s in `echo $samples | tr ":" "\n"`
			do
				id=$id"$type.$version.variants.$s.$run_num,"
			done
			$script_path/check_qstat.sh $limit
			qsub $args -N $type.$version.merge_raw_variants.$run_num -t 1-$numchrs:1 -hold_jid $id -l h_vmem=2G $script_path/merge_raw_variants.sh $output_dir $run_info
			$script_path/check_qstat.sh $limit
			qsub $args -N $type.$version.concat_raw_variants.$run_num -hold_jid $type.$version.merge_raw_variants.$run_num -l h_vmem=2G $script_path/concat_raw_variants.sh $output_dir $run_info	    
		fi	
		if [ $analysis != "alignment" ]
		then
			id=""
			for s in `echo $samples | tr ":" "\n"`
			do
				id=$id"$type.$version.sample_report.$s.$run_num,"
			done
			if [ $tool == "exome" ]
			then
				mem="-l h_vmem=5G"
			else
				mem="-l h_vmem=16G"
			fi    
			$script_path/check_qstat.sh $limit
			qsub $args -N $type.$version.merge_sample.$run_num -hold_jid $id $mem $script_path/merge_sample.sh $output_dir $run_info    
		fi
		id_igv=""
		id_numbers=""
		id_gene_summary=""
		id_coverage=""
		id_reads=
		for s1 in `echo $samples | tr ":" "\n"`
		do
			id_igv=$id_igv"$type.$version.igv_bam.$s1.$run_num,"
			id_numbers=$id_numbers"$type.$version.sample_numbers.$s1.$run_num,"
			id_gene_summary=$id_gene_summary"$type.$version.gene_summary.$s1.$run_num,"
			id_coverage=$id_coverage"$type.$version.getCoverage.$s1.$run_num,"
			id_reads=$id_reads"$type.$version.extract_reads_bam.$s1.$run_num,"
		done 
		if [ $analysis == "alignment" ]
		then
			hold="-hold_jid $id_numbers,$id_gene_summary"
		elif [ $analysis == "annotation" -o $analysis == "ontarget" ]
		then
			hold="-hold_jid $type.$version.merge_sample.$run_num,$id_numbers,$id_gene_summary"
		elif [[ $analysis == "mayo" || $analysis == "external" || $analysis == "variant" || $analysis == "realign-mayo" || $analysis == "realignment" ]]
		then
			if [ $tool == "whole_genome" ]
			then
				hold="-hold_jid $id_coverage,$id_igv,$id_numbers,$id_gene_summary,$type.$version.merge_sample.$run_num,$type.$version.annotation.CNV.sh.$run_num,$type.$version.annotation.SV.sh.$run_num,$id_reads"
			elif [[  $tool == "exome"  && $all_sites == "YES" ]]
			then
				hold="-hold_jid $id_coverage,$id_igv,$id_numbers,$id_gene_summary,$type.$version.merge_sample.$run_num,$type.$version.concat_raw_variants.$run_num,$id_reads"
			elif  [[  $tool == "exome"  && $all_sites == "NO" ]]
			then
				hold="-hold_jid $id_coverage,$id_igv,$id_numbers,$id_gene_summary,$type.$version.merge_sample.$run_num,$id_reads"
			fi
		fi
		## generate html page for all teh modules
		if [[ $upload_tb == "YES" ]]
		then
			mem="8G"
		else
			mem="2G"
		fi	
		$script_path/check_qstat.sh $limit
		qsub $args -l h_vmem=$mem -N $type.$version.generate_html.$run_num $hold $script_path/generate_html.sh $output_dir $run_info
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
				if [ $aligner == "novoalign" ]
				then
					echo "novoalign is used as aligner"
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.align_novo.$sample.$run_num -l h_vmem=4G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_novo.sh $sample $output_dir $run_info
					hold="$type.$version.align_novo.$sample.$run_num"
				elif [ $aligner == "bwa" ]
				then
					echo "bwa is used as aligner"
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.align_read_bwa.R1.$sample.$run_num -l h_vmem=1G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_read_bwa.sh $sample $output_dir 1 $run_info
					if [ $paired == 1 ]
					then
						$script_path/check_qstat.sh $limit
						qsub $args -N $type.$version.align_read_bwa.R2.$sample.$run_num -l h_vmem=1G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_read_bwa.sh $sample $output_dir 2 $run_info	
						hold="-hold_jid $type.$version.align_read_bwa.R2.$sample.$run_num,$type.$version.align_read_bwa.R1.$sample.$run_num"
					else
						hold="-hold_jid $type.$version.align_read_bwa.R1.$sample.$run_num"
					fi	
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.align_bwa.$sample.$run_num -l h_vmem=3G -pe threaded $threads $hold -t 1-$numfiles:1 $script_path/align_bwa.sh $sample $output_dir $run_info
					hold="$type.$version.align_bwa.$sample.$run_num"
				fi	    
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.processBAM.$sample.$run_num -pe threaded $threads -l h_vmem=4G -hold_jid $hold $script_path/processBAM.sh $align_dir $sample $run_info   
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.extract_reads_bam.$sample.$run_num -l h_vmem=8G -hold_jid $type.$version.processBAM.$sample.$run_num $script_path/extract_reads_bam.sh $align_dir $bamfile $run_info $output_dir/IGV_BAM		
			elif [[ $analysis == "realignment" || $analysis == "realign-mayo" ]]
			then
				infile=`cat $sample_info | grep -w ^BAM:${sample} | cut -d '=' -f2 `
				num_bams=`echo $infile | tr " " "\n" | wc -l`
				for ((i=1; i <=$num_bams; i++));
				do
					bam=`echo $infile | awk -v num=$i '{print $num}'`
					$samtools/samtools view -H $input/$bam 1>$input/$bam.header 2> $align_dir/$sample.$i.sorted.bam.log
					if [ `cat $align_dir/$sample.$i.sorted.bam.log | wc -l` -gt 0 ]
					then
						echo "$input/$bam : bam file is truncated or corrupted" 	
						exit 1;
					else
						rm $align_dir/$sample.$i.sorted.bam.log
					fi
					rm $input/$bam.header
					ln -s $input/$bam $align_dir/$sample.$i.sorted.bam
				done
				$script_path/check_qstat.sh $limit
				qsub $args -pe threaded $threads -N $type.$version.processBAM.$sample.$run_num -l h_vmem=4G $script_path/processBAM.sh $align_dir $sample $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.extract_reads_bam.$sample.$run_num -l h_vmem=8G -hold_jid $type.$version.processBAM.$sample.$run_num $script_path/extract_reads_bam.sh $align_dir $bamfile $run_info $output_dir/IGV_BAM
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
					infile=`cat $sample_info | grep -w ^BAM:${group} | cut -d '=' -f2`
					num_bams=`echo $infile | tr " " "\n" | wc -l`
					for ((i=1; i <=$num_bams; i++));
					do
						bam=`echo $infile | awk -v num=$i '{print $num}'`
						$samtools/samtools view -H $input/$bam 1>$input/$bam.header 2> $realign_dir/$group.$i.sorted.bam.log
						if [ `cat $realign_dir/$group.$i.sorted.bam.log | wc -l` -gt 0 ]
						then
							echo "$input/$bam : bam file is truncated or corrupted" 	
							exit 1;
						else
							rm $realign_dir/$group.$i.sorted.bam.log
						fi	
						rm $input/$bam.header
						ln -s $input/$bam $realign_dir/$group.$i.sorted.bam
					done
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.reformat_pairBAM.$group.$run_num -l h_vmem=8G $script_path/reformat_pairBAM.sh $realign_dir $group $run_info
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.extract_reads_bam.$group.$run_num -l h_vmem=8G -hold_jid $type.$version.reformat_pairBAM.$group.$run_num $script_path/extract_reads_bam.sh $realign_dir $group.sorted.bam $run_info $output_dir/IGV_BAM $group
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.split_bam_chr.$group.$run_num -hold_jid $type.$version.reformat_pairBAM.$group.$run_num -l h_vmem=2G -t 1-$numchrs:1 $script_path/split_bam_chr.sh $realign_dir $group $run_info
					variant_id="$type.$version.split_bam_chr.$group.$run_num"
				else        
					id=""
					for sample in $samples
					do
						id=$id"$type.$version.processBAM.$sample.$run_num,$type.$version.extract_reads_bam.$sample.$run_num"
					done    
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.realign_recal.$group.$run_num -hold_jid $id -l h_vmem=8G -t 1-$numchrs:1 $script_path/realign_recal.sh $input_dirs $bam_samples $names_samples $realign_dir $run_info 1
					variant_id="$type.$version.realign_recal.$group.$run_num"
				fi
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.split_sample_pair.$group.$run_num -hold_jid $variant_id -l h_vmem=2G -t 1-$numchrs:1 $script_path/split_sample_pair.sh $output_dir/realign $output_dir/IGV_BAM $group $output_dir/alignment $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.igv_bam.$group.$run_num -hold_jid $type.$version.split_sample_pair.$group.$run_num -l h_vmem=2G $script_path/igv_bam.sh $output_dir/realign $output_dir/IGV_BAM $group $output_dir/alignment $run_info 
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.variants.$group.$run_num -hold_jid $variant_id -pe threaded $threads -l h_vmem=4G -t 1-$numchrs:1 $script_path/variants.sh $realign_dir $names_samples $variant_dir 1 $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.merge_variant_group.$group.$run_num -l h_vmem=3G -pe threaded $threads -hold_jid $type.$version.variants.$group.$run_num $script_path/merge_variant_group.sh $output_dir/variants $group $output_dir/Reports_per_Sample/ $run_info 
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.OnTarget_BAM.$group.$run_num -hold_jid $type.$version.split_sample_pair.$group.$run_num -l h_vmem=2G -t 1-$numchrs:1 $script_path/OnTarget_BAM.sh $output_dir/IGV_BAM $output_dir/OnTarget $group $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.OnTarget_PILEUP.$group.$run_num -hold_jid $type.$version.split_sample_pair.$group.$run_num -l h_vmem=6G -t 1-$numchrs:1 $script_path/OnTarget_PILEUP.sh $realign_dir $output_dir/OnTarget $group $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.getCoverage.$group.$run_num -hold_jid $type.$version.OnTarget_PILEUP.$group.$run_num -l h_vmem=2G $script_path/getCoverage.sh $output_dir/OnTarget $output_dir/numbers $group $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.OnTarget_variant.$group.$run_num -l h_vmem=2G -hold_jid $type.$version.merge_variant_group.$group.$run_num -t 1-$numchrs:1 $script_path/OnTarget_variant.sh $output_dir/variants $output_dir/OnTarget $group $run_info
				hold_args="-hold_jid $type.$version.OnTarget_variant.$group.$run_num"
				## SIFT
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.sift.$group.$run_num $hold_args -t 1-$numchrs:1 -l h_vmem=4G $script_path/sift.sh $sift $output_OnTarget $group $run_info 
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.sift.$group.$run_num $hold_args -t 1-$numchrs:1 -l h_vmem=4G $script_path/sift.sh $sift $output_OnTarget $group $run_info TUMOR
				## SNPEFF
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.snpeff.$group.$run_num $hold_args -t 1-$numchrs:1 -l h_vmem=4G $script_path/snpeff.sh $snpeff $output_OnTarget $group $run_info
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.snpeff.$group.$run_num $hold_args -t 1-$numchrs:1 -l h_vmem=4G $script_path/snpeff.sh $snpeff $output_OnTarget $group $run_info TUMOR
				##POLYPHEN
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.polyphen.$group.$run_num $hold_args -t 1-$numchrs:1 -pe threaded $threads -l h_vmem=4G $script_path/polyphen.sh $polyphen $output_OnTarget $group $run_info  
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.polyphen.$group.$run_num $hold_args -t 1-$numchrs:1 -pe threaded $threads -l h_vmem=4G $script_path/polyphen.sh $polyphen $output_OnTarget $group $run_info TUMOR
				hold="$type.$version.sift.$group.$run_num,$type.$version.snpeff.$group.$run_num,$type.$version.polyphen.$group.$run_num"
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.sample_reports.$group.$run_num -hold_jid $hold -t 1-$numchrs:1 -l h_vmem=4G $script_path/sample_reports.sh $run_info $group $TempReports $output_OnTarget $sift $snpeff $polyphen $output_dir
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.sample_report.$sample.$run_num -l h_vmem=2G -hold_jid $type.$version.sample_reports.$group.$run_num $script_path/sample_report.sh $output_dir $TempReports $group $run_info
				hold="$type.$version.sift.${group}.$run_num,$type.$version.snpeff.$group.$run_num,$type.$version.polyphen.$group.$run_num"
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.sample_reports.$group.$run_num -hold_jid $hold -t 1-$numchrs:1 -l h_vmem=4G $script_path/sample_reports.sh $run_info $group $TempReports $output_OnTarget $sift $snpeff $polyphen $output_dir TUMOR
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.sample_report.$group.$run_num -l h_vmem=2G -hold_jid $type.$version.sample_reports.$group.$run_num $script_path/sample_report.sh $output_dir $TempReports $group $run_info TUMOR
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
						qsub $args -N $type.$version.run_crest_multi_cover.$group.$sam.$run_num -hold_jid $type.$version.split_sample_pair.$group.$run_num -l h_vmem=8G -t 1-$numchrs:1 $script_path/run_crest_multi_cover.sh $sam $group $output_dir/IGV_BAM $crest $run_info
						id=$id"$type.$version.run_crest_multi_cover.$group.$sam.$run_num,"
					done
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.run_crest_multi.$group.$run_num -hold_jid $id -l h_vmem=6G -t 1-$numchrs:1 $script_path/run_crest_multi.sh $group $output_dir/IGV_BAM $crest $run_info
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.run_segseq.$group.$run_num -hold_jid $type.$version.split_sample_pair.$group.$run_num -l h_vmem=16G -t 1-$numchrs:1 -l matlab_lic=1 $script_path/run_segseq.sh $group $output_dir/IGV_BAM $output_dir/cnv $run_info    
					let nump=$numchrs+1;    
					mkdir -p $break/$group
					id=""
					for sam in `cat $sample_info| grep -w "^$group" | cut -d '=' -f2`
					do
						$script_path/check_qstat.sh $limit
						qsub $args -N $type.$version.run_breakdancer.$group.$sam.$run_num -hold_jid $type.$version.split_sample_pair.$group.$run_num -l h_vmem=4G -t 1-$numchrs:1 $script_path/run_breakdancer.sh $sam $output_dir/IGV_BAM $break/$group $run_info $group
						$script_path/check_qstat.sh $limit
						qsub $args -N $type.$version.run_breakdancer_in.$group.$sam.$run_num -hold_jid $type.$version.igv_bam.$group.$run_num -l h_vmem=3G -t $nump-$nump:$nump $script_path/run_breakdancer.sh $sam $output_dir/IGV_BAM $break/$group $run_info $group
						id=$id"$type.$version.run_breakdancer.$group.$sam.$run_num,$type.$version.run_breakdancer_in.$group.$sam.$run_num,"
					done
					hhold="$id,$type.$version.run_segseq.$group.$run_num,$type.$version.run_crest_multi.$group.$run_num"
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.summaryze_struct_group.$group.$run_num -l h_vmem=5G -hold_jid $hhold $script_path/summaryze_struct_group.sh $group $output_dir $run_info
					mkdir -p $output_dir/circos;
					for i in $(seq 2 ${#sampleArray[@]})
					do  
						tumor=${sampleArray[$i]}
						$script_path/check_qstat.sh $limit
						qsub $args -N $type.$version.plot_circos_cnv_sv.$group.$tumor.$i.$run_num -hold_jid $type.$version.summaryze_struct_group.$group.$run_num -l h_vmem=2G $script_path/plot_circos_cnv_sv.sh $output_dir/struct/$group.$tumor.somatic.break $output_dir/struct/$group.$tumor.somatic.filter.crest $output_dir/cnv/$group/$tumor.cnv.filter.bed $group.$tumor $output_dir/circos $run_info
					done
				fi
			done
			mkdir -p $output_dir/Reports_per_Sample/ANNOT
			if [ $tool == "whole_genome" ]
			then
				id=""
				for group in `echo $groups | tr ":" "\n"`
				do
					samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" "\n")
					sampleArray=()
					i=1
					for sample in $sampleNames
					do
					sampleArray[$i]=$sample
					let i=i+1
					done
					for i in $(seq 2 ${#sampleArray[@]})
					do
						tumor=${sampleArray[$i]}
						id=$id"$type.$version.plot_circos_cnv_sv.$group.$tumor.$i.$run_num,"
					done
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.annotation_CNV.$group.$run_num -l h_vmem=2G -hold_jid $id $script_path/annotation_CNV.sh $output_dir/Reports_per_Sample/SV/ $run_info $output_dir/Reports_per_Sample/ANNOT $group
					$script_path/check_qstat.sh $limit
					qsub $args -N $type.$version.annotation_SV.$group.$run_num -l h_vmem=2G -hold_jid $id $script_path/annotation_SV.sh $output_dir $run_info $output_dir/Reports_per_Sample/ANNOT/ $group
				done
			fi
			### generate reports for all the samples
			id=""
			for group in `echo $groups | tr ":" "\n"`
			do
				id=$id"$type.$version.sample_report.$group.$run_num,"	
			done
			for group in `echo $groups | tr ":" "\n"`
			do
				samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" " ")
				for sample in $samples
				do
					id=$id"$type.$version.extract_reads_bam.$sample.$run_num,"
				done
			done

			if [ $tool == "whole_genome" ]
			then
				for group in `echo $groups | tr ":" "\n"`
				do
					id=$id"$type.$version.summaryze_struct_group.$group.$run_num,$type.$version.annotation_CNV.$group.$run_num,$type.$version.annotation_SV.$group.$run_num,"
				done
			fi    
			for group in `echo $groups | tr ":" "\n"`
			do
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.sample_numbers.$group.$run_num -l h_vmem=2G -hold_jid $id $script_path/sample_numbers.sh $output_dir $group $run_info $output_dir/numbers
				$script_path/check_qstat.sh $limit
				qsub $args -N $type.$version.gene_summary.$group.$run_num -l h_vmem=2G -hold_jid $id $script_path/gene_summary.sh $output_dir $group $run_info $output_dir/Reports_per_Sample
			done

			if [ $tool == "exome" ]
			then
				id=""
				for group in `echo $groups | tr ":" "\n"`
				do
					id=$id"$type.$version.getCoverage.$group.$run_num,$type.$version.sample_numbers.$group.$run_num,$type.$version.gene_summary.$group.$run_num,$type.$version.igv_bam.$group.$run_num,$type.$version.annotate_sample.$run_num,$type.$version.annotation_CNV.$run_num,$type.$version.annotation_SV.$run_num,"
				done          
			elif [ $tool == "whole_genome" ]
			then
				id=""
				for group in `echo $groups | tr ":" "\n"`
				do
					id=$id"$type.$version.getCoverage.$group.$run_num,$type.$version.sample_numbers.$group.$run_num,$type.$version.gene_summary.$group.$run_num,$type.$version.igv_bam.$group.$run_num,$type.$version.annotate_sample.$run_num,"
				done    
			fi
			if [[ $upload_tb == "YES" ]]
			then
				mem="8G"
			else
				mem="2G"
			fi	
			$script_path/check_qstat.sh $limit
			qsub $args -N $type.$version.generate_html.$run_num -l h_vmem=$mem -hold_jid $id $script_path/generate_html.sh $output_dir $run_info 	
		fi
	fi
	echo `date`
fi
