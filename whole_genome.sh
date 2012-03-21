#!/bin/sh

########################################################
###### 	MASTER SCRIPT FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			whole_genome_pipeline.sh
######		Date:				07/25/2011
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

	input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
	output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	dos2unix $sample_info
	dos2unix $tool_info
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
	numchrs=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
	paired=$( cat $run_info | grep -w '^PAIRED' | cut -d '=' -f2)
	threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2)
	variant_type=$(cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")   
    bed=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
    PATH=$bed/:$PATH
	perl $script_path/check_config.pl $run_info

	if [ -d $output/$PI ]
	then
		echo "$PI folder exists"
	else
		mkdir $output/$PI
	fi

	if [ -d $output/$PI/$tool ]
	then 	
		echo "$tool analysis folder exists"
	else
		mkdir $output/$PI/$tool
	fi

	if [ -d $output/$PI/$tool/$run_num ]
	then
		echo "ERROR : $run_num folder exists"
		exit 1;
	else 
		mkdir $output/$PI/$tool/$run_num
	fi

	mkdir $output/$PI/$tool/$run_num/logs
	output_dir=$output/$PI/$tool/$run_num
	$script_path/create_folder.sh $output_dir $run_info
	
	job_ids_dir=$output_dir/job_ids
	output_align=$output_dir/alignment
	if [ $analysis != "alignment" ]
    then
        output_OnTarget=$output_dir/OnTarget
        output_annot=$output_dir/annotation
        TempReports=$output_dir/TempReports
        sift=$output_annot/SIFT
        sseq=$output_annot/SSEQ
	fi

	touch $output_dir/tool_info.txt
	cat $tool_info > $output_dir/tool_info.txt
	touch $output_dir/sample_info.txt
	cat $sample_info > $output_dir/sample_info.txt
	touch $output_dir/run_info.txt
	cat $run_info > $output_dir/run_info.txt
	cp $script_path/${tool}_workflow.png $output_dir/${tool}_workflow.png
	cp $script_path/IGV_Setup.doc $output_dir/IGV_Setup.doc
	cp $script_path/ColumnDescription_Reports.xls $output_dir/
	
    if [ $tool == "whole_genome" ]
    then
        if [ $analysis != "alignment" ]
        then
            cat $master_gene_file | sort -n -k 1,12n -k 2,12n > $output_dir/bed_file.bed.temp
            $bed/mergeBed -i $output_dir/bed_file.bed.temp | awk '$1 !~ /random/ && $1 !~ /hap/ && $1 !~ /chrUn/' >  $output_dir/bed_file.bed
            rm $output_dir/bed_file.bed.temp	
        fi
    fi

    echo -e "${tool} analysis for ${run_num} for ${PI} " >> $output_dir/log.txt
	START=`date`
	echo -e "Analysis started at:" >> $output_dir/log.txt
	echo -e "${START}" >>  $output_dir/log.txt

	if [[ $analysis != "mayo" && $analysis != "external"  && $analysis != "realignment"  &&  $analysis != "variant" && $analysis != "alignment" && $analysis != "annotation" ]]
	then
		echo -e "\nPlease Specify the correct Analysis type(alignment,realignment,variant,external,mayo,annotation)\n"
		echo `date`
		exit 1;
	fi

	args="-V -wd $output_dir/logs -q $queue -m a -M $email -l h_stack=10M"
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
				mkdir -p $align_dir
			fi
			if [ $analysis == "mayo" -o $analysis == "external" -o $analysis == "alignment" ]
			then
				if [ $paired == 1 ]
				then
					let numfiles=(`cat $sample_info | grep -w "$sample" | cut -d '=' -f2| tr "\t" "\n" |wc -l`)/2
				else
					let numfiles=`cat $sample_info | grep -w "$sample" | cut -d '=' -f2| tr "\t" "\n" |wc -l`
				fi	
				
				if [ $aligner == "novoalign" ]
				then
					echo "novoalign is used as aligner"
					ALIGNMENT=`qsub $args -N $type.$version.align_novo.$sample.$run_num -l h_vmem=8G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_novo.sh $sample $output_dir $run_info`
				elif [ $aligner == "bwa" ]
				then
					echo "bwa is used as aligner"
					READ1=`qsub $args -N $type.$version.align_read_bwa.R1.$sample.$run_num -l h_vmem=8G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_read_bwa.sh $sample $output_dir 1 $run_info`
					r1_job=$(echo $READ1 | cut -d ' ' -f3 |tr "\n" "," | sed -e "s/\..*,//g")
					if [ $paired == 1 ]
					then
						READ2=`qsub $args -N $type.$version.align_read_bwa.R2.$sample.$run_num -l h_vmem=8G -pe threaded $threads -t 1-$numfiles:1 $script_path/align_read_bwa.sh $sample $output_dir 2 $run_info`
						r2_job=$(echo $READ2 | cut -d ' ' -f3 |tr "\n" "," | sed -e "s/\..*,//g")	
						hold="-hold_jid $r1_job,$r2_job"
					else
						hold="-hold_jid $r1_job"
					fi	
					ALIGNMENT=`qsub $args -N $type.$version.align_bwa.$sample.$run_num -l h_vmem=8G -pe threaded $threads $hold -t 1-$numfiles:1 $script_path/align_bwa.sh $sample $output_dir $run_info`
				else
					echo "Doesn't support the aligner"
				fi	
				job_id_align=`echo $ALIGNMENT | cut -d ' ' -f3| tr "\n" "," | sed -e "s/\..*,//g"`
				MERGE=`qsub $args -N $type.$version.processBAM.$sample.$run_num -l h_vmem=8G -hold_jid $job_id_align $script_path/processBAM.sh $align_dir $sample $run_info`
				job_id_convert=`echo $MERGE | cut -d ' ' -f3 `
				echo -e "$MERGE" >> $job_ids_dir/ALIGN	
			elif [ $analysis == "realignment" ]
			then
				infile=`cat $sample_info | grep -w "^$sample" | cut -d '=' -f2`
				num_bams=`echo $infile | tr " " "\n" | wc -l`
				for ((i=1; i <=$num_bams; i++));
				do
					bam=`echo $infile | awk -v num=$i '{print $num}'`
					ln -s $input/$bam $align_dir/$sample.$i.sorted.bam
				done  
				CONVERT=`qsub $args -N $type.$version.processBAM.$sample.$run_num -l h_vmem=8G $script_path/processBAM.sh $align_dir $sample $run_info`
				job_id_convert=`echo $CONVERT | cut -d ' ' -f3`
			fi    

			if [[ $analysis == "mayo" || $analysis == "external" || $analysis == "realignment" || $analysis == "variant" ]]
			then
				realign_dir=$output_dir/realign/$sample
				variant_dir=$output_dir/variants/$sample
				mkdir -p $realign_dir $variant_dir

				if [ $analysis == "variant" ]
				then
					infile=`cat $sample_info | grep -w "^$sample" | cut -d '=' -f2`
					num_bams=`echo $infile | tr " " "\n" | wc -l`
					for ((i=1; i <=$num_bams; i++));
					do
						bam=`echo $infile | awk -v num=$i '{print $num}'`
						ln -s $input/$bam $realign_dir/$sample.$i.sorted.bam
					done
					REFORMAT=`qsub $args -N $type.$version.reformat_BAM.$sample.$run_num -l h_vmem=8G $script_path/reformat_BAM.sh $realign_dir $sample $run_info`	
					job_id_convert=`echo $REFORMAT | cut -d ' ' -f3`
					VARIANT=`qsub $args -N $type.$version.split_bam_chr.$sample.$run_num -hold_jid $job_id_convert -l h_vmem=8G -t 1-$numchrs:1 $script_path/split_bam_chr.sh $realign_dir $sample $run_info`					
				else
					VARIANT=`qsub $args -N $type.$version.realign_recal.$sample.$run_num -hold_jid $job_id_convert -l h_vmem=8G -t 1-$numchrs:1 $script_path/realign_recal.sh $align_dir $bamfile $sample $realign_dir $variant_dir $run_info 1`	
				fi
				variant_id=`echo $VARIANT | cut -d ' ' -f3  | tr "\n" "," | sed -e "s/\..*,//g"`
				IGV=`qsub $args -N $type.$version.igv_bam.$sample.$run_num -l h_vmem=8G -hold_jid $variant_id $script_path/igv_bam.sh $output_dir/realign $output_dir/IGV_BAM $sample $output_dir/alignment $run_info`
                igv_id=`echo $IGV | cut -d ' ' -f3 ` 
				echo -e $IGV >> $job_ids_dir/IGV	
                CALLS=`qsub $args -N $type.$version.variants.$sample.$run_num -hold_jid $variant_id -pe threaded $threads -l h_vmem=8G -t 1-$numchrs:1 $script_path/variants.sh $realign_dir $sample $variant_dir 1 $run_info`
				calls_id=`echo $CALLS | cut -d ' ' -f3  | tr "\n" "," | sed -e "s/\..*,//g"`
				MERGEVARIANT=`qsub $args -N $type.$version.merge_variant_single.$sample.$run_num -l h_vmem=8G -hold_jid $calls_id $script_path/merge_variant_single.sh $output_dir/variants $sample $output_dir/Reports_per_Sample/ $run_info`
				job_id_mergevariant=`echo $MERGEVARIANT | cut -d ' ' -f3`
				ONTARGET_B=`qsub $args -N $type.$version.OnTarget_BAM.$sample.$run_num -hold_jid $variant_id -l h_vmem=8G -t 1-$numchrs:1 $script_path/OnTarget_BAM.sh $realign_dir $output_dir/OnTarget $sample $run_info`
				ONTARGET_P=`qsub $args -N $type.$version.OnTarget_PILEUP.$sample.$run_num -hold_jid $variant_id -l h_vmem=8G -t 1-$numchrs:1 $script_path/OnTarget_PILEUP.sh $realign_dir $output_dir/OnTarget $sample $run_info`
				job_ids_target_p=`echo $ONTARGET_P | cut -d ' ' -f3 |  tr "\n" "," | sed -e "s/\..*,//g"` 
				job_ids_target_b=`echo $ONTARGET_B | cut -d ' ' -f3 |  tr "\n" "," | sed -e "s/\..*,//g"` 
				COVERAGE=`qsub $args -N $type.$version.getCoverage.$sample.$run_num -hold_jid ${job_ids_target_p},${job_ids_target_b} $script_path/getCoverage.sh $output_dir/OnTarget $output_dir/numbers $sample $run_info`
				job_ids_coverage=`echo $COVERAGE | cut -d ' ' -f3| tr "\n" ","`
				ANNOTATE_SNV=`qsub $args -N $type.$version.OnTarget_variant.$sample.$run_num -l h_vmem=4G -hold_jid $job_id_mergevariant -t 1-$numchrs:1 $script_path/OnTarget_variant.sh $output_dir/variants $output_dir/OnTarget $sample $run_info`
				job_id_annotate_snv=`echo $ANNOTATE_SNV | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
			fi	
			if [ $analysis == "annotation" ]
			then
				if [ $variant_type == "BOTH" ]
				then
					REFORMAT=`qsub $args -N $type.$version.reformat_VARIANTs.$sample.$run_num -l h_vmem=4G $script_path/reformat_VARIANTs.sh $output_OnTarget $sample $run_info 2`
				elif [ $variant_type == "SNV" -o $variant_type == "INDEL" ]
				then
					REFORMAT=`qsub $arg -N $type.$version.reformat_VARIANTs.$sample.$run_num -t 1-$numsamples:1 -l h_vmem=4G $script_path/reformat_VARIANTs.sh $output_OnTarget $sample_info $run_info 1`
				fi
				job_ids_format=`echo $REFORMAT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
				hold_args="-hold_jid $job_ids_format"
			elif [ $analysis != "alignment" ]
			then
				hold_args="-hold_jid $job_id_annotate_snv"
			fi
			if [ $analysis != "alignment" ]
			then
				if [ $variant_type == "SNV" -o $variant_type == "BOTH" ]
				then
					SIFT=`qsub $args -N $type.$version.sift.$sample.$run_num $hold_args -t 1-$numchrs:1 -l h_vmem=4G $script_path/sift.sh $sift $output_OnTarget $sample $run_info` 
				fi
				SSEQ=`qsub $args -N $type.$version.sseq.$sample.$run_num $hold_args -l h_vmem=4G $script_path/sseq.sh $sseq $output_OnTarget $email $sample $run_info`
				echo -e $SSEQ >> $job_ids_dir/ANNOT	
				job_ids_sift=`echo $SIFT | cut -d ' ' -f3|  tr "\n" "," | sed -e "s/\..*,//g"`
				job_ids_sseq=`echo $SSEQ | cut -d ' ' -f3`
				ADD_ANOT=`qsub $args -N $type.$version.sample_reports.$sample.$run_num -hold_jid $job_ids_sseq,$job_ids_sift -t 1-$numchrs:1 -l h_vmem=8G $script_path/sample_reports.sh $run_info $sample $TempReports $output_OnTarget $sift $sseq $output_dir`
				job_ids_anot=`echo $ADD_ANOT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
				MERGE=`qsub $args -N $type.$version.sample_report.$sample.$run_num -hold_jid $job_ids_anot $script_path/sample_report.sh $output_dir $TempReports $sample $run_info`
				echo -e $MERGE >> $job_ids_dir/PER_SAMPLE

				if [[ $tool == "whole_genome"  && $analysis != "annotation" && $analysis != "alignment" ]]
				then
					crest=$output_dir/struct/crest
					break=$output_dir/struct/break
					cnv=$output_dir/cnv/$sample
					mkdir -p $break
					mkdir -p $crest
					mkdir -p $cnv
					### not smooth for all the samples so commenting out for time being
					SVCALL=`qsub $args -N $type.$version.run_single_crest.sh.$sample.$run_num -hold_jid $variant_id -t 1-$numchrs:1 -l h_vmem=8G $script_path/run_single_crest.sh $sample $realign_dir $crest $run_info`
					job_id_sv=`echo $SVCALL | cut -d ' ' -f3 |  tr "\n" "," | sed -e "s/\..*,//g"`
					CNVCALL=`qsub $args -N $type.$version.run_cnvnator.$sample.$run_num -hold_jid $variant_id -l h_vmem=8G -t 1-$numchrs:1 $script_path/run_cnvnator.sh $sample $realign_dir $cnv $run_info`
					job_id_cnv=`echo $CNVCALL | cut -d ' ' -f3 |  tr "\n" "," | sed -e "s/\..*,//g"`
					let nump=$numchrs+1;
					BREAKDANCER=`qsub $args -N $type.$version.run_breakdancer.$sample.$run_num -hold_jid $variant_id -l h_vmem=8G -t 1-$numchrs:1 $script_path/run_breakdancer.sh $sample $output_dir/realign $break $run_info`
					job_id_break=`echo $BREAKDANCER | cut -d ' ' -f3 |  tr "\n" "," | sed -e "s/\..*,//g"`
					BREAKDANCER_IN=`qsub $args -N $type.$version.run_breakdancer.$sample.$run_num -hold_jid $igv_id -l h_vmem=8G -t $nump-$nump:$nump $script_path/run_breakdancer.sh $sample $output_dir/IGV_BAM $break $run_info`
					job_id_break_in=`echo $BREAKDANCER_IN | cut -d ' ' -f3 |  tr "\n" "," | sed -e "s/\..*,//g"`
					### merge the structural variants
					MERGESTRUCT=`qsub $args -N $type.$version.summaryze_struct_single.$sample.$run_num -l h_vmem=8G -hold_jid $job_id_break_in,$job_ids_coverage,$job_id_sv,$job_id_cnv,$job_id_break $script_path/summaryze_struct_single.sh $sample $output_dir $run_info`
					echo -e $MERGESTRUCT >> $job_ids_dir/SV
					job_id_mergestruct=`echo $MERGESTRUCT | cut -d ' ' -f3`
					RUNCIRCOS=`qsub $args -N $type.$version.plot_circos_cnv_sv.$sample.$run_num -hold_jid $job_id_mergestruct -l h_vmem=8G $script_path/plot_circos_cnv_sv.sh $break/$sample/$sample.break $crest/$sample/$sample.filter.crest $cnv/$sample.cnv.filter.bed $sample $output_dir/circos $run_info` 
				fi
			fi	
		done
        if [ $analysis != "alignment" ]
		then
			job_ids_annn=$( cat $job_ids_dir/ANNOT |  cut -d ' ' -f3 | tr "\n" ",")
		fi
		if [[ $analysis != "annotation" && $analysis != "alignment" ]]
		then
			if [ $tool == "whole_genome" ]
			then
				job_ids=$( cat $job_ids_dir/SV | cut -d ' ' -f3  | tr "\n" "," )
				mkdir -p $output_dir/Reports_per_Sample/ANNOT
				ANNOTATE_CNV=`qsub $args -N $type.$version.annotation.CNV.sh.$run_num -l h_vmem=4G -hold_jid $job_ids -t 1-$numsamples:1 $script_path/annotation_CNV.sh $output_dir/Reports_per_Sample/SV/ $run_info $output_dir/Reports_per_Sample/ANNOT`
				ANNOTATE_SV=`qsub $args -N $type.$version.annotation.SV.sh.$run_num -l h_vmem=4G -hold_jid $job_ids -t 1-$numsamples:1 $script_path/annotation_SV.sh $output_dir $run_info $output_dir/Reports_per_Sample/ANNOT`
				job_id_annotate_sv=`echo $ANNOTATE_SV | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
				job_id_annotate_cnv=`echo $ANNOTATE_CNV | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
				hold_arg="-hold_jid $job_id_annotate_snv,$job_id_annotate_cnv,$job_id_annotate_sv,$job_ids_annn"    
			else
				hold_arg="-hold_jid $job_ids_annn"	
			fi	
		else
			hold_arg="-hold_jid $job_ids_annn"
		fi
		if [ $analysis != "alignment" ]
		then
			REPORT=`qsub $args -N $type.$version.merged_report.$run_num $hold_arg -t 1-$numchrs:1 -l h_vmem=8G $script_path/merged_report.sh $sift $sseq $TempReports $run_info $output_dir/OnTarget` 
			job_ids_report=`echo $REPORT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
			job_ids_sample=$( cat $job_ids_dir/PER_SAMPLE | cut -d ' ' -f3  | tr "\n" ",")
			ANNOT_SAMPLE=`qsub $args -N $type.$version.annotate_sample.$run_num -hold_jid $job_ids_sample -l h_vmem=8G $script_path/annotate_sample.sh $output_dir $run_info`   
			job_ids_annot_sample=`echo $ANNOT_SAMPLE | cut -d ' ' -f3 | tr "\n" ","`
			## annotate merged report
			MERGE_SAMPLE=`qsub $args -N $type.$version.annotate_merged_report.$run_num -hold_jid $job_ids_report -l h_vmem=8G $script_path/annotate_merged_report.sh $output_dir $TempReports $run_info`
			job_ids_merge_sample=`echo $MERGE_SAMPLE | cut -d ' ' -f3 | tr "\n" ","`
			if [ $analysis != "annotation" ]
			then
				job_ids_igv=$( cat $job_ids_dir/IGV | cut -d ' ' -f3  | tr "\n" ",")
				hold="-hold_jid $job_ids_report,$job_ids_annot_sample,$job_ids_merge_sample,$job_ids_sample,$job_ids_annn,$job_ids_igv"
			else
				hold="-hold_jid $job_ids_report,$job_ids_annot_sample,$job_ids_merge_sample,$job_ids_sample,$job_ids_annn"
			fi	
			### variant distance
			qsub $args -N $type.$version.variant_distance.$run_num -hold_jid $job_ids_report -l h_vmem=4G $script_path/variant_distance.sh $TempReports $output_dir $run_info
			NUMBERS=`qsub $args -N $type.$version.sample_numbers.$run_num $hold -t 1-$numsamples:1 $script_path/sample_numbers.sh $output_dir $run_info`
			GENE_SUMMARY=`qsub $args -N $type.$version.gene_summary.$run_num $hold -t 1-$numsamples:1 $script_path/gene_summary.sh $output_dir $run_info $output_dir/Reports_per_Sample`
			job_ids=`echo $NUMBERS | cut -d ' ' -f3 | cut -d '.' -f1 | tr "\n" ","`
			job_ids_summary=`echo $GENE_SUMMARY | cut -d ' ' -f3 | cut -d '.' -f1 | tr "\n" ","`
			HTML=`qsub $args -N $type.$version.generate_html.$run_num -hold_jid ${job_ids}${job_ids_summary} $script_path/generate_html.sh $output_dir $run_info`
			if [[  $tool == "exome"  && $all_sites == "YES" ]]
			then
				qsub $args -N $type.$version.merge_raw_variants.$run_num -hold_jid ${job_ids}${job_ids_summary} $script_path/merge_raw_variants.sh $output_dir $run_info
			fi	
			
		else
			job_ids=$( cat $job_ids_dir/ALIGN |  cut -d ' ' -f3 | tr "\n" ",")
			NUMBERS=`qsub $args -N $type.$version.sample_numbers.$run_num -hold_jid $job_ids -t 1-$numsamples:1 $script_path/sample_numbers.sh $output_dir $run_info`
            job_ids=`echo $NUMBERS | cut -d ' ' -f3 | cut -d '.' -f1 | tr "\n" ","`
            HTML=`qsub $args -N $type.$version.generate_html.$run_num -hold_jid ${job_ids} $script_path/generate_html.sh $output_dir $run_info`
		fi	
		echo `date`        
	else
		echo "Multi-sample"
		numgroups=$(cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
		for group in `echo $groups | tr ":" "\n"`
		do
			samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" "\n")
			job_ids_merge=""
			#bam_files=""
			bam_samples=""
			input_dirs=""
            names_samples=""
			for sample in $samples
			do
				align_dir=$output_dir/alignment/$sample;
				mkdir -p $align_dir

				if [[ $analysis == "mayo" || $analysis == "external" || $analysis == "alignment" ]]
				then
					if [ $paired == 1 ]
					then
						let numfiles=(`cat $sample_info | grep -w "$sample" | cut -d '=' -f2| tr "\t" "\n" |wc -l`)/2
					else
						let numfiles=`cat $sample_info | grep -w "$sample" | cut -d '=' -f2| tr "\t" "\n" |wc -l`
					fi
					
					if [ $aligner == "novoalign" ]
					then
						echo "novoalign is used as aligner"
						ALIGNMENT=`qsub $args -N $type.$version.align_novo.$sample.$run_num -l h_vmem=8G -pe threaded 4 -t 1-$numfiles:1 $script_path/align_novo.sh $sample $output_dir $run_info`
					elif [ $aligner == "bwa" ]
					then
						echo "bwa is used as aligner"
						READ1=`qsub $args -N $type.$version.align_read_bwa.R1.$sample.$run_num -l h_vmem=8G -pe threaded 4 -t 1-$numfiles:1 $script_path/align_read_bwa.sh $sample $output_dir 1 $run_info`
						r1_job=$(echo $READ1 | cut -d ' ' -f3 |tr "\n" "," | sed -e "s/\..*,//g")
						if [ $paired == 1 ]
						then
							READ2=`qsub $args -N $type.$version.align_read_bwa.R2.$sample.$run_num -l h_vmem=8G -pe threaded 4 -t 1-$numfiles:1 $script_path/align_read_bwa.sh $sample $output_dir 2 $run_info`
							r2_job=$(echo $READ2 | cut -d ' ' -f3 |tr "\n" "," | sed -e "s/\..*,//g")	
							hold="-hold_jid $r1_job,$r2_job"
						else
							hold="-hold_jid $r1_job"
						fi	
						ALIGNMENT=`qsub $args -N $type.$version.align_bwa.$sample.$run_num -l h_vmem=8G -pe threaded 4 $hold -t 1-$numfiles:1 $script_path/align_bwa.sh $sample $output_dir $run_info`
					else
						echo "ERROR : Doesn't support the aligner"
						exit 1
					fi	    
					job_id_align=`echo $ALIGNMENT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
					MERGE=`qsub $args -N $type.$version.processBAM.$sample.$run_num -l h_vmem=8G -hold_jid $job_id_align $script_path/processBAM.sh $align_dir $sample $run_info`
					job_id_convert=`echo $MERGE | cut -d ' ' -f3 `
					job_id_convert="$job_id_convert,$job_id_convert"                
				elif [[ $analysis == "realignment" ]]
				then
					infile=`cat $sample_info | grep -w "^$sample" | cut -d '=' -f2 `
					num_bams=`echo $infile | tr " " "\n" | wc -l`
					for ((i=1; i <=$num_bams; i++));
					do
						bam=`echo $infile | awk -v num=$i '{print $num}'`
						ln -s $input/$bam $align_dir/$sample.$i.sorted.bam
					done
					CONVERT=`qsub $args -N $type.$version.processBAM.$sample.$run_num -l h_vmem=8G $script_path/processBAM.sh $align_dir $sample $run_info`
					job_id_convert=`echo $CONVERT | cut -d ' ' -f3`
					job_ids_convert="$job_ids_convert,$job_id_convert"
				fi
				### setting the bams and its path for multiple sample anaylysis
				names_samples=$names_samples"$sample:"
				bam_samples=$bam_samples"$sample.sorted.bam:"
				input_dirs=$input_dirs"$output_dir/alignment/$sample:"
				#bam_files=$bam_files"$output_dir/alignment/$sample/$sample.sorted.bam:"
			done
			realign_dir=$output_dir/realign/$group
			variant_dir=$output_dir/variants/$group
			mkdir -p $realign_dir $variant_dir

			VARIANT=`qsub $args -N $type.$version.realign_recal.$group.$run_num -hold_jid $job_ids_convert -l h_vmem=8G -t 1-$numchrs:1 $script_path/realign_recal.sh $input_dirs $bam_samples $names_samples $realign_dir $variant_dir $run_info 1`
			variant_id=`echo $VARIANT | cut -d ' ' -f3  | tr "\n" "," | sed -e "s/\..*,//g"`
            SPLIT_IGV=`qsub $args -N $type.$version.split_sample_pair.$group.$run_num -hold_jid $variant_id -t 1-$numchrs:1 $script_path/split_sample_pair.sh $output_dir/realign $output_dir/IGV_BAM $group $output_dir/alignment $run_info` 
            split_igv_id=`echo $SPLIT_IGV | cut -d ' ' -f3  | tr "\n" "," | sed -e "s/\..*,//g"`
			IGV=`qsub $args -N $type.$version.igv_bam.$group.$run_num -hold_jid $split_igv_id $script_path/igv_bam.sh $output_dir/realign $output_dir/IGV_BAM $group $output_dir/alignment $run_info`  
			igv_id=`echo $IGV | cut -d ' ' -f3`
			CALLS=`qsub $args -N $type.$version.variants.$group.$run_num -hold_jid $variant_id -pe threaded $threads -l h_vmem=8G -t 1-$numchrs:1 $script_path/variants.sh $realign_dir $names_samples $variant_dir 1 $run_info`
			calls_id=`echo $CALLS | cut -d ' ' -f3  | tr "\n" "," | sed -e "s/\..*,//g"`
			MERGEVARIANT=`qsub $args -N $type.$version.merge_variant_group.$group.$run_num -l h_vmem=8G -hold_jid $calls_id $script_path/merge_variant_group.sh $output_dir/variants $group $output_dir/Reports_per_Sample/ $run_info`
			job_id_mergevariant=`echo $MERGEVARIANT | cut -d ' ' -f3`    
			ONTARGET_B=`qsub $args -N $type.$version.OnTarget_BAM.$group.$run_num -hold_jid $split_igv_id -l h_vmem=8G -t 1-$numchrs:1 $script_path/OnTarget_BAM.sh $output_dir/IGV_BAM $output_dir/OnTarget $group $run_info`
			ONTARGET_P=`qsub $args -N $type.$version.OnTarget_PILEUP.$group.$run_num -hold_jid $split_igv_id -l h_vmem=8G -t 1-$numchrs:1 $script_path/OnTarget_PILEUP.sh $output_dir/IGV_BAM $output_dir/OnTarget $group $run_info`
			job_ids_target_p=`echo $ONTARGET_P | cut -d ' ' -f3 |  tr "\n" "," | sed -e "s/\..*,//g"` 
			job_ids_target_b=`echo $ONTARGET_B | cut -d ' ' -f3 |  tr "\n" "," | sed -e "s/\..*,//g"`
			COVERAGE=`qsub $args -N $type.$version.getCoverage.$group.$run_num -hold_jid ${job_ids_target_p},${job_ids_target_b} $script_path/getCoverage.sh $output_dir/OnTarget $output_dir/numbers $group $run_info`
			ANNOTATE_SNV=`qsub $args -N $type.$version.OnTarget_variant.$sample.$run_num -l h_vmem=4G -hold_jid $job_id_mergevariant -t 1-$numchrs:1 $script_path/OnTarget_variant.sh $output_dir/variants $output_dir/OnTarget $group $run_info`
			job_id_annotate_snv=`echo $ANNOTATE_SNV | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
			hold_args="-hold_jid $job_id_annotate_snv"
			for sample in $samples
			do      
				SIFT=`qsub $args -N $type.$version.sift.$sample.$run_num $hold_args -t 1-$numchrs:1 -l h_vmem=4G $script_path/sift.sh $sift $output_OnTarget $sample $run_info` 
				SSEQ=`qsub $args -N $type.$version.sseq.$sample.$run_num $hold_args -l h_vmem=4G $script_path/sseq.sh $sseq $output_OnTarget $email $sample $run_info`
				echo -e $SSEQ >> $job_ids_dir/ANNOT 
				job_ids_sift=`echo $SIFT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
				job_ids_sseq=`echo $SSEQ | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
				ADD_ANOT=`qsub $args -N $type.$version.sample_reports.$sample.$run_num -hold_jid $job_ids_sseq,$job_ids_sift -t 1-$numchrs:1 -l h_vmem=8G $script_path/sample_reports.sh $run_info $sample $TempReports $output_OnTarget $sift $sseq $output_dir`
				job_ids_anot=`echo $ADD_ANOT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
				MERGE=`qsub $args -N $type.$version.sample_report.$sample.$run_num -hold_jid $job_ids_anot $script_path/sample_report.sh $output_dir $TempReports $sample $run_info`    
				echo -e $MERGE >> $job_ids_dir/PER_SAMPLE	
			done
			sampleNames=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2 )
			i=1
			for sample in $sampleNames
			do
				sampleArray[$i]=$sample
				let i=i+1
			done
			for i in $(seq 2 ${#sampleArray[@]})
			do  
				tumor=${sampleArray[$i]}
				SIFT=`qsub $args -N $type.$version.sift.$tumor.$run_num $hold_args -t 1-$numchrs:1 -l h_vmem=4G $script_path/sift.sh $sift $output_OnTarget $group.$tumor $run_info` 
				SSEQ=`qsub $args -N $type.$version.sseq.$tumor.$run_num $hold_args -l h_vmem=4G $script_path/sseq.sh $sseq $output_OnTarget $email $group.$tumor $run_info`
				echo -e $SSEQ >> $job_ids_dir/ANNOT 
				job_ids_sift=`echo $SIFT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
				job_ids_sseq=`echo $SSEQ | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
				ADD_ANOT=`qsub $args -N $type.$version.sample_reports.$tumor.$run_num  -hold_jid $job_ids_sseq,$job_ids_sift -t 1-$numchrs:1 -l h_vmem=8G $script_path/sample_reports.sh $run_info $group.$tumor $TempReports $output_OnTarget $sift $sseq $output_dir `
				job_ids_anot=`echo $ADD_ANOT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
				MERGE=`qsub $args -N $type.$version.sample_report.$tumor.$run_num -hold_jid $job_ids_anot $script_path/sample_report.sh $output_dir $TempReports $group.$tumor $run_info` 
				echo -e $MERGE >> $job_ids_dir/PER_SAMPLE
			done                    		
			crest=$output_dir/struct/crest
			break=$output_dir/struct/break
			mkdir -p $break
			mkdir -p $crest

			for sam in `cat $sample_info| grep -w "^$group" | cut -d '=' -f2`
            do
                CREST=`qsub $args -N $type.$version.run_crest_multi_cover.$sam.$run_num -hold_jid $split_igv_id -l h_vmem=8G -t 1-$numchrs:1 $script_path/run_crest_multi_cover.sh $sam $group $output_dir/IGV_BAM $crest $run_info`
                job_id_cover=`echo $CREST | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
                job_ids_cover="$job_ids_cover,$job_id_cover"
            done
            CRESTCALL=`qsub $args -N $type.$version.run_crest_multi.$group.$run_num -hold_jid $job_ids_cover -l h_vmem=8G -t 1-$numchrs:1 $script_path/run_crest_multi.sh $group $output_dir/IGV_BAM $crest $run_info`
			job_ids_crest=`echo $CRESTCALL | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g" `    
			SEGSEQCALL=`qsub $args -N $type.$version.run_segseq.$group.$run_num -hold_jid $split_igv_id -l h_vmem=12G -t 1-$numchrs:1 -l matlab_lic=1 $script_path/run_segseq.sh $group $output_dir/IGV_BAM $output_dir/cnv $run_info`
			job_ids_segseq=`echo $SEGSEQCALL | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g" `    
			let nump=$numchrs+1;    
			mkdir -p $break/$group
			for sam in `cat $sample_info| grep -w "^$group" | cut -d '=' -f2`
            do
                BREAKDANCER=`qsub $args -N $type.$version.run_breakdancer.$group.$run_num -hold_jid $split_igv_id -l h_vmem=8G -t 1-$numchrs:1 $script_path/run_breakdancer.sh $sam $output_dir/IGV_BAM $break/$group $run_info $group`
                job_id_break=`echo $BREAKDANCER | cut -d ' ' -f3 |  tr "\n" "," | sed -e "s/\..*,//g"`
                BREAKDANCER_IN=`qsub $args -N $type.$version.run_breakdancer.$group.$run_num -hold_jid $igv_id -l h_vmem=8G -t $nump-$nump:$nump $script_path/run_breakdancer.sh $sam $output_dir/IGV_BAM $break/$group $run_info`
                job_id_break_in=`echo $BREAKDANCER_IN | cut -d ' ' -f3 |  tr "\n" "," | sed -e "s/\..*,//g"`
                job_ids_break="$job_ids_break,$job_id_break"
                job_ids_break_in="$job_ids_break_in,$job_id_break_in"
            done
            MERGESTRUCT=`qsub $args -N $type.$version.summaryze_struct_group.$group.$run_num -hold_jid $job_ids_segseq,$job_ids_break_in,$job_ids_break,$job_ids_crest $script_path/summaryze_struct_group.sh $group $output_dir $run_info`
			job_id_mergestruct=`echo $MERGESTRUCT | cut -d ' ' -f3`
			echo -e $MERGESTRUCT >> $job_ids_dir/SV
			mkdir -p $output_dir/circos;
            for i in $(seq 2 ${#sampleArray[@]})
			do  
				tumor=${sampleArray[$i]}
				RUNCIRCOS=`qsub $args -N $type.$version.plot_circos_cnv_sv.$group.$run_num -hold_jid $job_id_mergestruct -l h_vmem=8G $script_path/plot_circos_cnv_sv.sh $output_dir/struct/$group.$tumor.somatic.break $output_dir/struct/$group.$tumor.somatic.filter.crest $output_dir/cnv/$group/$tumor.cnv.filter.bed $tumor $output_dir/circos $run_info` 
			done
		done

		job_ids_annn=$( cat $job_ids_dir/ANNOT |  cut -d ' ' -f3 | cut -d '.' -f1 | tr "\n" ",")
		job_ids=$( cat $job_ids_dir/SV | cut -d ' ' -f3  | tr "\n" "," )
		mkdir -p $output_dir/Reports_per_Sample/ANNOT
		ANNOTATE_CNV=`qsub $args -N $type.$version.annotation_CNV.$run_num -l h_vmem=4G -hold_jid $job_ids -t 1-$numgroups:1 $script_path/annotation_CNV.sh $output_dir/Reports_per_Sample/SV/ $run_info $output_dir/Reports_per_Sample/ANNOT`
        ANNOTATE_SV=`qsub $args -N $type.$version.annotation_SV.$run_num -l h_vmem=4G -hold_jid $job_ids -t 1-$numgroups:1 $script_path/annotation_SV.sh $output_dir $run_info $output_dir/Reports_per_Sample/ANNOT`
		job_id_annotate_sv=`echo $ANNOTATE_SV | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
		job_id_annotate_cnv=`echo $ANNOTATE_CNV | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
		hold_arg="-hold_jid $job_id_annotate_cnv,$job_id_annotate_sv,$job_ids_annn"  
		### generate reports for all the samples
        REPORT=`qsub $args -N $type.$version.merged_report.$run_num $hold_arg -t 1-$numchrs:1 -l h_vmem=8G $script_path/merged_report.sh $sift $sseq $TempReports $run_info $output_dir/OnTarget` 
		job_ids_report=`echo $REPORT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
		job_ids_sample=$( cat $job_ids_dir/PER_SAMPLE | cut -d ' ' -f3  | tr "\n" ",")
		ANNOT_SAMPLE=`qsub $args -N $type.$version.annotate_sample.$run_num -hold_jid $job_ids_sample -l h_vmem=8G $script_path/annotate_sample.sh $output_dir $run_info`   
		job_ids_annot_sample=`echo $ANNOT_SAMPLE | cut -d ' ' -f3 | tr "\n" ","`
		## annotate merged report
		MERGE_SAMPLE=`qsub $args -N $type.$version.annotate_merged_report.$run_num -hold_jid $job_ids_report -l h_vmem=8G $script_path/annotate_merged_report.sh $output_dir $TempReports $run_info`
		job_ids_merge_sample=`echo $MERGE_SAMPLE | cut -d ' ' -f3 | tr "\n" ","`
		### variant distance
		qsub $args -N $type.$version.variant_distance.$run_num -hold_jid $job_ids_report -l h_vmem=4G $script_path/variant_distance.sh $TempReports $output_dir $run_info
		hold="-hold_jid $job_ids_report,$job_ids_annot_sample,$job_ids_merge_sample"
		
        NUMBERS=`qsub $args -N $type.$version.sample_numbers.$run_num $hold -t 1-$numgroups:1 $script_path/sample_numbers.sh $output_dir $run_info`
#			GENE_SUMMARY=`qsub $args -N $type.$version.gene_summary.$run_num $hold -t 1-$numsamples:1 $script_path/gene_summary.sh $output_dir $run_info $output_dir/Reports_per_Sample`
			job_ids=`echo $NUMBERS | cut -d ' ' -f3 | cut -d '.' -f1 | tr "\n" ","`
			job_ids_summary=`echo $GENE_SUMMARY | cut -d ' ' -f3 | cut -d '.' -f1 | tr "\n" ","`
#			qsub $args -N $type.$version.generate_html.$run_num -hold_jid ${job_ids}${job_ids_summary} $script_path/generate_html.sh $output_dir $run_info 	
		#done
		echo `date`
	fi
fi
