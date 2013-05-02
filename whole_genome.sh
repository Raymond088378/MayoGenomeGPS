#!/bin/sh

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
###
### Version 1.2.1
###  
### TODO Refactor all workflow steps to imperative logic (started 1/29/2013 - Christian Ross)
### TODO Valid Analysis Types
### 	alignment
###		external
###		mayo
###		realignment
###		realign-mayo
###		variant
### Script Dependencies (not complete - CR 2/3/2013)
###		create_folder.sh
###		copy_config.sh
### 	dashboard.sh
###		check_qstat.sh
###
### [[ $analysis == "alignment" || $analysis == "annotation" || $analysis == "external" || $analysis == "mayo" || $analysis == "ontarget" || $analysis == "realignment" || $analysis == "realign-mayo" || $analysis == "variant" ]]

### functions
####
# to validate the local variables
function check_variable()	{
	message=$1
	if [[ "$2" == "" ]] 
	then 
		echo "$message is not set correctly."
		exit 1
	fi		
}
#
####

###  
# check for full path
function check_dir()	{
	if [ $2 == "." ]
	then
		echo -e "$message : should be specified as complete path"
		exit 1;
	fi	
}
#
###				

### 
# to check and validate the configuration file presence
function check_config()	{
	message=$1
	if [ ! -s $2 ]
	then
		echo -e "$message : doesn't exist"
		exit 1;
	fi	
	dir_info=`dirname $2`
	if [ $dir_info == "." ]
	then
		echo -e "$message : should be specified as complete path"
		exit 1;
	fi	
	dos2unix $2
	cat $2 | sed 's/^[ \t]*//;s/[ \t]*$//' > $2.tmp
	mv $2.tmp $2														
}	
#
###

if [ $# != 1 ]
then	
	echo -e "Wrapper script to submit all the jobs for dna-seq workflow \
			\nUsage: ./whole_genome.sh <specify full /path/to/run_info file>";
	exit 1;
fi

#### start of the main wrapper

run_info=`readlink -f $1`
check_config "RunInfo:$run_info" $run_info

### get the path for all the scripts
if [ $JOB_ID ]
then
	script=`qstat -j $JOB_ID | grep -w script_file| awk '{print $NF}'`
	script_path=`dirname $script`
else
	script_path=`dirname $0`
fi
check_dir "WORKFLOW_PATH:$script_path" $script_path 			



## extract paths and set local variables for configuration files
tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
check_config "RunInfo:$tool_info" $tool_info
check_config "RunInfo:$sample_info" $sample_info

### validate the config file
#
$script_path/check_config.pl $run_info > $run_info.configuration_errors.log
if [ `cat $run_info.configuration_errors.log | wc -l` -gt 0 ]
then
	echo "Configuration files are not configured properly: look at the errors \
			in $run_info.configuration_errors.log "
	exit 1;
else
	rm $run_info.configuration_errors.log
fi	
#
###

input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2)
groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2)
tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
type=$( cat $run_info | grep -w '^TOOL' | cut -d '=' -f2|tr "[a-z]" "[A-Z]")
version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
queue=$( cat $tool_info | grep -w '^QUEUE' | sed -e '/QUEUE=/s///g')
gatkqueue=$( cat $tool_info | grep -w '^GATKQUEUE' | sed -e '/GATKQUEUE=/s///g')
run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
all_sites=$( cat $tool_info | grep -w '^EMIT_ALL_SITES' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
aligner=$( cat $run_info | grep -w '^ALIGNER' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]")
numchrs=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
paired=$( cat $run_info | grep -w '^PAIRED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2)
variant_type=$(cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")   
somatic_caller=$(cat $run_info | grep -w '^SOMATIC_CALLER' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]") 
samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
limit=$( cat $tool_info | grep -w '^JOB_LIMIT' | cut -d '=' -f2 )
info=$(cat $run_info | grep -w '^SAMPLEINFORMATION' | cut -d '=' -f2 )
workflow=$( cat $run_info | grep '^TOOL=' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
somatic_calling=$( cat $tool_info | grep -w '^SOMATIC_CALLING' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
stop_after_realignment=$( cat $tool_info | grep -w '^STOP_AFTER_REALIGNMENT' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
 	    	     
### create folders
$script_path/create_folder.sh $run_info
output_dir=$output/$PI/$tool/$run_num
config=$output_dir/config
if [ -f $output_dir/folder_exist.log ]
then
	echo "ERROR: folder : $output_dir already exist"
	exit 1;
fi	

## copy config files
$script_path/copy_config.sh $output_dir $run_info
run_info=$output_dir/run_info.txt
### modify the run info file to use configurations in the output folder and assigning local variable for all the configuration files
add=`date +%D`
cat $run_info | grep -w -v -E '^TOOL_INFO|^SAMPLE_INFO|^MEMORY_INFO' > $run_info.tmp
echo -e "\nTOOL_INFO=$config/tool_info.txt\nSAMPLE_INFO=$config/sample_info.txt\nMEMORY_INFO=$config/memory_info.txt\nDATE=$add" \
	| cat $run_info.tmp - > $config/run_info.txt
rm $run_info.tmp
run_info=$config/run_info.txt
tool_info=$config/tool_info.txt
sample_info=$config/sample_info.txt
memory_info=$config/memory_info.txt
###

### creating local variables
output_align=$output_dir/alignment
if [ $analysis != "alignment" ]
then
	output_realign=$output_dir/realign
	output_variant=$output_dir/variants
	igv=$output_dir/IGV_BAM
	RSample=$output_dir/Reports_per_Sample
	if [ $tool == "whole_genome" ]
	then
		sv=$output_dir/Reports_per_Sample/SV
		struct=$output_dir/struct
		cnv=$output_dir/cnv
		circos=$output_dir/circos
	fi
	numbers=$output_dir/numbers
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
email=`finger $USER | awk -F ';' '{print $2}' | head -n1`
args="-V -wd $output_dir/logs -q $queue -m a -M $email -l h_stack=10M"
### Test if gatkqueue is empty; if so, fail through to using standard args for queue
if [[ "${gatkqueue+xxx}" = "xxx" ]]
then
	gatk_args=$args
	echo "Failed to find GATK specific queue, using $queue"
else 
	gatk_args="-V -wd $output_dir/logs -q $gatkqueue -m a -M $email -l h_stack=10M"
fi

echo -e "\nRCF arguments used : $args\n$gatk_args\n" >> $output_dir/log.txt
echo -e "Started the ${tool} analysis for ${run_num} for ${PI}\n\n${info}\n\nCourtesy: $workflow $version" | mailx -v -s "Analysis Started" "$email"
#############################################################

### Single sample workflow
if [ $multi_sample != "YES" ]
then
	echo "Single sample"
	numsamples=$(cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
	for sample in `echo $samples | tr ":" "\n"`
	do     
		### 	alignment external mayo realignment realign-mayo variant
		align_dir=$output_dir/alignment/$sample
		bamfile=$sample.sorted.bam
		if [[ $analysis == "mayo" || $analysis == "external" || $analysis == "alignment" ]]
		then
			mkdir -p $align_dir
			if [[ $paired == 1  || $paired == "YES" ]]
			then
				let numfiles=(`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" |wc -l`)/2
			else
				let numfiles=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" |wc -l`
			fi	
			if [ $numfiles -eq 0 ]
			then
				echo "Sample info :$sample_info file is not properly configured, number of files in sample_info is not defined for $sample"
				exit 1;
			fi	
			### Invoke dashboard for mayo internal job tracking only
			if [ $analysis == "mayo" ]
			then
				for ((i=1; i <=$numfiles; i++));
				do
					$script_path/dashboard.sh $sample $run_info Beginning started $i
				done
			fi
			#### alignment	
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
				echo "Workflow doesn't support the aligner: $aligner"
			fi	
			if [ $aligner == "bwa" ]
			then
				hold="-hold_jid $type.$version.align_bwa.$sample.$run_num.$identify"
			elif [ $aligner == "novoalign" ]
			then
				hold="-hold_jid $type.$version.align_novo.$sample.$run_num.$identify"
			fi    
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^convert_sam_bam' | cut -d '=' -f2)
			for i in $(seq 1 $numfiles)
			do 
				qsub_args="-N $type.$version.convert_sam_bam.$sample.$run_num.$identify $hold -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/convert_sam_bam.sh $output_dir/$sample $sample.$i.bam $sample.$i $run_info $i
			done
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^processBAM' | cut -d '=' -f2)
			qsub_args="-N $type.$version.processBAM.$sample.$run_num.$identify $type.$version.convert_sam_bam.$sample.$run_num.$identify \
						-pe threaded $threads -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/processBAM.sh $align_dir $sample $run_info
			### mayo, external  	
			### not just doing alignment so split the bam files and put unmapped reads in $bamfile.extra.bam 
			if [ $analysis != "alignment" ]
			then
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^extract_reads_bam' | cut -d '=' -f2)
				qsub_args="-N $type.$version.extract_reads_bam.$sample.$run_num.$identify \
						-hold_jid $type.$version.processBAM.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/extract_reads_bam.sh $align_dir $bamfile $run_info $igv
			fi
		#########################################
		### 	Realignment Entry Point		#####
		#########################################
		elif [[ $analysis == "realignment" || $analysis == "realign-mayo" ]]
		then
			mkdir -p $align_dir
			infile=`cat $sample_info | grep -w ^BAM:${sample} | cut -d '=' -f2`
			num_bams=`echo $infile | tr " " "\n" | wc -l`
			if [ $num_bams -eq 0 ]
			then
				echo "sample info : $sample_info file is not properly configured"
				exit 1;
			fi
			### Mayo internal, call the dashboard
			if [ $analysis == "realign-mayo" ]
			then
				$script_path/dashboard.sh $sample $run_info Beginning started 
			fi
			
			### sort each bam file
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
				#### run fastc on bam files
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^fastqc' | cut -d '=' -f2)
				qsub_args="-N $type.$version.fastqc.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/fastqc.sh $fastqc $align_dir/$sample.$i.sorted.bam $tool_info
			done  
			### Sort, index, deduplicate the bam files 
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^processBAM' | cut -d '=' -f2)
			qsub_args="-N $type.$version.processBAM.$sample.$run_num.$identify \
					-hold_jid $type.$version.fastqc.$sample.$run_num.$identify -l h_vmem=$mem pe threaded $threads"
			qsub $args $qsub_args $script_path/processBAM.sh $align_dir $sample $run_info
			### Split on chromsomes and store unmapped reads for igv
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^extract_reads_bam' | cut -d '=' -f2)
			qsub_args="-N $type.$version.extract_reads_bam.$sample.$run_num.$identify \
						-hold_jid $type.$version.processBAM.$sample.$run_num.$identify -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/extract_reads_bam.sh $align_dir $bamfile $run_info $igv
		fi    
		############################
		### VARIANT  ENTRY POINT ###
		############################
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
					echo "sample info :$sample_info file is not properly configured"
					exit 1;
				fi
				### extract headers for each bam file and create symbolic link to the sorted bams (that will be?) present in the realign directory
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
					#### run fastc on bam files
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^fastqc' | cut -d '=' -f2)
					qsub_args="-N $type.$version.fastqc.$sample.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/fastqc.sh $fastqc $realign_dir/$sample.$i.sorted.bam $tool_info				
				done
				### Merge bam files within sample group, and create sample.sorted.bam in the realign dir
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^reformat_BAM' | cut -d '=' -f2)
				qsub_args="-N $type.$version.reformat_BAM.$sample.$run_num.$identify \
							-hold_jid $type.$version.fastqc.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/reformat_BAM.sh $realign_dir $sample $run_info	
				$script_path/check_qstat.sh $limit
				### Extract Unmapped Reads for IGV
				mem=$( cat $memory_info | grep -w '^extract_reads_bam' | cut -d '=' -f2)
				qsub_args="-N $type.$version.extract_reads_bam.$sample.$run_num.$identify \
							-hold_jid $type.$version.reformat_BAM.$sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/extract_reads_bam.sh $realign_dir $bamfile $run_info $igv
				$script_path/check_qstat.sh $limit
				#### Split each sample.sorted.bam file by chromosome and store in the realign dir
				mem=$( cat $memory_info | grep -w '^split_bam_chr' | cut -d '=' -f2)
				qsub_args="-N $type.$version.split_bam_chr.$sample.$run_num.$identify \
						-hold_jid $type.$version.reformat_BAM.$sample.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/split_bam_chr.sh $realign_dir $sample $run_info
				variant_id="$type.$version.split_bam_chr.$sample.$run_num.$identify"
			else
				### Not variant, must be mayo, external, realignment or realign-mayo -- so run realignment and recalibration
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^realign_recal' | cut -d '=' -f2)
				qsub_args="-N $type.$version.realign_recal.$sample.$run_num.$identify \
						-hold_jid $type.$version.processBAM.$sample.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $gatk_args $qsub_args $script_path/realign_recal.sh $align_dir $bamfile $sample $realign_dir $run_info 1	
				variant_id="$type.$version.realign_recal.$sample.$run_num.$identify"
			fi
			### IGV_BAM
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^igv_bam' | cut -d '=' -f2)
			qsub_args="-N $type.$version.igv_bam.$sample.$run_num.$identify \
					-hold_jid $variant_id,$type.$version.extract_reads_bam.$sample.$run_num.$identify -pe threaded $threads -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/igv_bam.sh $output_realign $igv $sample $output_align $run_info
			### If we're just doing realignment, then jump to the next sample without doing anything else
			if [ $stop_after_realignment == "YES" ]
            then
            	continue;
            fi
            
            ###############################
			### VARIANT CALLING SECTION ###
			###############################
            $script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^variants' | cut -d '=' -f2)
			qsub_args="-N $type.$version.variants.$sample.$run_num.$identify -hold_jid $variant_id -pe threaded $threads -t 1-$numchrs:1 -l h_vmem=$mem"
			qsub $gatk_args $qsub_args $script_path/variants.sh $realign_dir $sample $variant_dir 1 $run_info
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^merge_variant_single' | cut -d '=' -f2)
			qsub_args="-N $type.$version.merge_variant_single.$sample.$run_num.$identify -pe threaded $threads \
					-hold_jid $type.$version.variants.$sample.$run_num.$identify -l h_vmem=$mem"
			qsub $gatk_args $qsub_args $script_path/merge_variant_single.sh $output_variant $sample $RSample $run_info
			$script_path/check_qstat.sh $limit
			### ODOT 

			mem=$( cat $memory_info | grep -w '^OnTarget_BAM' | cut -d '=' -f2)
			qsub_args="-N $type.$version.OnTarget_BAM.$sample.$run_num.$identify -hold_jid $variant_id -t 1-$numchrs:1 -l h_vmem=$mem"
			qsub $args $qsub_args  $script_path/OnTarget_BAM.sh $realign_dir $output_OnTarget $sample $run_info
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^OnTarget_PILEUP' | cut -d '=' -f2)
			qsub_args="-N $type.$version.OnTarget_PILEUP.$sample.$run_num.$identify -hold_jid $variant_id -t 1-$numchrs:1 -l h_vmem=$mem"
			qsub $gatk_args $qsub_args $script_path/OnTarget_PILEUP.sh $realign_dir $output_OnTarget $sample $run_info
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^getCoverage' | cut -d '=' -f2)
			qsub_args="-N $type.$version.getCoverage.$sample.$run_num.$identify -hold_jid $type.$version.OnTarget_PILEUP.$sample.$run_num.$identify -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/getCoverage.sh $output_OnTarget $numbers $sample $run_info
        ### END mayo external realignment variant realign-mayo SECTION
		fi

		
	### done ### FOR EACH SAMPLE 
	### Closing the main loop, reopen after backfill

	### best dependency at this point
	### $type.$version.variants.$sample (*).$run_num.$identify
	

	#################################
	### Single-Sample Backfilling ###
	#################################
	
	### TODO Set up synchronized merge task dependent on variants
 
	###	that we haven't broken the pipeline dependencies between variants and OnTarget				
	
		
	### Reopening the job submission loop	
	### Loop and submit merge variant single, set future dependencies on it
	
	### for sample in `echo $samples | tr ":" "\n"`
	### do 
	### 	if [[ $analysis == "mayo" || $analysis == "external" || $analysis == "realignment" || $analysis == "variant" || $analysis == "realign-mayo" ]]
	### 	then
	###			### Depend on backfill tasks
	###			qsub_args="-N $type.$version.merge_variant_single.$sample.$run_num.$identify -pe threaded $threads -hold_jid $type.$version.variants.$sample.$run_num.$identify -l h_vmem=$mem"
	###			qsub $gatk_args $qsub_args $script_path/merge_variant_single.sh $output_variant $sample $RSample $run_info
	###			$script_path/check_qstat.sh $limit
		
				

				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^snpeff' | cut -d '=' -f2)
				qsub_args="-N $type.$version.snpeff.$sample.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $gatk_args $qsub_args $script_path/snpeff.sh $snpeff $output_OnTarget $sample $run_info germline		
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
			qsub $args $qsub_args $script_path/run_cnvnator.sh $sample $realign_dir $cnv $run_info single
			let nump=$numchrs+1;
			
			## Breakdancer Removed for 2.0 
				
			# $script_path/check_qstat.sh $limit
			# mem=$( cat $memory_info | grep -w '^run_breakdancer' | cut -d '=' -f2)
			# qsub_args="-N $type.$version.run_breakdancer.$sample.$run_num.$identify -hold_jid $variant_id -t 1-$numchrs:1 -l h_vmem=$mem"
			# qsub $args $qsub_args $script_path/run_breakdancer.sh $sample $output_realign $break $run_info single
			
			# $script_path/check_qstat.sh $limit
			# qsub_args="-N $type.$version.run_breakdancer_in.$sample.$run_num.$identify -hold_jid $type.$version.igv_bam.$sample.$run_num.$identify -t $nump-$nump:$nump -l h_vmem=$mem"
			# qsub $args $qsub_args $script_path/run_breakdancer.sh $sample $igv $break $run_info single
			
			### merge the structural variants
			hold="-hold_jid $type.$version.run_single_crest.$sample.$run_num.$identify,$type.$version.run_cnvnator.$sample.$run_num.$identify,"
	
			# hold=$hold"$type.$version.run_breakdancer.$sample.$run_num.$identify,$type.$version.run_breakdancer_in.$sample.$run_num.$identify"
			mkdir -p $sv
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^summaryze_struct_single' | cut -d '=' -f2)
			qsub_args="-N $type.$version.summaryze_struct_single.$sample.$run_num.$identify -l h_vmem=$mem $hold"
			qsub $gatk_args $qsub_args $script_path/summaryze_struct_single.sh $sample $output_dir $run_info
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^plot_circos_cnv_sv' | cut -d '=' -f2)
			qsub_args="-N $type.$version.plot_circos_cnv_sv.$sample.$run_num.$identify -hold_jid $type.$version.summaryze_struct_single.$sample.$run_num.$identify -l h_vmem=$mem"
			
			break_file=$break/$sample/$sample.break
			crest_file=$crest/$sample/$sample.final.crest
			cnv_file=$cnv/$sample.cnv.final.bed
			qsub $args $qsub_args $script_path/plot_circos_cnv_sv.sh $break_file $crest_file $cnv_file $sample $circos $run_info	
		fi
		if [ $annot_flag == "YES" ]
		then
			if [[ $tool == "whole_genome" && $analysis != "alignment" && $analysis != "annotation" && $analysis != "ontarget" ]]
			then
				mkdir -p $output_dir/Reports_per_Sample/ANNOT
				### Concatenante CNV outputs, awk.
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
		if [ $annot_flag == "YES" ]
		then
        	### For all analyses, run sample_numbers
        	$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^sample_numbers' | cut -d '=' -f2)
			qsub_args="-N $type.$version.sample_numbers.$sample.$run_num.$identify $hold_args -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/sample_numbers.sh $output_dir $sample $run_info $numbers
		fi
	done ### FOR EACH SAMPLE 
	
	#########################
	### MERGE CHROMOSOMES ###
	#########################

	### concat raw varaints
	if [ $stop_after_realignment == "NO" ]
    then
    	if [[  $tool == "exome"  && $all_sites == "YES" ]]
		then
			### Build hold on variants for merge job
			id=""
			for s in `echo $samples | tr ":" "\n"`
			do
				id=$id"$type.$version.variants.$s.$run_num.$identify,"
			done
			### array on chromosomes, merge each chromosome across all samples
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^merge_raw_variants' | cut -d '=' -f2)
			qsub_args="-N $type.$version.merge_raw_variants.$run_num.$identify -t 1-$numchrs:1 -hold_jid $id -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/merge_raw_variants.sh $output_dir $run_info
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^concat_raw_variants' | cut -d '=' -f2)
			qsub_args="-N $type.$version.concat_raw_variants.$run_num.$identify -hold_jid $type.$version.merge_raw_variants.$run_num.$identify -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/concat_raw_variants.sh $output_dir $run_info	    
		fi	
    fi
	if [ $stop_after_realignment == "NO" ]
	then
		if [ $analysis != "alignment" ]
		then
			id=""
			for s in `echo $samples | tr ":" "\n"`
			do
				id=$id"$type.$version.sample_report.$s.$run_num.$identify,"
			done 
			if [ $annot_flag == "YES" ]
			then
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^merge_sample' | cut -d '=' -f2)
				qsub_args="-N $type.$version.merge_sample.$run_num.$identify -hold_jid $id -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/merge_sample.sh $output_dir $run_info    
			fi
		fi
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
	if [ $stop_after_realignment == "NO" ]
	then
        if [ $annot_flag == "YES" ]
		then
		$script_path/check_qstat.sh $limit
		mem=$( cat $memory_info | grep -w '^generate_html' | cut -d '=' -f2)
		qsub_args="-N $type.$version.generate_html.$run_num.$identify $hold -l h_vmem=$mem"
		qsub $args $qsub_args $script_path/generate_html.sh $output_dir $run_info
	fi
	fi
    #####################################
    ### End of Single Sample Workflow ###
	#####################################
else
	### Multiple Samples
	echo "Multi-sample"
	numgroups=$(cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
	
	##### alignment and processing of the reads per sample is started
	for sample in `echo $samples | tr ":" "\n"`
	do
		bamfile=$sample.sorted.bam
		align_dir=$output_dir/alignment/$sample;
		mkdir -p $align_dir
		if [[ $analysis == "mayo" || $analysis == "external" || $analysis == "alignment" ]]
		then
			if [ $paired == 1 ]
			then
				let numfiles=(`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" |wc -l`)/2
			else
				let numfiles=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" |wc -l`
			fi
			
			if [ $numfiles -eq 0 ]
			then
				echo "sample info : $sample_info file is not properly configured, please correct it and start the workflow again"
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
			else 
				echo "Workflow doesn't support the aligner: $aligner"
				### ENDIF ALIGNERS
			fi	    
			
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^convert_sam_bam' | cut -d '=' -f2)
			for i in $(seq 1 $numfiles)
			do 
				qsub_args="-N $type.$version.convert_sam_bam.$sample.$run_num.$identify $hold -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/convert_sam_bam.sh $output_dir/$sample $sample.$i.bam $sample.$i $run_info $i
			done
			
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^processBAM' | cut -d '=' -f2)
			qsub_args="-N $type.$version.processBAM.$sample.$run_num.$identify -hold_jid $type.$version.convert_sam_bam.$sample.$run_num.$identify -l h_vmem=$mem"
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
				echo "sample info : $sample_info file is not properly configured, no bam files found for realignment"
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
	##### alignment and processing the reads per sample is completed and bam file is ready
	
	### For All Other Analyses, other than alignment
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
				### setting the bams and its path for multiple sample analysis
				names_samples=$names_samples"$sample:"
				bam_samples=$bam_samples"$sample.sorted.bam:"
				input_dirs=$input_dirs"$output_dir/alignment/$sample:"
			done
            
            bam_samples=`echo $bam_samples | sed '$s/.$//'`
            input_dirs=`echo $input_dirs | sed '$s/.$//'`
            names_samples=`echo $names_samples | sed '$s/.$//'`
                        
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
					echo "sample info file is not properly configured, there are no bam files for $group"
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
				samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" " ")
				for ss in $sample
				do
					qsub $args $qsub_args $script_path/extract_reads_bam.sh $realign_dir $group.sorted.bam $run_info $igv $group $ss
				done
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^split_bam_chr' | cut -d '=' -f2)
				qsub_args="-N $type.$version.split_bam_chr.$group.$run_num.$identify -hold_jid $type.$version.reformat_pairBAM.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/split_bam_chr.sh $realign_dir $group $run_info
				
				variant_id="$type.$version.split_bam_chr.$group.$run_num.$identify,$type.$version.extract_reads_bam.$group.$run_num.$identify"
				vvid="$type.$version.split_bam_chr.$group.$run_num.$identify,$type.$version.extract_reads_bam.$group.$run_num.$identify"
			else ### analysis is  mayo, mayo-realign, external, annotation       
				
				### set holds on all the extract read processes 
				vvid=""
				for sample in $samples
				do
					vvid=$vvid"$type.$version.processBAM.$sample.$run_num.$identify,$type.$version.extract_reads_bam.$sample.$run_num.$identify"
				done    
				
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^realign_recal' | cut -d '=' -f2)
				qsub_args="-N $type.$version.realign_recal.$group.$run_num.$identify -hold_jid $vvid -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $gatk_args $qsub_args $script_path/realign_recal.sh $input_dirs $bam_samples $names_samples $realign_dir $run_info 1
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
			
			if [ $stop_after_realignment == "YES" ]
            then
            	continue;
			fi    
			
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^variants' | cut -d '=' -f2)
			qsub_args="-N $type.$version.variants.$group.$run_num.$identify -hold_jid $variant_id -pe threaded $threads -t 1-$numchrs:1 -l h_vmem=$mem"
			qsub $gatk_args $qsub_args $script_path/variants.sh $realign_dir $names_samples $variant_dir 1 $run_info
			
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^merge_variant_group' | cut -d '=' -f2)
			qsub_args="-N $type.$version.merge_variant_group.$group.$run_num.$identify -pe threaded $threads -hold_jid $type.$version.variants.$group.$run_num.$identify -l h_vmem=$mem"
			qsub $gatk_args $qsub_args $script_path/merge_variant_group.sh $output_variant $group $RSample $run_info 
			
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^OnTarget_BAM' | cut -d '=' -f2)
			qsub_args="-N $type.$version.OnTarget_BAM.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/OnTarget_BAM.sh $igv $output_OnTarget $group $run_info
			
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^OnTarget_PILEUP' | cut -d '=' -f2)
			qsub_args="-N $type.$version.OnTarget_PILEUP.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
			qsub $gatk_args $qsub_args $script_path/OnTarget_PILEUP.sh $realign_dir $output_OnTarget $group $run_info
			
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^getCoverage' | cut -d '=' -f2)
			qsub_args="-N $type.$version.getCoverage.$group.$run_num.$identify -hold_jid $type.$version.OnTarget_PILEUP.$group.$run_num.$identify -l h_vmem=$mem"
			qsub $args $qsub_args $script_path/getCoverage.sh $output_OnTarget $numbers $group $run_info
			
			$script_path/check_qstat.sh $limit
			mem=$( cat $memory_info | grep -w '^OnTarget_variant' | cut -d '=' -f2)
			qsub_args="-N $type.$version.OnTarget_variant.$group.$run_num.$identify -hold_jid $type.$version.merge_variant_group.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
			qsub $gatk_args $qsub_args $script_path/OnTarget_variant.sh $output_variant $output_OnTarget $group $run_info
			
			hold_args="-hold_jid $type.$version.OnTarget_variant.$group.$run_num.$identify"
			
			### annotation for the called variants started
			## SIFT
			if [ $annot_flag == "YES" ]
            then
            	$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^sift' | cut -d '=' -f2)
				qsub_args="-N $type.$version.sift.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/sift.sh $sift $output_OnTarget $group $run_info germline
				if [ $somatic_calling == "YES" ]
				then	
					$script_path/check_qstat.sh $limit
					qsub_args="-N $type.$version.sift.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/sift.sh $sift $output_OnTarget $group $run_info somatic
				fi
			
				## SNPEFF
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^snpeff' | cut -d '=' -f2)
				qsub_args="-N $type.$version.snpeff.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $gatk_args $qsub_args $script_path/snpeff.sh $snpeff $output_OnTarget $group $run_info germline
				if [ $somatic_calling == "YES" ]
				then
					$script_path/check_qstat.sh $limit
					qsub_args="-N $type.$version.snpeff.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/snpeff.sh $snpeff $output_OnTarget $group $run_info somatic
				fi
			
				##POLYPHEN
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^polyphen' | cut -d '=' -f2)
				qsub_args="-N $type.$version.polyphen.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $gatk_args $qsub_args $script_path/polyphen.sh $polyphen $output_OnTarget $group $run_info germline  
				if [ $somatic_calling == "YES" ]
				then
					$script_path/check_qstat.sh $limit
					qsub_args="-N $type.$version.polyphen.$group.$run_num.$identify $hold_args -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/polyphen.sh $polyphen $output_OnTarget $group $run_info somatic
				fi ### somatic calling
				
				hold="$type.$version.sift.$group.$run_num.$identify,$type.$version.snpeff.$group.$run_num.$identify,$type.$version.polyphen.$group.$run_num.$identify"
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^reports' | cut -d '=' -f2)
				qsub_args="-N $type.$version.reports.$group.$run_num.$identify -hold_jid $hold -t 1-$numchrs:1 -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/reports.sh $run_info $group $TempReports $output_OnTarget $sift $snpeff $polyphen $output_dir germline
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^sample_report' | cut -d '=' -f2)
				qsub_args="-N $type.$version.sample_report.$group.$run_num.$identify -hold_jid $type.$version.reports.$group.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/sample_report.sh $output_dir $TempReports $group $run_info germline
				if [ $somatic_calling == "YES" ]
				then
					hold="$type.$version.sift.$group.$run_num.$identify,$type.$version.snpeff.$group.$run_num.$identify,$type.$version.polyphen.$group.$run_num.$identify"
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^reports' | cut -d '=' -f2)
					qsub_args="-N $type.$version.reports.$group.$run_num.$identify -hold_jid $hold -t 1-$numchrs:1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/reports.sh $run_info $group $TempReports $output_OnTarget $sift $snpeff $polyphen $output_dir somatic
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^sample_report' | cut -d '=' -f2)
					qsub_args="-N $type.$version.sample_report.$group.$run_num.$identify -hold_jid $type.$version.reports.$group.$run_num.$identify -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/sample_report.sh $output_dir $TempReports $group $run_info somatic
				fi ### somatic calling  
			fi ### annotation flag
			### structural variant calls
			if [ $tool == "whole_genome" ]
			then
				crest=$output_dir/struct/crest
				break=$output_dir/struct/break
				mkdir -p $break
				mkdir -p $crest
				for sam in `cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" " "`
				do
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^run_crest_multi_cover' | cut -d '=' -f2)
					if [ $analysis == "variant" ]
					then
						qsub_args="-N $type.$version.run_crest_multi_cover.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
						qsub $args $qsub_args $script_path/run_crest_multi_cover.sh $sam $group $igv $crest $run_info
					else
						qsub_args="-N $type.$version.run_crest_multi_cover.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
						qsub $args $qsub_args $script_path/run_crest_multi_cover.sh $sam $group $output_align/$sam/ $crest $run_info
					fi	
				done
				
				if [[ $somatic_calling == "YES" ]]
				then
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^run_crest_multi' | cut -d '=' -f2)
					qsub_args="-N $type.$version.run_crest_multi.$group.$run_num.$identify -hold_jid $type.$version.run_crest_multi_cover.$group.$run_num.$identify -t 1-$numchrs:1 -pe threaded 2 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/run_crest_multi.sh $group $igv $crest $run_info
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^run_segseq' | cut -d '=' -f2)
					qsub_args="-N $type.$version.run_segseq.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l matlab_lic=1 -l h_vmem=$mem"
					qsub $args $qsub_args $script_path/run_segseq.sh $group $igv $cnv $run_info    
				else
					num=`cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" "\n" | wc -l`
					for flag in $(seq 1 $num)
					do
						$script_path/check_qstat.sh $limit
						mem=$( cat $memory_info | grep -w '^run_cnvnator' | cut -d '=' -f2)
						qsub_args="-N $type.$version.run_cnvnator.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
						qsub $args $qsub_args $script_path/run_cnvnator.sh $group $igv $cnv/$group $run_info $flag
					done
				fi
				
				let nump=$numchrs+1;    
				
				# Breakdancer removed for 2.0
				# mkdir -p $break/$group
				# for sam in `cat $sample_info| grep -w "^$group" | cut -d '=' -f2`
				# do
				#	$script_path/check_qstat.sh $limit
				#	mem=$( cat $memory_info | grep -w '^run_breakdancer' | cut -d '=' -f2)
				# 	qsub_args="-N $type.$version.run_breakdancer.$group.$run_num.$identify -hold_jid $type.$version.split_sample_pair.$group.$run_num.$identify -t 1-$numchrs:1 -l h_vmem=$mem"
				#	qsub $args $qsub_args $script_path/run_breakdancer.sh $sam $igv $break/$group $run_info $group
				#	$script_path/check_qstat.sh $limit
				#	qsub_args="-N $type.$version.run_breakdancer.$group.$run_num.$identify -hold_jid $type.$version.igv_bam.$group.$run_num.$identify -t $nump-$nump:$nump -l h_vmem=$mem"
				#	qsub $args $qsub_args $script_path/run_breakdancer.sh $sam $igv $break/$group $run_info $group
				#done

				if [ $somatic_calling == "YES" ]
				then
					hhold="$type.$version.run_segseq.$group.$run_num.$identify,$type.$version.run_crest_multi.$group.$run_num.$identify,"
					# hhold=$hhold"$type.$version.run_breakdancer.$group.$run_num.$identify,"
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^summaryze_struct_group' | cut -d '=' -f2)
					qsub_args="-N $type.$version.summaryze_struct_group.$group.$run_num.$identify -hold_jid $hhold -l h_vmem=$mem"
					qsub $gatk_args $qsub_args $script_path/summaryze_struct_group.sh $group $output_dir $run_info
				else
					hhold="$type.$version.run_cnvnator.$group.$run_num.$identify,$type.$version.run_breakdancer.$group.$run_num.$identify,$type.$version.run_crest_multi_cover.$group.$run_num.$identify,"
					$script_path/check_qstat.sh $limit
					mem=$( cat $memory_info | grep -w '^summaryze_struct_single' | cut -d '=' -f2)
					qsub_args="-N $type.$version.summaryze_struct_single.$group.$run_num.$identify -hold_jid $hhold -l h_vmem=$mem"
					for sam in `cat $sample_info| grep -w "^$group" | cut -d '=' -f2`
					do
						qsub $args $qsub_args $script_path/summaryze_struct_single.sh $sam $output_dir $run_info $group
					done
				fi
				
				mkdir -p $output_dir/circos;
				if [ $somatic_calling == "YES" ]
				then
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
						qsub_args="-N $type.$version.plot_circos_cnv_sv.$group.$run_num.$identify -hold_jid $type.$version.summaryze_struct_group.$group.$run_num.$identify -l h_vmem=$mem"
						crest_file=$struct/$group.$tumor.somatic.final.crest
						cnv_file=$cnv/$group/$tumor.cnv.final.bed
						# Breakdancer removed 2.0
						# break_file=$struct/$group.$tumor.somatic.break
						# qsub $args $qsub_args $script_path/plot_circos_cnv_sv.sh $break_file $crest_file $cnv_file $group.$tumor $circos $run_info
						qsub $args $qsub_args $script_path/plot_circos_cnv_sv.sh $crest_file $cnv_file $group.$tumor $circos $run_info
					done
				else
					for sample in `cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" " "`
					do
						$script_path/check_qstat.sh $limit
						mem=$( cat $memory_info | grep -w '^plot_circos_cnv_sv' | cut -d '=' -f2)
						qsub_args="-N $type.$version.plot_circos_cnv_sv.$group.$run_num.$identify -hold_jid $type.$version.summaryze_struct_single.$group.$run_num.$identify -l h_vmem=$mem"
						
						crest_file=$struct/crest/$group/$sample.final.crest
						cnv_file=$cnv/$group/$sample.cnv.final.bed
						## Breakdancer removed 2.0
						# break_file=$struct/break/$group/$sample.break
						# qsub $args $qsub_args $script_path/plot_circos_cnv_sv.sh $break_file $crest_file $cnv_file $sample $circos $run_info
						qsub $args $qsub_args $script_path/plot_circos_cnv_sv.sh $crest_file $cnv_file $sample $circos $run_info
					done
				fi	
			fi
		done
		
		mkdir -p $annot
		if [ $stop_after_realignment == "NO" ]
		then
			if [ $annot_flag == "YES" ]
        then
        	if [ $tool == "whole_genome" ]
        	then
					if [ $somatic_calling == "YES" ]
					then
						for group in `echo $groups | tr ":" "\n"`
						do
							$script_path/check_qstat.sh $limit
							mem=$( cat $memory_info | grep -w '^annotation_CNV' | cut -d '=' -f2)
							qsub_args="-N $type.$version.annotation_CNV.$group.$run_num.$identify -hold_jid $type.$version.plot_circos_cnv_sv.$group.$run_num.$identify -l h_vmem=$mem"
							qsub $args $qsub_args $script_path/annotation_CNV.sh $sv $run_info $annot $group
							$script_path/check_qstat.sh $limit
							mem=$( cat $memory_info | grep -w '^annotation_SV' | cut -d '=' -f2)
							qsub_args="-N $type.$version.annotation_SV.$group.$run_num.$identify -hold_jid $type.$version.plot_circos_cnv_sv.$group.$run_num.$identify -l h_vmem=$mem"
							qsub $args $qsub_args $script_path/annotation_SV.sh $output_dir $run_info $annot $group
						done
					else
						### SOMATIC_CALLING != "YES" 
						for group in `echo $groups | tr ":" "\n"`
						do
							for sample in `cat $sample_info| grep -w "^$group" | cut -d '=' -f2 | tr "\t" " "`
							do
								$script_path/check_qstat.sh $limit
								mem=$( cat $memory_info | grep -w '^annotation_CNV' | cut -d '=' -f2)
								qsub_args="-N $type.$version.annotation_CNV.$group.$run_num.$identify -hold_jid $type.$version.plot_circos_cnv_sv.$group.$run_num.$identify -l h_vmem=$mem"
								qsub $args $qsub_args $script_path/annotation_CNV.sh $sv $run_info $annot $sample
								$script_path/check_qstat.sh $limit
								mem=$( cat $memory_info | grep -w '^annotation_SV' | cut -d '=' -f2)
								qsub_args="-N $type.$version.annotation_SV.$group.$run_num.$identify -hold_jid $type.$version.plot_circos_cnv_sv.$group.$run_num.$identify -l h_vmem=$mem"
								qsub $args $qsub_args $script_path/annotation_SV.sh $output_dir $run_info $annot $sample
							done
						done		
					fi	
					### ENDIF SOMATIC_CALLING
				fi
				### ENDIF WHOLE GENOME
		fi
			### ENDIF ANNOTATION FLAG
		
			### generate reports for all the samples
			id=""
			for group in `echo $groups | tr ":" "\n"`
			do
				id=$id"$type.$version.sample_report.$group.$run_num.$identify,$type.$version.OnTarget_PILEUP.$group.$run_num.$identify,$type.$version.OnTarget_BAM.$group.$run_num.$identify,"	
			done
			
			if [ $tool == "whole_genome" ]
			then
				for group in `echo $groups | tr ":" "\n"`
				do
					if [ $somatic_calling == "YES" ]
					then
						id=$id"$type.$version.summaryze_struct_group.$group.$run_num.$identify,$type.$version.annotation_CNV.$group.$run_num.$identify,$type.$version.annotation_SV.$group.$run_num.$identify,"
					else
						id=$id"$type.$version.summaryze_struct_single.$group.$run_num.$identify,$type.$version.annotation_CNV.$group.$run_num.$identify,$type.$version.annotation_SV.$group.$run_num.$identify,"
					fi	
				done
			fi    
			
			if [ $annot_flag == "YES" ]
            then
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
				### ENDIF Exome or WholeGenome Tool

				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^merge_sample' | cut -d '=' -f2)
				qsub_args="-N $type.$version.merge_sample.$run_num.$identify -hold_jid $id -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/merge_sample.sh $output_dir $run_info	
				$script_path/check_qstat.sh $limit
				mem=$( cat $memory_info | grep -w '^generate_html' | cut -d '=' -f2)
				qsub_args="-N $type.$version.generate_html.$run_num.$identify -hold_jid $type.$version.merge_sample.$run_num.$identify -l h_vmem=$mem"
				qsub $args $qsub_args $script_path/generate_html.sh $output_dir $run_info
			fi
			### ENDIF Annotation Flag Set "YES"
		fi
		### ENDIF Stop After Realignment is "NO" 
	fi
	### ENDIF Analysis is not "alignment" 
fi
### ENDIF Multiple Samples
echo `date`
#### end of the script
