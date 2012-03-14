#!/bin/sh

########################################################
###### 	ALIGNMENT SCRIPT FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			bwa.splitted.PR.sh
######		Date:				07/25/2011
######		Summary:			Alignment done using BWA on PE FASTQ files and realignment of the
######                                          soft masked reads with novoalign
######							$1	=	/path/to/output directory
######							$2	=	sample name
######							$3	=	/path/to/run_info.txt
######		Output files:		BAM files
######          now checking the quality scores for fastq's and using the specific paramter for illumina and sanger quality in nova align
########################################################

if [ $# != 4 ];
then
    echo -e "Usage: wrapper script to run the alignment using NOVO ALIGN \n align_split_thread.sh <sample name> <output_dir> </path/to/run_info.txt>";
else	
    set -x 
    echo `date`
    sample=$1
    output_dir=$2
    read=$3
    run_info=$4

    
########################################################	
######	Reading run_info.txt and assigning to variables
    seq_file=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2)
    analysis=`echo "$analysis" | tr "[A-Z]" "[a-z]"`
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    center=$( cat $run_info | grep -w '^CENTER' | cut -d '=' -f2 )
    platform=$( cat $run_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
    fastqc=$( cat $tool_info | grep -w '^FASTQC' | cut -d '=' -f2)
    genome_bwa=$( cat $tool_info | grep -w '^BWA_REF' | cut -d '=' -f2)
    genome_novo=$( cat $tool_info | grep -w '^NOVO_REF' | cut -d '=' -f2)
    bwa=$( cat $tool_info | grep -w '^BWA' | cut -d '=' -f2)
    novoalign=$( cat $tool_info | grep -w '^NOVOALIGN' | cut -d '=' -f2)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
    filenames=$(cat $sample_info | grep -w "$sample" | cut -d '=' -f2| tr "\t" "\n")
    output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
    paired=$( cat $run_info | grep -w '^PAIRED' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)

########################################################	
######		Check FASTQ for Illumina or Sanger quality scrore
    
    output_dir_sample=$output_dir/alignment/$sample
    fastq=$output_dir/fastq
    fastqc=$output_dir/fastqc
    
    pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sample | cut -d ":" -f1)
    lane=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tr "," "\n" | head -n $SGE_TASK_ID | tail -n 1)
    index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tr "," "\n" | head -n $SGE_TASK_ID | tail -n 1)
    if [ $analysis == "mayo" ]
    then
        if [ $index == "-" ]
        then
            $java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -r $run_num -s Alignment -a WholeGenome -v $version
        else
            $java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -i $index -r $run_num -s Alignment -a WholeGenome -v $version
        fi    
    fi
    
    if [ $read -eq 1 ]
    then
	let fidx=($SGE_TASK_ID*2)-1 
    else
	let fidx=($SGE_TASK_ID*2)
    fi	
    R1=`cat $sample_info | grep -w "$sample" | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
    extension=$(echo $R1 | sed 's/.*\.//')
    filename1=$(echo $R1 | sed 's/\.[^\.]*$//')
    if [ $extension == "gz" ]
    then
        $script_path/fastq.sh $R1 $seq_file $filename1 $fastq $run_info $fastqc
    else
        filename1=$R1
        $script_path/fastq.sh $R1 $seq_file $filename1 $fastq $run_info $fastqc
    fi
## check if the fastqs are zipped or not
    
    if [ ! -s $fastq/$filename1 ]
    then
        echo "ERROR: $0 File $fastq/$filename1 does not exist" 
        exit 1
    fi
    ILL2SANGER1=`perl $script_path/checkFastqQualityScores.pl $fastq/$filename1 1000`
   

########################################################	
######		Run bwa alignemnt module
    echo `date`
    if [ $ILL2SANGER1 -gt 65 ]
    then
        $bwa/bwa aln -l 32 -t 4 -I $genome_bwa $fastq/$filename1 > $output_dir_sample/$sample.$SGE_TASK_ID.R$read.sai
    else
        $bwa/bwa aln -l 32 -t 4 $genome_bwa $fastq/$filename1 > $output_dir_sample/$sample.$SGE_TASK_ID.R$read.sai
    fi        
  
    
    if [ ! -s $output_dir_sample/$sample.$SGE_TASK_ID.R$read.sai ]
    then
        echo "ERROR: $0 File $output_dir_sample/$sample.$SGE_TASK_ID.sam not generated"
        exit 1 
    fi  
fi
	
   
