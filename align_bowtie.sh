#!/bin/sh



### DO NOT USE 


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

if [ $# != 3 ];
then
    echo -e "Usage: wrapper script to run the alignment using NOVO ALIGN \n align_split_thread.sh <sample name> <output_dir> </path/to/run_info.txt>";
else	
    set -x 
    echo `date`
    sample=$1
    output_dir=$2
    run_info=$3

    
########################################################	
######	Reading run_info.txt and assigning to variables
    seq_file=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    center=$( cat $run_info | grep -w '^CENTER' | cut -d '=' -f2 )
    platform=$( cat $run_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
    fastqc=$( cat $tool_info | grep -w '^FASTQC' | cut -d '=' -f2)
    readlen=$( cat $run_info | grep -w '^READLENGTH' | cut -d '=' -f2)
    genome_bowtie=$( cat $tool_info | grep -w '^BOWTIE_REF' | cut -d '=' -f2)
 
    bowtie=$( cat $tool_info | grep -w '^BOWTIE' | cut -d '=' -f2)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    filenames=$(cat $sample_info | grep -w "$sample" | cut -d '=' -f2| tr "\t" "\n")
    output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
	threads=$( cat $tool_info | grep -w '^THREADS'| cut -d '=' -f2)
    paired=$( cat $run_info | grep -w '^PAIRED' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)

########################################################	
######		Check FASTQ for Illumina or Sanger quality scrore
    
    output_dir_sample=$output_dir/alignment/$sample
    fastq=$output_dir/fastq
    fastqc=$output_dir/fastqc
    
	$script_path/dashboard.sh $sample $run_info Alignment started $SGE_TASK_ID
   
    if [ $paired == 1 ]
    then
        let fidx=($SGE_TASK_ID*2)-1 
        let sidx=$SGE_TASK_ID*2
        R1=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
        R2=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $sidx | tail -n 1`
		## run fastqc depending on flag and convert the zip to unzip fastq
        j=1
        for i in $R1 $R2
        do
            extension=$(echo $i | sed 's/.*\.//')
            if [ $extension == "gz" ]
            then
                eval filename${j}=$(echo $i | sed 's/\.[^\.]*$//')
                $script_path/fastq.sh $i $seq_file $(eval echo \$filename$j) $fastq $run_info $fastqc
            else
                eval filename${j}=$i
                $script_path/fastq.sh $i $seq_file $(eval echo \$filename$j) $fastq $run_info $fastqc
            fi
            let j=j+1
        done
    elif [ $paired == 0 ]
    then
        let fidx=$SGE_TASK_ID
        R1=`cat $sample_info | grep -w ^FASTQ:$sample | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
        extension=$(echo $R1 | sed 's/.*\.//')
        filename1=$(echo $R1 | sed 's/\.[^\.]*$//')
        if [ $extension == "gz" ]
        then
            $script_path/fastq.sh $R1 $seq_file $filename1 $fastq $run_info $fastqc
        else
            filename1=$R1
            $script_path/fastq.sh $R1 $seq_file $filename1 $fastq $run_info $fastqc
        fi
    else
        echo "ERROR: $run info file is not synced with $sample_info file (sequencing information is not correct)"
        exit 1;
    fi    
    ## check if the fastqs are zipped or not
    
    if [ ! -s $fastq/$filename1 ]
    then
        echo "ERROR: align_novo.sh File $fastq/$filename1 does not exist" 
        exit 1
    fi

    if [ $paired == 1 ]
    then
        if [ ! -s $fastq/$filename2 ]
        then
            echo "ERROR: align_novo.sh File $fastq/$filename2 does not exist"
            exit 1
        fi
    fi    
    ILL2SANGER1=`perl $script_path/checkFastqQualityScores.pl $fastq/$filename1 1000`
    if [ $paired == 1 ]
    then
        ILL2SANGER2=`perl $script_path/checkFastqQualityScores.pl $fastq/$filename2 1000`
    fi

########################################################	
######		Run novoalign for PE or SR
    echo `date`
    if [ $paired == 1 ]
    then
        if [ $ILL2SANGER1 -gt 65 ] && [ $ILL2SANGER2 -gt 65 ]
        then
            $novoalign --hdrhd off -v 120 -c $threads -i PE 425,80 -x 5 -r Random -d $genome_novo -F ILMFQ -f $fastq/$filename1 $fastq/$filename2 \
            -o SAM "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:$platform\tCN:$center" > $output_dir_sample/$sample.$SGE_TASK_ID.sam
        
        else
            $novoalign --hdrhd off -v 120 -c $threads -i PE 425,80 -x 5 -r Random -d $genome_novo -F STDFQ -f $fastq/$filename1 $fastq/$filename2 \
            -o SAM "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:$platform\tCN:$center" > $output_dir_sample/$sample.$SGE_TASK_ID.sam   
        fi
    else
        if [ $ILL2SANGER1 -gt 65 ]
        then
            $novoalign --hdrhd off -v 120 -c $threads -i PE 425,80 -x 5 -r Random -d $genome_novo -F ILMFQ -f $fastq/$filename1 \
            -o SAM "@RG\tID:$sample\tSM:$sample\tLB:$GenomeBuild\tPL:$platform\tCN:$center" > $output_dir_sample/$sample.$SGE_TASK_ID.sam
        else
            $novoalign --hdrhd off -v 120 -c $threads -i PE 425,80 -x 5 -r Random -d $genome_novo -F STDFQ -f $fastq/$filename1 \
            -o SAM "@RG\tID:$sample\tSM:$sample\tLB:$GenomeBuild\tPL:$platform\tCN:$center" > $output_dir_sample/$sample.$SGE_TASK_ID.sam   
        fi        
    fi    
    
    if [ ! -s $output_dir_sample/$sample.$SGE_TASK_ID.sam ]
    then
        echo "ERROR: align_novo.sh File $output_dir_sample/$sample.$SGE_TASK_ID.sam not generated"
        exit 1
    else
        if [ $paired == 0 ]
        then
            rm $fastq/$filename1
        else
            rm $fastq/$filename1 $fastq/$filename2
        fi    
    fi  

    $samtools/samtools view -bS $output_dir_sample/$sample.$SGE_TASK_ID.sam > $output_dir_sample/$sample.$SGE_TASK_ID.bam  
    if [ ! -s $output_dir_sample/$sample.$SGE_TASK_ID.bam ]
    then
        echo "Error: align_novo.sh File $output_dir_sample/$sample.$SGE_TASK_ID.bam not generated"
        exit 1
    else
        rm $output_dir_sample/$sample.$SGE_TASK_ID.sam  
    fi


########################################################	
######		Sort BAM, adds RG & remove duplicates

    $script_path/convert.bam.sh $output_dir_sample $sample.$SGE_TASK_ID.bam $sample.$SGE_TASK_ID $SGE_TASK_ID $run_info
fi