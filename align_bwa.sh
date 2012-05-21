if [ $# != 3 ];
then
    echo -e "Usage: wrapper script to run the alignment using NOVO ALIGN \n align_split_thread.sh <sample name> <output_dir> </path/to/run_info.txt>";
else	
    set -x 
    echo `date`
    sample=$1
    output_dir=$2
    run_info=$3
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	genome_bwa=$( cat $tool_info | grep -w '^BWA_REF' | cut -d '=' -f2)
	bwa=$( cat $tool_info | grep -w '^BWA' | cut -d '=' -f2)
	center=$( cat $run_info | grep -w '^CENTER' | cut -d '=' -f2 )
	platform=$( cat $run_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
	GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	output_dir_sample=$output_dir/alignment/$sample
    fastq=$output_dir/fastq
	paired=$( cat $run_info | grep -w '^PAIRED' | cut -d '=' -f2)
	
	if [ $paired == 1 ]
	then
		let fidx=($SGE_TASK_ID*2)-1 
		let sidx=($SGE_TASK_ID*2)
		R1=`cat $sample_info | grep -w "$sample" | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
        R2=`cat $sample_info | grep -w "$sample" | cut -d '=' -f2| tr "\t" "\n" | head -n $sidx | tail -n 1`
		j=1
        for i in $R1 $R2
        do
            extension=$(echo $i | sed 's/.*\.//')
			if [ $extension == "gz" ]
            then
                eval filename${j}=$(echo $i | sed 's/\.[^\.]*$//')
            else
                eval filename${j}=$i
            fi
            let j=j+1
        done
		
		$bwa/bwa sampe -r "@RG\tID:$sample\tSM:$sample\tLB:$GenomeBuild\tPL:$platform\tCN:$center" $genome_bwa $output_dir_sample/$sample.$SGE_TASK_ID.R1.sai $output_dir_sample/$sample.$SGE_TASK_ID.R2.sai $fastq/$filename1 $fastq/$filename2 > $output_dir_sample/$sample.$SGE_TASK_ID.sam 
	else
		let fidx=($SGE_TASK_ID*2)-1 
		R1=`cat $sample_info | grep -w "$sample" | cut -d '=' -f2| tr "\t" "\n" | head -n $fidx | tail -n 1`
		extension=$(echo $R1 | sed 's/.*\.//')
		if [ $extension == "gz" ]
		then
			filename1=$(echo $R1 | sed 's/\.[^\.]*$//')
		else
			filename1=$R1
		fi	
		$bwa/bwa samse -r "@RG\tID:$sample\tSM:$sample\tLB:$GenomeBuild\tPL:$platform\tCN:$center" $genome_bwa $output_dir_sample/$sample.$SGE_TASK_ID.R1.sai $fastq/$filename1 > $output_dir_sample/$sample.$SGE_TASK_ID.sam 	
	fi
	
	if [ ! -s $output_dir_sample/$sample.$SGE_TASK_ID.sam ]
    then
        $script_path/errorlog.sh align_bwa.sh $output_dir_sample/$sample.$SGE_TASK_ID.sam ERROR empty
        exit 1
    else
        if [ $paired == 0 ]
        then
            rm $fastq/$filename1
            rm $output_dir_sample/$sample.$SGE_TASK_ID.R1.sai
        else
            rm $fastq/$filename1 $fastq/$filename2
            rm $output_dir_sample/$sample.$SGE_TASK_ID.R1.sai $output_dir_sample/$sample.$SGE_TASK_ID.R2.sai
        fi    
    fi  

    $samtools/samtools view -bt $ref.fai $output_dir_sample/$sample.$SGE_TASK_ID.sam > $output_dir_sample/$sample.$SGE_TASK_ID.bam  
    if [ ! -s $output_dir_sample/$sample.$SGE_TASK_ID.bam ]
    then
        $script_path/errorlog.sh align_bwa.sh $output_dir_sample/$sample.$SGE_TASK_ID.bam ERROR empty
        exit 1
    else
        rm $output_dir_sample/$sample.$SGE_TASK_ID.sam  
    fi
	
	########################################################	
######		Sort BAM, adds RG 

    $script_path/convert.bam.sh $output_dir_sample $sample.$SGE_TASK_ID.bam $sample.$SGE_TASK_ID $SGE_TASK_ID $run_info
    echo `date`
fi	