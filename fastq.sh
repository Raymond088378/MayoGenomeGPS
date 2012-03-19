#!/bin/sh

if [ $# != 6 ]
then
    echo -e "Usage: SCRIPT unzip and run fastqc\n fastq.sh <input fastq> <input dir> <output fastq> <output dir> <run info> <FASTQC directory>"
else
    set -x
    echo `date`
    in_fastq=$1
    input=$2
    out_fastq=$3
    output=$4
    run_info=$5
    fastqc_dir=$6
    extension=$(echo $in_fastq | sed 's/.*\.//')
    filename=$(echo $in_fastq | sed 's/\.[^\.]*$//')
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    fastqc_path=$( cat $tool_info | grep -w '^FASTQC' | cut -d '=' -f2)
    FASTQC=$( cat $run_info | grep -w '^FASTQC' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    FOLDER_FASTQC=$( cat $run_info | grep -w '^FOLDER_FASTQC' | cut -d '=' -f2 )
    
    if [ $extension == "gz" ]
    then
        gunzip -c $input/$in_fastq > $output/${out_fastq}
    else
        ln -s $input/$in_fastq $output/${out_fastq}
    fi	
    
    if [ ! -s $output/${out_fastq} ]
    then
        echo "ERROR: fastq.sh $output/${out_fastq} is empty "
    fi
    
    ##FASTQC
    if [ $FASTQC == "YES" ]
    then
        $fastqc_path/fastqc -Xmx256m -Dfastqc.output_dir=$fastqc_dir/ $output/${out_fastq}
        ext=$(echo $out_fastq | sed 's/.*\.//')
		if [ $ext == "fastq" ]
		then
			file=$(echo $out_fastq | sed 's/\.[^\.]*$//')
		else
			file=$out_fastq
		fi	
    else
        read=${out_fastq}
		ext=$(echo $read | sed 's/.*\.//')
		if [ $ext == "fastq" ]
		then
			file=$(echo $read | sed 's/\.[^\.]*$//')
		else
			file=$read
		fi		
		ln -s $FOLDER_FASTQC/${file}_fastqc $fastqc_dir/
    fi		
    
    ## error checking
    if [ -s $fastqc_dir/${file}_fastqc ]
    then
        echo "FASTQC generation is successful "
    else
        echo "ERROR : FASTQC generation failed "
    fi	
    echo `date`
fi	