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

if [ $# != 5 ]
then
    echo -e "Usage: Script to run crest on a paired sample \n <sample name> <group name> </path/to/input directory> </path/to/output directory> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    sample=$1
    group=$2
    input=$3
    output_dir=$4
    run_info=$5



    ########################################################	
    ######		Reading run_info.txt and assigning to variables
    #SGE_TASK_ID=1
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
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
    samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2)
    #bam=$input/chr$chr.cleaned.bam

    min_read=$( cat $tool_info | grep -w '^STRUCT_MIN_SUPPORT' | cut -d '=' -f2)
    min_id=$( cat $tool_info | grep -w '^STRUCT_MIN_IDENTITY' | cut -d '=' -f2)
    blacklist_sv=$( cat $tool_info | grep -w '^BLACKLIST_SV' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    PATH=$bedtools/:$PATH

    echo `date`
    PERL5LIB=$perllib
    PATH=$PATH:$blat:$crest
    mkdir -p $output_dir/$group

    mkdir -p $output_dir/$group/log




    ln -s $input/$group.$sample.chr$chr.bam $output_dir/$group/$sample.chr$chr.bam
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
        $samtools/samtools index $file
    fi

    $crest/extractSClip.pl -i $file -r chr$chr --ref_genome $ref_genome -o $output_dir/$group -p $sample
    echo `date`
fi
