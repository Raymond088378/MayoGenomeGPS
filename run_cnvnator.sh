#!/bin/sh

########################################################
###### 	CNV CALLER FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			run_cnvnator.sh
######		Date:				09/26/2011
######		Summary:			Calls CNVnator
######		Input 
######		$1	=	samplename
######		$2	=	input_bam
######		$3	=	/path/to/output directory
######		$4	=	/path/to/run_info.txt
######		Output files:	BAM files. 
########################################################

if [ $# != 4 ]
then
    echo "\nUsage: samplename </path/to/output directory> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    sample=$1
    input=$2
    output_dir=$3
    run_info=$4
    #SGE_TASK_ID=1

########################################################	
######		Reading run_info.txt and assigning to variables

    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    split_genome=$( cat $tool_info | grep -w '^SPLIT_GENOME' | cut -d '=' -f2)
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
    lqueue=$( cat $run_info | grep -w '^LQUEUE' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    gap=$( cat $tool_info | grep -w '^GAP_GENOME' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    distgap=$( cat $tool_info | grep -w '^DISTGAP' | cut -d '=' -f2 )
    cnvnator=$( cat $tool_info | grep -w '^CNVNATOR' | cut -d '=' -f2 )
    rootlib=$( cat $tool_info | grep -w '^ROOTLIB' | cut -d '=' -f2 )
    bin_size=$( cat $tool_info | grep -w '^CNVNATOR_BINSIZE' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)

    output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	PATH=$bedtools/:$PATH
########################################################	
######		

    input_bam=$input/chr${chr}.cleaned.bam
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$rootlib
    root_file=$output_dir/$sample.$chr.root
    
	if [ -f $split_genome/chr$chr.fa ]
    then
        ln -s $split_genome/chr$chr.fa $output_dir/chr$chr.fa
    else
        echo "ERROR : Cnvnator: Chromosome reference file: $split_genome/chr$chr.fa does not exist"
        exit 1;
    fi
    cd $output_dir

    $cnvnator/cnvnator -root $root_file -chrom $chr -tree $input_bam

    if [ -s $rootfile ]
    then
        $cnvnator/cnvnator -root $root_file -chrom $chr -his $bin_size
        $cnvnator/cnvnator -root $root_file -chrom $chr -stat $bin_size
        $cnvnator/cnvnator -root $root_file -chrom $chr -partition $bin_size
        $cnvnator/cnvnator -root $root_file -chrom $chr -call $bin_size > $output_dir/$sample.$chr.cnv.txt

        if [ -f $output_dir/$sample.$chr.cnv.txt ]
        then
            grep "^deletion" $output_dir/$sample.$chr.cnv.txt | $script_path/cnv2bed.pl $sample | sed -e "s/^/chr/" > $output_dir/$sample.$chr.del.bed
            perl $script_path/CNVnator2VCF.pl -i $output_dir/$sample.$chr.del.bed -f $ref -o $output_dir/$sample.$chr.del.vcf -s $sample -t $samtools
            if [ ! -s $output_dir/$sample.$chr.del.vcf.fail ]
            then
                rm $output_dir/$sample.$chr.del.vcf.fail
            fi    
            grep "^duplication" $output_dir/$sample.$chr.cnv.txt | $script_path/cnv2bed.pl $sample | sed -e "s/^/chr/" > $output_dir/$sample.$chr.dup.bed
            perl $script_path/CNVnator2VCF.pl -i $output_dir/$sample.$chr.dup.bed -f $ref -o $output_dir/$sample.$chr.dup.vcf -s $sample -t $samtools
            if [ ! -s $output_dir/$sample.$chr.dup.vcf.fail ]
            then
                rm $output_dir/$sample.$chr.dup.vcf.fail
            fi   
            $bedtools/closestBed -a $output_dir/$sample.$chr.del.bed -b $gap -d | awk "\$13>$distgap" | cut -f 1-6 > $output_dir/$sample.$chr.filter.del.bed
            perl $script_path/CNVnator2VCF.pl -i $output_dir/$sample.$chr.filter.del.bed -f $ref -o $output_dir/$sample.$chr.filter.del.vcf -s $sample -t $samtools
            if [ ! -s $output_dir/$sample.$chr.filter.del.vcf.fail ]
            then
                rm $output_dir/$sample.$chr.filter.del.vcf.fail
            fi  
            $bedtools/closestBed -a $output_dir/$sample.$chr.dup.bed -b $gap -d | awk "\$13>$distgap" | cut -f 1-6 > $output_dir/$sample.$chr.filter.dup.bed
            perl $script_path/CNVnator2VCF.pl -i $output_dir/$sample.$chr.filter.dup.bed -f $ref -o $output_dir/$sample.$chr.filter.dup.vcf -s $sample -t $samtools
            if [ ! -s $output_dir/$sample.$chr.filter.dup.vcf.fail ]
            then
                rm $output_dir/$sample.$chr.filter.dup.vcf.fail
            fi
            rm $root_file $output_dir/chr$chr.fa $output_dir/$sample.$chr.cnv.txt
        else
           echo "ERROR : Cnvnator: $output_dir/$sample.$chr.cnv.txt not"
           exit 1;
        fi
    else
        echo "ERROR : Cnvnator: rootfile not properly created: $root_file"
        exit 1;
    fi
    echo `date`
fi
