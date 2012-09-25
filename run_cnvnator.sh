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

if [ $# -le 3 ]
then
    echo "\nUsage: samplename </path/to/output directory> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    sample=$1
    input=$2
    output_dir=$3
    run_info=$4
    if [ $5 ]
    then
    	SGE_TASK_ID=$5
    fi	

########################################################	
######		Reading run_info.txt and assigning to variables

    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    split_genome=$( cat $tool_info | grep -w '^SPLIT_GENOME' | cut -d '=' -f2)
    gap=$( cat $tool_info | grep -w '^GAP_GENOME' | cut -d '=' -f2 )
    blacklist_sv=$( cat $tool_info | grep -w '^BLACKLIST_SV' | cut -d '=' -f2 )
    pct_overlap=$(cat $tool_info | grep -w '^STRUCT_PCT_BLACKLIST' | cut -d "=" -f 2)
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    distgap=$( cat $tool_info | grep -w '^DISTGAP' | cut -d '=' -f2 )
    cnvnator=$( cat $tool_info | grep -w '^CNVNATOR' | cut -d '=' -f2 )
    rootlib=$( cat $tool_info | grep -w '^ROOTLIB' | cut -d '=' -f2 )
    bin_size=$( cat $tool_info | grep -w '^CNVNATOR_BINSIZE' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
########################################################	
######		

    input_bam=$input/chr${chr}.cleaned.bam
	
	$samtools/samtools view -H $input_bam 1>$input_bam.cnv.header 2>$input_bam.fix.cnv.log
	if [ `cat $input_bam.fix.cnv.log | wc -l` -gt 0 ]
	then
		$script_path/errorlog.sh $input_bam run_cnvnator.sh ERROR "truncated or corrupt bam"
		exit 1;
	else
		rm $input_bam.fix.cnv.log
	fi	
	rm $input_bam.cnv.header
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$rootlib
    root_file=$output_dir/$sample.$chr.root
    
	if [ -f $split_genome/chr$chr.fa ]
    then
        ln -s $split_genome/chr$chr.fa $output_dir/chr$chr.fa
    else
        $script_path/errorlog.sh $split_genome/chr$chr.fa run_cnvnator.sh ERROR "does not exist"
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
            $script_path/CNVnator2VCF.pl -i $output_dir/$sample.$chr.del.bed -f $ref -o $output_dir/$sample.$chr.del.vcf -s $sample -t $samtools
            if [ ! -s $output_dir/$sample.$chr.del.vcf.fail ]
            then
                rm $output_dir/$sample.$chr.del.vcf.fail
            fi    
            grep "^duplication" $output_dir/$sample.$chr.cnv.txt | $script_path/cnv2bed.pl $sample | sed -e "s/^/chr/" > $output_dir/$sample.$chr.dup.bed
            $script_path/CNVnator2VCF.pl -i $output_dir/$sample.$chr.dup.bed -f $ref -o $output_dir/$sample.$chr.dup.vcf -s $sample -t $samtools
            if [ ! -s $output_dir/$sample.$chr.dup.vcf.fail ]
            then
                rm $output_dir/$sample.$chr.dup.vcf.fail
            fi   

            if [ -s $output_dir/$sample.$chr.del.bed ]
			then
				$bedtools/closestBed -a $output_dir/$sample.$chr.del.bed -b $gap -d | awk "\$13>$distgap" | cut -f 1-6 |\
				$bedtools/intersectBed -a stdin -b $blacklist_sv -v -f $pct_overlap -wa  > $output_dir/$sample.$chr.filter.del.bed
			else
				touch $output_dir/$sample.$chr.filter.del.bed
			fi
			$script_path/CNVnator2VCF.pl -i $output_dir/$sample.$chr.filter.del.bed -f $ref -o $output_dir/$sample.$chr.filter.del.vcf -s $sample -t $samtools
            if [ ! -s $output_dir/$sample.$chr.filter.del.vcf.fail ]
            then
                rm $output_dir/$sample.$chr.filter.del.vcf.fail
            fi  
            if [ -s $output_dir/$sample.$chr.dup.bed ]
			then
				$bedtools/closestBed -a $output_dir/$sample.$chr.dup.bed -b $gap -d | awk "\$13>$distgap" | cut -f 1-6 |\
				$bedtools/intersectBed -a stdin -b $blacklist_sv -v -f $pct_overlap -wa > $output_dir/$sample.$chr.filter.dup.bed
			else
				touch $output_dir/$sample.$chr.filter.dup.bed
			fi
			$script_path/CNVnator2VCF.pl -i $output_dir/$sample.$chr.filter.dup.bed -f $ref -o $output_dir/$sample.$chr.filter.dup.vcf -s $sample -t $samtools
            if [ ! -s $output_dir/$sample.$chr.filter.dup.vcf.fail ]
            then
                rm $output_dir/$sample.$chr.filter.dup.vcf.fail
            fi
            rm $root_file $output_dir/chr$chr.fa $output_dir/$sample.$chr.cnv.txt
        else
           $script_path/errorlog.sh $output_dir/$sample.$chr.cnv.txt run_cnvnator.sh ERROR "does not exist"
           exit 1;
        fi
    else
		$script_path/errorlog.sh $root_file run_cnvnator.sh ERROR "does not exist"
        exit 1;
    fi
    echo `date`
fi
