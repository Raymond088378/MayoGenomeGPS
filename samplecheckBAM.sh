#!/bin/bash

if [ $# != 7 ]
then
    echo -e "wrapper to validate bam file\nUsage: samplecheckBAM.sh <input dir><input bam><outputdir><run_info><sample><1 or 0 if bam is per chr><which _chr>";
else	
    set -x
    echo `date`
    input=$1
    bam=$2
    output=$3
    run_info=$4
    sample=$5
    chopped=$6
    chr=$7
    
    if [ -d $output/temp ]
    then
        echo "already there"
    else
        mkdir -p $output/temp
    fi
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)	
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    
    $samtools/samtools view -H $input/$bam 1>$input/$bam.sc.$chr.header 2>$input/$bam.fix.sc.$chr.log
	
	if [ `cat $input/$bam.fix.sc.$chr.log | wc -l` -gt 0 ]
	then
		$script_path/errorlog.sh $input/$bam samplecheckBAM.sh ERROR "is truncated or corrupted"
		exit 1;
	else
		rm $input/$bam.fix.sc.$chr.log
	fi
	rm $input/$bam.sc.$chr.header	
	
	check=`[ -f $input/$bam.bai ] && echo "1" || echo "0"`
    if [ $check -eq 0 ]
    then
		ln -s $input/$bam $output/$bam.$chr.bam
        $samtools/samtools index $output/$bam.$chr.bam
    else
        ln -s $input/$bam $output/$bam.$chr.bam
        ln -s $input/$bam.bai $output/$bam.$chr.bam.bai	
    fi	
            
    ## extracting the BAM for specific chromosome
    if [ $chopped == 0 ]
    then
        $samtools/samtools view -b $output/$bam.$chr.bam chr${chr} > $output/$sample.chr${chr}.bam
        $samtools/samtools index $output/$sample.chr${chr}.bam
    else
        ln -s  $output/$bam.$chr.bam $output/$sample.chr${chr}.bam
        $samtools/samtools index $output/$sample.chr${chr}.bam
    fi
    
    ## check if BAM is sorted
    SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $output/$sample.chr${chr}.bam -s $samtools`
    if [ $SORT_FLAG == 1 ]
    then
        ln -s $output/$sample.chr${chr}.bam $output/$sample.chr${chr}-sorted.bam
    else
        $script_path/sortbam.sh $output/$sample.chr${chr}.bam $output/$sample.chr${chr}-sorted.bam $output/temp/ coordinate true $run_info
    fi
    
    ## check if read group and platform is availbale in BAM
    RG_FLAG=`perl $script_path/checkBAMreadGroup.pl -i $output/$sample.chr${chr}-sorted.bam -s $samtools`
    if [ $RG_FLAG == 0 ]
    then
        $script_path/addReadGroup.sh $output/$sample.chr${chr}-sorted.bam $output/$sample.chr${chr}-sorted.bam.rg.bam $output/temp/ $run_info $sample   
    else
        ### after this point BAM is good to go with the GATk tool (any module)
		$samtools/samtools index $output/$sample.chr${chr}-sorted.bam    
    fi
    echo `date`
fi	
	