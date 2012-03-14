#!/bin/sh
##	INFO
#	To Intersect pileup with OnTarget Kit by splitting the bam file into 200 files

######################################
#		$1		=	input folder (realignment sample folder)
#		$2		=	chromsome index
#		$3		=	Ontarget output folder
#		$4		=	sample name
#		$5		=	run info file
#########################################

if [ $# != 4 ];
then
    echo -e "Usage: SCRIPT to get Ontarget pileup \n<input dir> <pileup> <chromsome> <output Ontarget> <sample> <run ifno>";
else	
    set -x
    echo `date`
    input=$1
    output=$2
    sample=$3
    run_info=$4
    #SGE_TASK_ID=10
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    bed=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    CaptureKit=$( cat $tool_info | grep -w '^CAPTUREKIT' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    tool=`echo "$tool" | tr "[A-Z]" "[a-z]"`
    out=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    out_dir=$out/$PI/$tool/$run_num
    multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2)
    multi=`echo $multi | tr "[a-z]" "[A-Z]"`
    PATH=$bed/:$PATH
    #cd $output
    if [ $tool == "whole_genome" ]
    then
        kit=$out_dir/bed_file.bed
    else
        kit=$CaptureKit
    fi
    
    bam=$input/chr$chr.cleaned.bam
    if [ ! -s $bam ]
    then
        echo "ERROR : OnTarget.PILEUP.sh $bam not found"
    fi    
    #make bed format pileup
    
    if [ $multi == "YES" ]
    then
        pair=$( cat $sample_info | grep -w "$sample" | cut -d '=' -f2)
        for i in $pair
        do
            $samtools/samtools view -b -r $i $bam | $samtools/samtools pileup -s -f $ref - | awk '{if ($1 ~ /chr/) {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}}'  > $output/$i.chr$chr.pileup.bed
            total=`cat $output/$i.chr$chr.pileup.bed | wc -l`
            perl $script_path/split.a.file.into.n.parts.pl 25 $output/$i.chr$chr.pileup.bed $total
            rm $output/$i.chr$chr.pileup.bed
            cd $output
            # intersect the pileup with the intersect kit
            for j in $i.chr$chr.pileup.bed.*.txt
            do
                $bed/intersectBed -a $kit -b $output/$j -wa -wb > $output/$j.i
                rm $output/$j	
            done
        
            #merge all the interscted pileup
            for((j=0; j<=39; j++))
            do
                total=0
                for k in $i.chr$chr.pileup.bed.*.txt.i
                do
                    a=0
                    a=`awk '$(NF-1)>'$j'' $k | wc -l`
                    total=`expr $total "+" $a`
                done
                echo $total >> $output/$i.chr$chr.pileup.i.out
            done    
            rm $output/$i.chr$chr.pileup.bed.*.txt.i
        done    
    else
        if [ -s $input/chr$chr.pileup ]
		then
			cat $input/chr$chr.pileup | awk '{if ($1 ~ /chr/) {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}}'  > $output/$sample.chr$chr.pileup.bed
		else	
			$samtools/samtools pileup -s -f $ref $bam | awk '{if ($1 ~ /chr/) {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}}'  > $output/$sample.chr$chr.pileup.bed
        fi
		#split the file into 25 parts to use less memory
        total=`cat $output/$sample.chr$chr.pileup.bed | wc -l`
        perl $script_path/split.a.file.into.n.parts.pl 25 $output/$sample.chr$chr.pileup.bed $total
        rm $output/$sample.chr$chr.pileup.bed
        cd $output
        # intersect the pileup with the intersect kit
        for i in $sample.chr$chr.pileup.bed.*.txt
        do
                $bed/intersectBed -a $kit -b $output/$i -wa -wb > $output/$i.i
                rm $output/$i	
        done
        
        #merge all the interscted pileup
        for((j=0; j<=39; j++))
        do
            total=0
            for i in $sample.chr$chr.pileup.bed.*.txt.i
            do
                a=0
                a=`awk '$(NF-1)>'$j'' $i | wc -l`
                total=`expr $total "+" $a`
            done
            echo $total >> $output/$sample.chr$chr.pileup.i.out
        done    
        rm $output/$sample.chr$chr.pileup.bed.*.txt.i
    fi
    echo `date`
fi	
    
	
	
	