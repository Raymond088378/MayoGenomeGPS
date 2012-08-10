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

if [ $# -le 3 ];
then
    echo -e "Usage: SCRIPT to get Ontarget pileup \n<input dir> <pileup> <chromsome> <output Ontarget> <sample> <run ifno>";
else	
    set -x
    echo `date`
    input=$1
    output=$2
    sample=$3
    run_info=$4
    
    if [ $5 ]
    then
	SGE_TASK_ID=$5
    fi
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
    
    if [ ! -s $bam ]
    then
        $script_path/errorlog.sh $bam OnTarget_PILEUP.sh ERROR "not exist"
		exit 1;
    fi    
    #make bed format pileup
    
    if [ $multi == "YES" ]
    then
        pair=$( cat $sample_info | grep -w "^$sample" | cut -d '=' -f2 | tr "\t" " ")
        for i in $pair
        do
			$samtools/samtools index $input/$sample.$i.chr$chr.bam
			$samtools/samtools mpileup -A -s -f $ref $input/$sample.$i.chr$chr.bam | awk '{if ($1 ~ /chr/) {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}}'  > $output/$sample.$i.chr$chr.pileup.bed
            rm $input/$sample.$i.chr$chr.bam.bai
			total=`cat $output/$sample.$i.chr$chr.pileup.bed | wc -l`
            perl $script_path/split.a.file.into.n.parts.pl 25 $output/$sample.$i.chr$chr.pileup.bed $total
            rm $output/$sample.$i.chr$chr.pileup.bed
            cd $output
            # intersect the pileup with the intersect kit
            for ((j=1;j<=26; j++))
            do
                if [ -f $output/$sample.$i.chr$chr.pileup.bed.$j.txt ]
				then
					if [ -s $output/$sample.$i.chr$chr.pileup.bed.$j.txt ]
					then
						$bed/intersectBed -a $kit -b $output/$sample.$i.chr$chr.pileup.bed.$j.txt -wa -wb > $output/$sample.$i.chr$chr.pileup.bed.$j.txt.i
					fi
					rm $output/$sample.$i.chr$chr.pileup.bed.$j.txt	
				fi
			done
        
            #merge all the interscted pileup
            for((j=0; j<=99; j++))
            do
                total=0
                for((k=1; k<=26; k++))
				do
                    a=0
                    if [ -f $sample.$i.chr$chr.pileup.bed.$k.txt.i ]
					then
						a=`awk '$(NF-1)>'$j'' $sample.$i.chr$chr.pileup.bed.$k.txt.i | wc -l`
						total=`expr $total "+" $a`
					fi	
                done
                echo $total >> $output/$sample.$i.chr$chr.pileup.i.out
            done    
            for((k=1; k<=26; k++))
			do
				if [ -f $output/$sample.$i.chr$chr.pileup.bed.$k.txt.i ]
				then
					rm $output/$sample.$i.chr$chr.pileup.bed.$k.txt.i
				fi
			done	
        done    
    else	
		bam=$input/chr$chr.cleaned.bam
		$samtools/samtools mpileup -A -s -f $ref $bam | awk '{if ($1 ~ /chr/) {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}}'  > $output/$sample.chr$chr.pileup.bed
		#split the file into 25 parts to use less memory
        total=`cat $output/$sample.chr$chr.pileup.bed | wc -l`
        perl $script_path/split.a.file.into.n.parts.pl 25 $output/$sample.chr$chr.pileup.bed $total
        rm $output/$sample.chr$chr.pileup.bed
        cd $output
        # intersect the pileup with the intersect kit
        for ((i=1;i<=26; i++))
        do
			if [ -f $output/$sample.chr$chr.pileup.bed.$i.txt ]
			then
				if [ -s $output/$sample.chr$chr.pileup.bed.$i.txt ]
				then
					$bed/intersectBed -a $kit -b $output/$sample.chr$chr.pileup.bed.$i.txt -wa -wb > $output/$sample.chr$chr.pileup.bed.$i.txt.i
				fi
				rm $output/$sample.chr$chr.pileup.bed.$i.txt	
			fi	
        done
        
        #merge all the interscted pileup
        for((j=0; j<=99; j++))
        do
            total=0
            for ((i=1;i<=26; i++))
            do
                a=0
                if [ -f $sample.chr$chr.pileup.bed.$i.txt.i ]
				then
					a=`awk '$(NF-1)>'$j'' $sample.chr$chr.pileup.bed.$i.txt.i | wc -l`
					total=`expr $total "+" $a`
				fi
			done
            echo $total >> $output/$sample.chr$chr.pileup.i.out
        done    
        
		for((k=1; k<=26; k++))
		do
			if [ -f $output/$sample.chr$chr.pileup.bed.$k.txt.i ]
			then
				rm $output/$sample.chr$chr.pileup.bed.$k.txt.i
			fi
		done	
	fi
    echo `date`
fi	
    
	
	
	