#!/bin/sh

########################################################
###### 	CV CALLER FOR TUMOR/NORMAL PAIR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			run_cnvnator.sh
######		Date:				09/26/2011
######		Summary:			Calls Crest
######		Input 
######		$1	=	group name
######		$2	=	bam list (first is normal)
######		$3	=	/path/to/output directory
######		$4	=	/path/to/run_info.txt
########################################################


if [ $# != 4 ]
then
	echo -e "Usage: \n groupname </path/to/input directory> </path/to/output directory> </path/to/run_info.txt>";
else
	set -x
	echo `date`
	group=$1	#from run ifo and sample info file (default=pair)
	input=$2	#colon separated list of full paths of bams to use
	output_dir=$3	#where you want it to go
	run_info=$4	#run_info
	## This script should be run by chromosome so submit as an array job -t 1-24
	## or else run all together and uncomment the following
#SGE_TASK_ID=2
	########################################################	
	######		Reading run_info.txt and assigning to variables

	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2)
	email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
	queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	segseq=$( cat $tool_info | grep -w '^SEGSEQ' | cut -d '=' -f2 )
	matlab=$( cat $tool_info | grep -w '^MATLAB' | cut -d '=' -f2 )
	chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
	ref_genome=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
	min_fold=$( cat $tool_info | grep -w '^MINFOLD' | cut -d '=' -f2 )
	max_fold=$( cat $tool_info | grep -w '^MAXFOLD' | cut -d '=' -f2 )
	gap=$( cat $tool_info | grep -w '^GAP_GENOME' | cut -d '=' -f2 )
	bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
	distgap=$( cat $tool_info | grep -w '^DISTGAP' | cut -d '=' -f2 )
	blacklist_sv=$( cat $tool_info | grep -w '^BLACKLIST_SV' | cut -d '=' -f2 )
	pct_overlap=$(cat $tool_info | grep -w '^STRUCT_PCT_BLACKLIST' | cut -d "=" -f2)
	PATH=$bedtools/:$PATH

	outdir=$output_dir/$group
	mkdir -p $outdir
	mkdir -p $outdir/info
	mkdir -p $outdir/map
    #input_bam=$input/chr${chr}.cleaned.bam
	let num_tumor=`echo $samples|tr " " "\n"|wc -l`-1
	normal_bam=`echo $samples| tr " " "\n" | head -n 1 `
	tumor_list=`echo $samples | tr " " "\n" | tail -$num_tumor`
	sample_normal=`echo $normal_bam`
	
	for tumor_bam in $tumor_list
	do
		sample_tumor=`basename $tumor_bam | sed -e "s/\.sorted\.bam//g"`

		#Create info and map files

		infofile=$outdir/info/$sample_tumor.$chr.info
		echo -e "File\tSample\tType" > $infofile
		#Create info and map files for normal
		# $samtools/samtools view -b -r $sample_normal $input_bam > $outdir/$sample_normal.bam
		ln -s $input/$group.$sample_normal.chr$chr.bam $outdir/$sample_normal.$chr.bam
		$samtools/samtools index $outdir/$sample_normal.$chr.bam
		for lane in `$samtools/samtools view $outdir/$sample_normal.$chr.bam | cut -f1 | sed -e "s/:/\t/g" | cut -f 1-2 | awk '{print $1":"$2}'| head -n 1000000 | sort | uniq | tail -n 2`
		do
			echo "Lane:$lane"
			echo -e "$outdir/map/$sample_normal.$lane.$chr.map\t$sample_normal\tNormal" >> $infofile
			$samtools/samtools view $outdir/$sample_normal.$chr.bam | awk "\$1~/$lane/" | awk '{print $3,"\t",$4,"\t",$9}'\
			|sed -e "s/chr//g" | sed -e "s/X/23/g" | sed -e "s/Y/24/g" | sed -e "s/M/25/g" \
			|awk 'BEGIN{OFS="\t"} {if ($3>0) print $1,$2,0; else print $1,$2,1}' > $outdir/map/$sample_normal.$lane.$chr.map
		done
        rm $outdir/$sample_normal.$chr.bam $outdir/$sample_normal.$chr.bam.bai
		# $samtools/samtools view -b -r $sample_tumor $input_bam > $outdir/$sample_tumor.bam
		ln -s $input/$group.$sample_tumor.chr$chr.bam $outdir/$sample_tumor.$chr.bam
		$samtools/samtools index $outdir/$sample_tumor.$chr.bam 
		#Create info and map files for tumor
		for lane in `$samtools/samtools view $outdir/$sample_tumor.$chr.bam | cut -f1 | sed -e "s/:/\t/g" | cut -f 1-2 | awk '{print $1":"$2}' | head -n 1000000 | sort | uniq | tail -n 2`
		do
			echo -e "$outdir/map/$sample_tumor.$lane.$chr.map\t$sample_tumor\tTumor" >> $infofile;

			$samtools/samtools view $outdir/$sample_tumor.$chr.bam | awk "\$1~/$lane/" |awk '{print $3,"\t",$4,"\t",$9}'\
			|sed -e "s/chr//g" | sed -e "s/X/23/g" | sed -e "s/Y/24/g" | sed -e "s/M/25/g" \
			|awk 'BEGIN{OFS="\t"} {if ($3>0) print $1,$2,0; else print $1,$2,1}' > $outdir/map/$sample_tumor.$lane.$chr.map
		done
		
        rm $outdir/$sample_tumor.$chr.bam $outdir/$sample_tumor.$chr.bam.bai
		export MATLABPATH=$MATLABPATH:$segseq
		cd $outdir
		$matlab/matlab -nodisplay -nojvm -nodesktop -nosplash -r "addpath $segseq;SegSeq -a 500 -i $infofile -s $sample_tumor.$chr;exit;"


		if [ ! -f $outdir/$sample_tumor.$chr.txt ]
		then
			echo "ERROR SegSeq.sh: File $outdir/$sample_tumor.$chr.*.txt not created"
			exit 1;
		fi

		outfile="$outdir/$sample_tumor.$chr.txt"

		cat $outfile | $script_path/filter_fold_segseq.pl $sample_tumor $min_fold > $outdir/$sample_tumor.$chr.del.bed
		cat $outfile | $script_path/filter_fold_segseq.pl $sample_tumor $max_fold > $outdir/$sample_tumor.$chr.dup.bed
		
		## convert to vcf
		perl $script_path/SegSeq2vcf.pl -i $outdir/$sample_tumor.$chr.del.bed -s $sample_tumor -o $outdir/$sample_tumor.$chr.del.vcf -r $ref_genome -t $samtools
		if [ ! -s $outdir/$sample_tumor.$chr.del.vcf.fail ]
		then
			rm $outdir/$sample_tumor.$chr.del.vcf.fail
		fi		
		perl $script_path/vcfsort.pl $ref_genome.fai $outdir/$sample_tumor.$chr.del.vcf > $outdir/$sample_tumor.$chr.del.vcf.sort
		mv $outdir/$sample_tumor.$chr.del.vcf.sort $outdir/$sample_tumor.$chr.del.vcf
		perl $script_path/SegSeq2vcf.pl -i $outdir/$sample_tumor.$chr.dup.bed -s $sample_tumor -o $outdir/$sample_tumor.$chr.dup.vcf -r $ref_genome -t $samtools
		if [ ! -s $outdir/$sample_tumor.$chr.dup.vcf.fail ]
		then
			rm $outdir/$sample_tumor.$chr.dup.vcf.fail
		fi
		perl $script_path/vcfsort.pl $ref_genome.fai $outdir/$sample_tumor.$chr.dup.vcf > $outdir/$sample_tumor.$chr.dup.vcf.sort
		mv $outdir/$sample_tumor.$chr.dup.vcf.sort $outdir/$sample_tumor.$chr.dup.vcf
		
		## filter the raw file
		
		$bedtools/closestBed -a $outdir/$sample_tumor.$chr.del.bed -b $gap -d | awk "\$13>$distgap" | cut -f 1-6 |\
		    $bedtools/intersectBed -a stdin -b $blacklist_sv -v -f $pct_overlap -wa > $outdir/$sample_tumor.$chr.filter.del.bed

		$bedtools/closestBed -a $outdir/$sample_tumor.$chr.dup.bed -b $gap -d | awk "\$13>$distgap" | cut -f 1-6 |\
		    $bedtools/intersectBed -a stdin -b $blacklist_sv -v -f $pct_overlap -wa > $outdir/$sample_tumor.$chr.filter.dup.bed

		### convert to vcf files
		perl $script_path/SegSeq2vcf.pl -i $outdir/$sample_tumor.$chr.filter.del.bed -s $sample_tumor -o $outdir/$sample_tumor.$chr.filter.del.vcf -r $ref_genome -t $samtools
		if [ ! -s $outdir/$sample_tumor.$chr.filter.del.vcf.fail ]
		then
			rm $outdir/$sample_tumor.$chr.filter.del.vcf.fail
		fi		
		perl $script_path/vcfsort.pl $ref_genome.fai $outdir/$sample_tumor.$chr.filter.del.vcf > $outdir/$sample_tumor.$chr.filter.del.vcf.sort
		mv $outdir/$sample_tumor.$chr.filter.del.vcf.sort $outdir/$sample_tumor.$chr.filter.del.vcf
		perl $script_path/SegSeq2vcf.pl -i $outdir/$sample_tumor.$chr.filter.dup.bed -s $sample_tumor -o $outdir/$sample_tumor.$chr.filter.dup.vcf -r $ref_genome -t $samtools
		if [ ! -s $outdir/$sample_tumor.$chr.filter.dup.vcf.fail ]
		then
			rm $outdir/$sample_tumor.$chr.filter.dup.vcf.fail
		fi
		perl $script_path/vcfsort.pl $ref_genome.fai $outdir/$sample_tumor.$chr.filter.dup.vcf > $outdir/$sample_tumor.$chr.filter.dup.vcf.sort
		mv $outdir/$sample_tumor.$chr.filter.dup.vcf.sort $outdir/$sample_tumor.$chr.filter.dup.vcf
		rm $outfile
		rm $outdir/$sample_tumor.$chr.mat
		rm $outdir/$sample_tumor.${chr}_aligned_reads.mat
	done
	echo `date`
fi

